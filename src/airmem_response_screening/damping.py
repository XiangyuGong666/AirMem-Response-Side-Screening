"""EFDD and RDT estimators operating on one measured acceleration window.

These estimators do not select groups, load archives or access result files.
"""
from __future__ import annotations

import math
from typing import Any, Dict, Mapping, Tuple

import numpy as np
from scipy import signal



# ---------- EFDD 辅助 ----------

def _complex_mac(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.vdot(a, a).real * np.vdot(b, b).real)
    return float(abs(np.vdot(a, b)) ** 2 / denom) if denom > np.finfo(float).eps else 0.0


def _leading_csd_curves(
    data: np.ndarray,
    fs: float,
    nperseg: int,
    nfft: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """单段 Welch CSD 主导奇异曲线。"""
    if data.shape[1] != nperseg:
        raise ValueError(f"时窗样本数 {data.shape[1]} 与 Welch nperseg={nperseg} 不一致")
    detrended = signal.detrend(data, axis=1, type="constant")
    window = signal.get_window("hann", nperseg, fftbins=True)
    fft = np.fft.rfft(detrended * window[None, :], n=nfft, axis=1)
    scale = 1.0 / (fs * float(np.sum(window**2)))
    freqs = np.fft.rfftfreq(nfft, 1.0 / fs)
    power = np.sum(np.abs(fft) ** 2, axis=0)
    one_sided = np.ones(freqs.size, dtype=float)
    if freqs.size > 2:
        one_sided[1:-1 if nfft % 2 == 0 else None] = 2.0
    s1 = power * scale * one_sided
    norms = np.sqrt(np.maximum(power, np.finfo(float).eps))
    vectors = (np.conjugate(fft) / norms[None, :]).T
    return freqs, s1, vectors


def _efdd_failure(reason: str, f_refined: float = float("nan"), line_count: int = 0) -> Dict[str, Any]:
    return {
        "f_refined_Hz": f_refined,
        "DampingRatio": float("nan"),
        "Fit_R2": float("nan"),
        "BellLineCount": line_count,
        "PassedQualityGate": False,
        "FailureReason": reason,
    }


def estimate_efdd_window(
    data: np.ndarray,
    f0: float,
    group_low: float,
    group_high: float,
    fs: float,
    cfg: Mapping[str, Any],
) -> Dict[str, Any]:
    """对单窗数据做 EFDD 阻尼比估计。"""
    window_sec = float(cfg.get("window_sec", 60.0))
    nperseg = int(cfg.get("nperseg") or round(window_sec * fs))
    nfft = nperseg * int(cfg["zero_padding_factor"])
    freqs, s1, vecs = _leading_csd_curves(data, fs, nperseg, nfft)
    fraction = float(cfg["local_band_fraction"])
    local = np.flatnonzero(
        (freqs >= max(0.0, f0 * (1.0 - fraction))) & (freqs <= f0 * (1.0 + fraction))
    )
    if local.size == 0:
        return _efdd_failure("局部频带无频率线")

    occurrence_band = local[(freqs[local] >= group_low) & (freqs[local] <= group_high)]
    reference_candidates = occurrence_band if occurrence_band.size else local
    ref = int(reference_candidates[np.argmax(s1[reference_candidates])])
    ref_vec = vecs[ref]
    macs = np.asarray([_complex_mac(vecs[k], ref_vec) for k in local])
    accepted = macs >= float(cfg["bell_mac_threshold"])
    position = int(np.where(local == ref)[0][0])
    left = position
    right = position
    while left > 0 and accepted[left - 1]:
        left -= 1
    while right + 1 < local.size and accepted[right + 1]:
        right += 1
    bell_idx = local[left : right + 1]
    line_count = int(bell_idx.size)
    if line_count < int(cfg["minimum_bell_lines"]):
        return _efdd_failure("钟形区频率线不足", line_count=line_count)

    weights = s1[bell_idx]
    if not np.isfinite(weights).all() or float(weights.sum()) <= 0:
        return _efdd_failure("钟形区第一奇异值无效", line_count=line_count)
    f_refined = float(np.sum(freqs[bell_idx] * weights) / np.sum(weights))
    bell = np.zeros_like(s1)
    bell[bell_idx] = weights
    autocorr = np.fft.irfft(bell, n=nfft)
    autocorr /= max(abs(float(autocorr[0])), np.finfo(float).eps)
    envelope = np.abs(signal.hilbert(autocorr))

    fit_cfg = cfg["fit"]
    times = np.arange(nfft, dtype=float) / fs
    mask = (
        (times >= float(fit_cfg["start_time_s"]))
        & (times <= min(float(fit_cfg["maximum_time_s"]), times[-1]))
        & (envelope >= float(fit_cfg["minimum_envelope_fraction"]) * envelope[0])
        & np.isfinite(envelope)
        & (envelope > 0)
    )
    fit_t = times[mask]
    fit_y = np.log(envelope[mask])
    if fit_t.size < int(fit_cfg["minimum_points"]):
        return _efdd_failure("衰减拟合点不足", f_refined, line_count)

    slope, intercept = np.polyfit(fit_t, fit_y, 1)
    predicted = slope * fit_t + intercept
    ss_res = float(np.sum((fit_y - predicted) ** 2))
    ss_tot = float(np.sum((fit_y - fit_y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > np.finfo(float).eps else 0.0
    alpha = -float(slope)
    zeta = (
        alpha / math.sqrt((2.0 * math.pi * f_refined) ** 2 + alpha**2)
        if alpha > 0
        else float("nan")
    )
    reasons = []
    if not np.isfinite(zeta):
        reasons.append("包络未衰减")
    if r2 < float(fit_cfg["r2_min"]):
        reasons.append("拟合R²低于门限")
    if np.isfinite(zeta) and not (
        float(cfg["quality"]["zeta_min"]) <= zeta <= float(cfg["quality"]["zeta_max"])
    ):
        reasons.append("阻尼比超出合理区间")

    return {
        "f_refined_Hz": f_refined,
        "DampingRatio": zeta,
        "Fit_R2": r2,
        "BellLineCount": line_count,
        "PassedQualityGate": not reasons,
        "FailureReason": ";".join(reasons),
    }


# ---------- RDT 辅助 ----------

def _phase_aligned_real_weights(mode_vector: np.ndarray) -> np.ndarray:
    """将复模态向量相位对齐并转成实模态坐标权重。"""
    vector = np.asarray(mode_vector, dtype=np.complex128).reshape(-1)
    if vector.size != 12 or not np.isfinite(vector.real).all() or not np.isfinite(vector.imag).all():
        raise ValueError("实测主响应向量无效")
    dominant = int(np.argmax(np.abs(vector)))
    if abs(vector[dominant]) <= np.finfo(float).eps:
        raise ValueError("实测主响应向量为零")
    aligned = vector * np.exp(-1j * np.angle(vector[dominant]))
    weights = aligned.real
    norm = float(np.linalg.norm(weights))
    if norm <= np.finfo(float).eps:
        raise ValueError("相位对齐后的实测主响应向量为零")
    return weights / norm


def _rdt_failure(
    reason: str,
    possible_crosstalk: bool,
    trigger_count: int = 0,
    fit_peak_count: int = 0,
) -> Dict[str, Any]:
    return {
        "f_refined_Hz": float("nan"),
        "DampingRatio": float("nan"),
        "Fit_R2": float("nan"),
        "TriggerCount": trigger_count,
        "FitPeakCount": fit_peak_count,
        "PassedQualityGate": False,
        "FailureReason": reason,
        "PossibleCrosstalk": bool(possible_crosstalk),
        "BeatDistortion": False,
        "PeakSpacingCV": float("nan"),
        "MaxEnvelopeIncrease": float("nan"),
    }


def estimate_rdt_window(
    data: np.ndarray,
    f0: float,
    mode_vector: np.ndarray,
    fs: float,
    cfg: Mapping[str, Any],
    possible_crosstalk: bool = False,
    return_curve: bool = False,
) -> Dict[str, Any]:
    """对单窗数据做 RDT 阻尼比估计。"""
    rdt = cfg["rdt"]
    try:
        weights = _phase_aligned_real_weights(mode_vector)
    except Exception as exc:
        return _rdt_failure(str(exc), possible_crosstalk)

    matrix = np.asarray(data, dtype=float)
    if matrix.shape[0] != 12 or matrix.shape[1] < 20:
        return _rdt_failure("时窗数据不足", possible_crosstalk)

    matrix = signal.detrend(matrix, axis=1, type="constant")
    modal_coordinate = weights @ matrix

    half_width = max(
        float(rdt["band_half_width_min_hz"]),
        float(rdt["band_half_width_fraction"]) * f0,
    )
    low = max(1e-6, f0 - half_width)
    high = min(0.999 * fs / 2.0, f0 + half_width)
    if not 0 < low < high < fs / 2.0:
        return _rdt_failure("RDT带通频带无效", possible_crosstalk)

    sos = signal.butter(int(rdt["filter_order"]), [low, high], btype="bandpass", fs=fs, output="sos")
    try:
        filtered = signal.sosfiltfilt(sos, modal_coordinate)
    except Exception as exc:
        return _rdt_failure(f"零相位滤波失败:{exc}", possible_crosstalk)

    trim = int(math.ceil(float(rdt["edge_trim_periods"]) * fs / f0))
    segment_length = int(round(float(rdt["segment_periods"]) * fs / f0)) + 1
    if 2 * trim + segment_length >= filtered.size:
        return _rdt_failure("去边界后不足一个RDT片段", possible_crosstalk)

    valid = filtered[trim : filtered.size - trim]
    sigma = float(np.std(valid, ddof=1))
    threshold = float(rdt["trigger_sigma"]) * sigma
    if not np.isfinite(sigma) or sigma <= np.finfo(float).eps:
        return _rdt_failure("模态坐标标准差无效", possible_crosstalk)

    crossings = np.flatnonzero((valid[:-1] < threshold) & (valid[1:] >= threshold)) + 1
    crossings = crossings[crossings + segment_length <= valid.size]
    trigger_count = int(crossings.size)
    if trigger_count < int(rdt["minimum_segments"]):
        return _rdt_failure("RDT片段数不足", possible_crosstalk, trigger_count=trigger_count)

    decay = np.mean(np.stack([valid[i : i + segment_length] for i in crossings]), axis=0)
    decay = decay / max(abs(float(decay[0])), np.finfo(float).eps)
    times = np.arange(decay.size, dtype=float) / fs

    min_distance = max(1, int(round(0.65 * fs / f0)))
    peaks, _ = signal.find_peaks(decay, distance=min_distance)
    peaks = peaks[peaks > 0]
    max_peaks = int(rdt["fit_peak_count_max"])
    peaks = peaks[:max_peaks]
    fit_peak_count = int(peaks.size)
    if fit_peak_count < int(rdt["fit_peak_count_min"]):
        return _rdt_failure(
            "自由衰减峰数不足",
            possible_crosstalk,
            trigger_count=trigger_count,
            fit_peak_count=fit_peak_count,
        )

    amplitudes = np.abs(decay[peaks])
    valid_amp = np.isfinite(amplitudes) & (amplitudes > np.finfo(float).eps)
    peaks = peaks[valid_amp]
    amplitudes = amplitudes[valid_amp]
    fit_peak_count = int(peaks.size)
    if fit_peak_count < int(rdt["fit_peak_count_min"]):
        return _rdt_failure(
            "有效自由衰减峰数不足",
            possible_crosstalk,
            trigger_count=trigger_count,
            fit_peak_count=fit_peak_count,
        )

    candidates = []
    for count in range(int(rdt["fit_peak_count_min"]), fit_peak_count + 1):
        candidate_peaks = peaks[:count]
        candidate_t = times[candidate_peaks]
        candidate_y = np.log(amplitudes[:count])
        candidate_slope, candidate_intercept = np.polyfit(candidate_t, candidate_y, 1)
        candidate_predicted = candidate_slope * candidate_t + candidate_intercept
        candidate_res = float(np.sum((candidate_y - candidate_predicted) ** 2))
        candidate_tot = float(np.sum((candidate_y - candidate_y.mean()) ** 2))
        candidate_r2 = 1.0 - candidate_res / candidate_tot if candidate_tot > np.finfo(float).eps else 0.0
        candidates.append((candidate_r2, count, candidate_slope, candidate_intercept))

    r2, fit_peak_count, slope, intercept = max(candidates, key=lambda item: (item[0], -item[1]))
    peaks = peaks[:fit_peak_count]
    amplitudes = amplitudes[:fit_peak_count]
    fit_t = times[peaks]

    zeros = np.flatnonzero(np.signbit(decay[:-1]) != np.signbit(decay[1:]))
    if zeros.size >= 3:
        zero_times = []
        for idx in zeros:
            y0, y1 = float(decay[idx]), float(decay[idx + 1])
            fraction = abs(y0) / max(abs(y0) + abs(y1), np.finfo(float).eps)
            zero_times.append((idx + fraction) / fs)
        half_period = float(np.median(np.diff(zero_times)))
        f_refined = 1.0 / (2.0 * half_period) if half_period > 0 else float("nan")
    else:
        f_refined = float("nan")

    alpha = -float(slope)
    zeta = (
        alpha
        / math.sqrt((2 * math.pi * f_refined) ** 2 + alpha**2)
        if alpha > 0 and np.isfinite(f_refined) and f_refined > 0
        else float("nan")
    )

    intervals = np.diff(fit_t)
    spacing_cv = (
        float(np.std(intervals, ddof=1) / np.mean(intervals))
        if intervals.size > 1 and np.mean(intervals) > 0
        else 0.0
    )
    increases = amplitudes[1:] / np.maximum(amplitudes[:-1], np.finfo(float).eps) - 1.0
    max_increase = float(np.max(increases)) if increases.size else 0.0
    beat = bool(
        spacing_cv > float(rdt.get("beat_peak_spacing_cv_max", 0.20))
        or max_increase > float(rdt.get("beat_amplitude_increase_max", 0.30))
    )

    reasons = []
    if not np.isfinite(zeta):
        reasons.append("自由衰减未得到正阻尼")
    if r2 < float(rdt["r2_min"]):
        reasons.append("拟合R²低于门限")
    if np.isfinite(zeta) and not (
        float(rdt["zeta_min"]) <= zeta <= float(rdt["zeta_max"])
    ):
        reasons.append("阻尼比超出合理区间")
    if beat:
        reasons.append("拍频畸变")

    result = {
        "f_refined_Hz": f_refined,
        "DampingRatio": zeta,
        "Fit_R2": r2,
        "TriggerCount": trigger_count,
        "FitPeakCount": fit_peak_count,
        "PassedQualityGate": not reasons,
        "FailureReason": ";".join(reasons),
        "PossibleCrosstalk": bool(possible_crosstalk),
        "BeatDistortion": beat,
        "PeakSpacingCV": spacing_cv,
        "MaxEnvelopeIncrease": max_increase,
    }
    if return_curve:
        result["_curve"] = decay
        result["_times"] = times
        result["_fit_peaks"] = peaks
        result["_fit_line"] = np.exp(intercept + slope * times)
    return result


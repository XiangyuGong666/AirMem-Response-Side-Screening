"""V19.0 reproducible response-side analysis. Identification never receives a target frequency.

Only monitoring records and the uncorrected FE model enter identification.
Published tables are verification references, never a computation fallback.
"""
from __future__ import annotations

import hashlib
import json
import platform
from contextlib import nullcontext
import zipfile
import io
from pathlib import Path
from datetime import datetime, timedelta, timezone

import h5py
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.signal import find_peaks, get_window

from . import __version__
from .damping import estimate_efdd_window, estimate_rdt_window

ASSETS = Path(__file__).parent / 'model_inputs'
CONFIG = Path(__file__).parent / 'config'
ACC = 'analyacc_22_101112_dropbaseline_downsample10hz_12vertical.mat'
WIND = 'analywindspeed_downsample_10hz.mat'
CALENDAR = 'paper1_calendar_2022Q4.csv'
SENSORS = np.array([1217,1218,2499,4527,4519,4526,6668,6741,6667,8610,9740,9741])-1


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''): h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                                   default=lambda x: x.item() if hasattr(x,'item') else str(x)), encoding='utf-8')


def load_cfg():
    return json.loads((CONFIG/'analysis.json').read_text(encoding='utf-8'))



def monitoring_window(acceleration, window_id, samples, *, end_aligned=False):
    """Return the same float64 demeaned window for identification, damping and examples."""
    end = acceleration.shape[1] if end_aligned else window_id * samples
    x = np.asarray(acceleration[:, end-samples:end], dtype=np.float64)
    if x.shape != (12, samples) or not np.isfinite(x).all():
        raise ValueError(f'Invalid acceleration samples in window {window_id}')
    x -= x.mean(axis=1, keepdims=True)
    return x


def peak_csd(spectrum, frequencies, peak_hz, halfwidth_hz):
    """Dominant eigenvector of the CSD band, up to its common PSD/line-count scale."""
    z = spectrum[:, abs(frequencies-peak_hz)<=halfwidth_hz]
    matrix = z.conj() @ z.T
    _, vectors = np.linalg.eigh(matrix)
    u = vectors[:, -1]
    u *= np.exp(-1j*np.angle(u[np.argmax(abs(u))]))
    return matrix,u


def group_frequencies(freq, tol):
    """Observed-frequency groups using the current occurrence-weighted median."""
    vals, counts = np.unique(freq, return_counts=True)
    batches=[]; current=[]; cc=[]; med=None
    for v,c in zip(vals,counts):
        if med is not None and abs(v-med)>tol+1e-10:
            batches.append(current); current=[]; cc=[]
        current.append(v); cc.append(c)
        cum=np.cumsum(cc); n=int(cum[-1]); vv=np.asarray(current)
        med=(vv[np.searchsorted(cum,(n-1)//2+1)]+vv[np.searchsorted(cum,n//2+1)])/2
    if current: batches.append(current)
    labels=np.empty(len(freq),int)
    for gid,vs in enumerate(batches,1): labels[np.isin(freq,vs)]=gid
    return labels


def group_metrics(events, vectors, window_count, labels, scratch):
    rows=[]; reps=[]
    for gid in np.unique(labels):
        ix=np.flatnonzero(labels==gid); s=events.iloc[ix]; u=vectors[ix]
        _,v=np.linalg.eigh(u.T@u.conj()/len(ix)); rep=v[:,-1]
        rep*=np.exp(-1j*np.angle(rep[np.argmax(abs(rep))])); reps.append(rep)
        count=len(ix)*(len(ix)-1)//2
        # Exact population median; disk-backed storage bounds RAM for large groups.
        if count:
            tmp=Path(scratch)/'pairwise_mac.tmp'
            pairs=np.memmap(tmp,mode='w+',dtype='float64',shape=(count,)); pos=0
            for start in range(0,len(ix),128):
                block=abs(u[start:start+128].conj()@u.T)**2
                for j,line in enumerate(block):
                    tail=line[start+j+1:]; pairs[pos:pos+len(tail)]=tail; pos+=len(tail)
            shape=float(np.median(pairs,overwrite_input=True)); del pairs; tmp.unlink()
        else: shape=0.
        valid=s[s.FE_Mode>0]; votes=valid.FE_Mode.value_counts()
        mode=int(votes.index[0]) if len(votes) else 0
        consistency=float(votes.iloc[0]/len(valid)) if len(valid) else 0.
        nwin=int(s.WindowID.nunique()); occ=nwin/window_count
        rows.append(dict(PeakGroupID=int(gid),MedianPeakFreq_Hz=float(s.MeasFreq_Hz.median()),
            StdPeakFreq_Hz=float(s.MeasFreq_Hz.std()) if len(s)>1 else 0.,
            MinPeakFreq_Hz=float(s.MeasFreq_Hz.min()),MaxPeakFreq_Hz=float(s.MeasFreq_Hz.max()),
            OccurrenceCount=nwin,PeakEventCount=len(s),OccurrenceRate=occ,DominantFE_Mode=mode,
            ValidMatchedEvents=len(valid),FE_ModeConsistency=consistency,
            MedianMAC=float(valid.MAC.median()),MedianScore=float(valid.Score.median()),
            ShapePairCount=count,ShapeConsistencyMedian=shape,ThreeMetricScore=float(np.cbrt(occ*consistency*shape)),
            Votes=json.dumps({str(k):int(v) for k,v in votes.items()})))
    return pd.DataFrame(rows),np.array(reps)


def select_fundamental(groups,cfg,min_mode):
    ok=(groups.OccurrenceCount>=cfg['repetition_min_windows']) & groups.MedianPeakFreq_Hz.between(
        cfg['analysis_min_frequency_hz'],cfg['analysis_max_frequency_hz']) & (groups.ValidMatchedEvents>0)
    ranked=groups[ok].sort_values(['ThreeMetricScore','OccurrenceCount','MedianPeakFreq_Hz'],
        ascending=[False,False,True],kind='stable').copy()
    ranked['OverallRank']=np.arange(1,len(ranked)+1)
    eligible=ranked[ranked.DominantFE_Mode==min_mode]
    if eligible.empty: raise ValueError('No eligible group dominated by the lowest FE mode. No forced frequency fallback.')
    return ranked,eligible.iloc[0]


def load_model():
    nodes={}; mem=[]; cables=[]; typ=0
    for line in (ASSETS/'full_model_pure.mac').read_text(encoding='utf-8').splitlines():
        a=line.split(',')
        if a[0]=='n': nodes[int(a[1])]=list(map(float,a[2:5]))
        elif a[0]=='type': typ=int(a[1])
        elif a[0]=='e': (mem if typ==1 else cables).append([int(x) for x in a[1:] if x.strip()])
    xyz=np.array([nodes[i] for i in range(1,len(nodes)+1)])
    tri=np.array(mem)-1; seg=np.array(cables)-1
    cross=np.cross(xyz[tri[:,1]]-xyz[tri[:,0]],xyz[tri[:,2]]-xyz[tri[:,0]])
    area=np.linalg.norm(cross,axis=1)/2; normal=cross/(2*area[:,None])
    normal*=np.where(normal[:,2]<0,-1,1)[:,None]
    return xyz,tri,seg,area,normal


def model_mass(out, f0):
    xyz,tri,seg,area,normal=load_model()
    with np.load(ASSETS/'modes600.npz',allow_pickle=False) as z:
        U=z['U'][:585]
    bare=np.load(ASSETS/'my_freqs_585z.npy',allow_pickle=False)
    mu=np.empty(len(bare)); beta=np.empty(len(bare))
    for i,u in enumerate(U):
        em=abs(u[tri]).mean(1); mu[i]=np.sum(np.sum(em**2,axis=1)*area)
        signed=np.sum(u[tri].mean(1)*normal,axis=1)
        beta[i]=abs(np.sum(area*signed))/np.sum(area*abs(signed))
    ma=(bare[0]**2/f0**2-1)/mu[0]; pred=bare/np.sqrt(1+ma*mu)
    pd.DataFrame(dict(Mode=np.arange(1,len(bare)+1),Bare_Hz=bare,Mu=mu,Beta=beta,
                      Corrected_Hz=pred)).to_csv(out/'model_mass.csv',index=False)
    return dict(ma=float(ma),mu1=float(mu[0]),beta1=float(beta[0]),
                classic_ma=.68*1.225*115,beta_ma=float(.68*1.225*115*beta[0]))


def resolve_monitoring(path):
    root = Path(path).expanduser().resolve()
    for candidate in (root, root / '2022Q4'):
        if (candidate / ACC).is_file() and (candidate / WIND).is_file():
            return candidate
    raise FileNotFoundError(f'Select 2022Q4 or its parent. Missing {ACC} or {WIND} in {root}.')


def write_daily_coverage(intervals, out):
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            raise ValueError('Calendar coverage contains a nonpositive interval')
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    pd.DataFrame(merged, columns=['start', 'end']).to_csv(out/'coverage_intervals.csv', index=False)
    rows = []
    for start in pd.date_range('2022-10-01', '2022-12-31'):
        end = start + timedelta(days=1)
        seconds = sum((min(b, end) - max(a, start)).total_seconds()
                      for a, b in merged if b > start and a < end)
        rows.append(dict(date=start.date().isoformat(), coverage_h=seconds / 3600))
    pd.DataFrame(rows).to_csv(out/'coverage_daily.csv', index=False)


def calendar_coverage(path, out):
    calendar = pd.read_csv(path)
    starts = pd.to_datetime(calendar['calendar_start'], format='mixed')
    ends = pd.to_datetime(calendar['calendar_end'], format='mixed')
    if starts.empty or starts.isna().any() or ends.isna().any():
        raise ValueError('The public calendar index has missing timestamps')
    write_daily_coverage(zip(starts, ends), out)


def coverage(raw_root, out):
    """Rebuild calendar coverage from each original file's recorded time vector."""
    raw_root = Path(raw_root)
    archive_path = raw_root.is_file() and raw_root.suffix.lower() == '.zip'
    with zipfile.ZipFile(raw_root) if archive_path else nullcontext(None) as archive:
        files = sorted(f for f in archive.namelist() if f.lower().endswith('.npz')) if archive else sorted(raw_root.glob('*.mat'))
        if not files:
            raise FileNotFoundError('Expected timestamped NPZ time vectors or original MAT files.')
        rows=[]; intervals=[]
        for f in files:
            name=Path(f).name
            start=datetime.strptime(Path(name).stem,'%Y-%m-%d %H-%M-%S')
            if archive:
                with np.load(io.BytesIO(archive.read(f)), allow_pickle=False) as z:
                    t=np.asarray(z['time'],float).ravel()
            else:
                t=np.asarray(loadmat(f,variable_names=['time'])['time'],float).ravel()
            dd=np.diff(t); positive=dd[dd>0]
            if not np.isfinite(t).all() or not len(positive):
                raise ValueError(f'Invalid original time vector: {name}')
            dt=float(np.median(positive)); bad=np.flatnonzero((dd<=0)|(dd>1.5*dt))
            cuts=np.r_[0,bad+1,len(t)]
            for a,z in zip(cuts[:-1],cuts[1:]):
                intervals.append((start+timedelta(seconds=float(t[a])),start+timedelta(seconds=float(t[z-1]+dt))))
            rows.append(dict(file=name,samples=len(t),sample_dt_s=dt,time_first=float(t[0]),time_last=float(t[-1]),
                             time_sha256=hashlib.sha256(t.tobytes()).hexdigest()))
    pd.DataFrame(rows).to_csv(out/'coverage_files.csv',index=False)
    write_daily_coverage(intervals, out)


def run(data_root, out, *, original_root=None, progress=None):
    """Recompute the manuscript analysis from the two 2022Q4 monitoring MAT inputs."""
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=1, user_api='blas'):
        return _run(data_root, out, original_root, progress or _print_progress)


def _print_progress(percent, message):
    print(f'{percent:3d}% {message}', flush=True)


def _run(data_root, out, original_root, progress):
    out=Path(out).resolve(); out.mkdir(parents=True,exist_ok=True)
    # Never mix an incomplete run with previous figures or tables.
    if any(out.iterdir()): raise FileExistsError(f'Use an empty output directory: {out}')
    data_root = resolve_monitoring(data_root)
    cfg=load_cfg(); N=cfg['window_samples']; fs=cfg['sampling_hz']
    dc=json.loads((CONFIG/'damping.json').read_text(encoding='utf-8'))
    provenance = dict(manuscript='V19.0', period='2022Q4', software_version=__version__,
                      run_schema=2,
                      started_utc=datetime.now(timezone.utc).isoformat(),
                      status='running', config=cfg, damping_config=dc,
                      python=platform.python_version(), numpy=np.__version__, inputs={},
                      analysis_variant='per_window_demeaned')
    write_json(out/'run.json',provenance)
    try:
        for p in [data_root/ACC, data_root/WIND, data_root/CALENDAR,
                  *(CONFIG / name for name in ('analysis.json','damping.json','energy.json')),
                  *(ASSETS / name for name in ('FE_modal_sensor_z585.csv','full_model_pure.mac','modes600.npz','my_freqs_585z.npy'))]:
            if p.is_file():
                progress(2, f'Checking input {p.name}')
                provenance['inputs'][p.name]=dict(path=str(p.resolve()),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=digest(p))
        if original_root and Path(original_root).is_file():
            p=Path(original_root)
            provenance['inputs']['original_time_vectors']=dict(path=str(p.resolve()),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=digest(p))
        provenance['code_sha256'] = {p.name: digest(p) for p in Path(__file__).parent.glob('*.py')}
        write_json(out/'run.json',provenance)
        fe=pd.read_csv(ASSETS/'FE_modal_sensor_z585.csv'); ff=fe.freq_hz.to_numpy().reshape(-1,12)[:,0]
        phi=fe.uz.to_numpy().reshape(-1,12); phi/=np.linalg.norm(phi,axis=1)[:,None]
        wind=loadmat(data_root/WIND,variable_names=['analywindspeed_downsample'])['analywindspeed_downsample'].ravel()
        with h5py.File(data_root/ACC) as h:
            acc=h['analyacc_downsample']
            if acc.ndim != 2 or acc.shape[0] != 12 or acc.shape[1] != len(wind):
                raise ValueError('Expected 12 acceleration channels aligned with the wind record')
            sample_count=acc.shape[1]; total,remainder=divmod(sample_count,N)
            if not total or not np.isfinite(wind).all():
                raise ValueError('Wind data must be finite and contain at least one complete window')
            mean_wind=float(np.mean(wind,dtype=np.float64))
            wmean=wind[:total*N].reshape(total,N).mean(1)
            starts=np.arange(total,dtype=np.int64)*N
            if remainder:
                wmean=np.r_[wmean,wind[-N:].mean()]
                starts=np.r_[starts,sample_count-N]
            del wind
            selected_windows=wmean>=cfg['wind_threshold_mps']; ids=np.flatnonzero(selected_windows)+1
            pd.DataFrame(dict(WindowID=np.arange(1,len(wmean)+1),StartSample=starts,
                EndAligned=np.arange(len(wmean))>=total,WindMean_mps=wmean,
                Selected=selected_windows)).to_csv(out/'window_screening.csv',index=False)
            if not len(ids): raise ValueError('No high-wind windows')
            wind=wmean[ids-1]
            freq=np.fft.rfftfreq(N,1/fs); hann=get_window('hann',N)
            fac=np.full(len(freq),2/(fs*np.sum(hann**2))); fac[[0,-1]]/=2
            rows=[]; vectors=[]; psds=[]
            for k,(wid,ws) in enumerate(zip(ids,wind)):
                x=monitoring_window(acc,wid,N,end_aligned=wid>total)
                z=np.fft.rfft(x*hann,axis=1); p=(abs(z)**2*fac).mean(0); psds.append(p)
                pk,_=find_peaks(p,prominence=0); pk=pk[np.argsort(-p[pk],kind='stable')[:cfg['peaks_per_window']]]
                for j,a in enumerate(pk,1):
                    _,u=peak_csd(z,freq,freq[a],cfg['csd_halfwidth_hz'])
                    mac=abs(phi@u)**2; delta=abs(ff-freq[a])/freq[a]
                    score=np.where(delta<=cfg['max_relative_fe_error'],mac*np.exp(-delta**2/(2*cfg['frequency_weight_sigma']**2)),-np.inf)
                    best=int(score.argmax()); valid=np.isfinite(score[best])
                    rows.append(dict(WindowID=int(wid),PeakID=j,MeasFreq_Hz=float(freq[a]),WindMean_mps=float(ws),
                        PeakAmplitude=float(p[a]),FE_Mode=best+1 if valid else 0,MAC=float(mac[best]) if valid else np.nan,
                        Score=float(score[best]) if valid else np.nan)); vectors.append(u)
                if k % 1000 == 0:
                    progress(5 + 42 * k // len(ids), f'Scanning monitoring windows {k:,}/{len(ids):,}')
        events=pd.DataFrame(rows); U=np.array(vectors); PSD=np.array(psds)
        labels=group_frequencies(events.MeasFreq_Hz.to_numpy(),cfg['group_tolerance_hz']); events['PeakGroupID']=labels
        progress(49, 'Grouping peaks and calculating exact pairwise MAC medians')
        groups,reps=group_metrics(events,U,len(ids),labels,out)
        ranked,choice=select_fundamental(groups,cfg,1); gid=int(choice.PeakGroupID); base=events[events.PeakGroupID==gid]
        groups['OverallRank']=groups.PeakGroupID.map(ranked.set_index('PeakGroupID').OverallRank)
        groups['SelectedFundamental']=groups.PeakGroupID==gid
        events.to_csv(out/'peak_events.csv',index=False); groups.to_csv(out/'groups.csv',index=False)
        ranked.to_csv(out/'ranking.csv',index=False); base.to_csv(out/'selected_members.csv',index=False)
        np.savez_compressed(out/'spectra_vectors.npz',U=U,PSD=PSD,WindowID=ids,WindMean_mps=wind,freq=freq,representatives=reps)
        selected=choice.to_dict()
        selected['RepresentativeMode1MAC']=float(abs(np.vdot(phi[0],reps[gid-1]))**2)
        write_json(out/'identification.json',selected)
        progress(62, f"Selected {selected['MedianPeakFreq_Hz']:.3f} Hz after ranking")
        # Identification artifacts above are committed before model correction or damping.
        identification_hashes={p.name:digest(p) for p in [out/'peak_events.csv',out/'groups.csv',out/'ranking.csv',out/'identification.json']}
        f0=float(choice.MedianPeakFreq_Hz); mass=model_mass(out,f0)
        mass['std_relative']=2*(1+mass['ma']*mass['mu1'])/(mass['ma']*mass['mu1'])*float(choice.StdPeakFreq_Hz)/f0
        mass['sigma_interval']=[mass['ma']*(1-mass['std_relative']),mass['ma']*(1+mass['std_relative'])]
        model=pd.read_csv(out/'model_mass.csv'); comp=[]
        for _,q in groups[groups.MedianPeakFreq_Hz.between(f0,1.3)].iterrows():
            f=float(q.MedianPeakFreq_Hz); j=int(np.argmin(abs(model.Corrected_Hz-f))); b=int(np.argmin(abs(model.Bare_Hz-f)))
            comp.append(dict(PeakGroupID=int(q.PeakGroupID),Measured=f,Windows=int(q.OccurrenceCount),
                BareMode=b+1,Mode=j+1,Bare=float(model.Bare_Hz.iloc[b]),Corrected=float(model.Corrected_Hz.iloc[j]),
                BeforeError=100*(model.Bare_Hz.iloc[b]/f-1),AfterError=100*(model.Corrected_Hz.iloc[j]/f-1),Calibration=bool(q.SelectedFundamental)))
        pd.DataFrame(comp).to_csv(out/'frequency_comparison.csv',index=False)
        members = base.assign(Distance=abs(base.MeasFreq_Hz - f0))
        members = members.sort_values(['WindowID', 'Distance', 'PeakAmplitude'],
                                      ascending=[True, True, False], kind='stable').drop_duplicates('WindowID')
        dr=[]; dc['window_sec']=N/fs; dc['nperseg']=N
        with h5py.File(data_root/ACC) as h:
            for k, (ix, member) in enumerate(members.iterrows()):
                wid = int(member.WindowID)
                x=monitoring_window(h['analyacc_downsample'],wid,N,end_aligned=wid>total)
                ef=estimate_efdd_window(x,f0,float(choice.MinPeakFreq_Hz),float(choice.MaxPeakFreq_Hz),fs,dc)
                rd=estimate_rdt_window(x,f0,U[ix],fs,dc,True)
                row=dict(WindowID=wid, WindMean_mps=float(member.WindMean_mps),
                         PeakID=int(member.PeakID), PeakFreq_Hz=float(member.MeasFreq_Hz))
                for name,result in [('EFDD',ef),('RDT',rd)]:
                    row[name]=result['DampingRatio']; row[name+'_gate']=bool(result['PassedQualityGate'])
                    row[name+'_R2']=result.get('Fit_R2',np.nan); row[name+'_reason']=result.get('FailureReason','')
                dr.append(row)
                if k % 50 == 0:
                    progress(67 + 24 * k // len(members), f'EFDD and RDT {k:,}/{len(members):,}')
        dd=pd.DataFrame(dr); dd.to_csv(out/'damping_windows.csv',index=False); paired=dd[dd.EFDD_gate&dd.RDT_gate]
        rng=np.random.default_rng(0); stats={}
        for name in ['EFDD','RDT']:
            v=100*dd.loc[dd[name+'_gate'],name].to_numpy(); pv=100*paired[name].to_numpy()
            ci=np.percentile(np.median(rng.choice(pv,(10000,len(pv)),replace=True),axis=1),[5,95]).tolist() if len(pv) else [None,None]
            stats[name]=dict(passing=len(v),median=float(np.median(v)) if len(v) else None,
                std=float(np.std(v,ddof=1)) if len(v)>1 else None,
                paired_at_most_one_percent=int((pv<=1).sum()),
                paired_median=float(np.median(pv)) if len(pv) else None,CI90=ci)
        bins=[]
        for lo,hi in [(3,4),(4,5),(5,6),(6,7),(7,8),(8,np.inf)]:
            mask = (dd.WindMean_mps >= lo) & (dd.WindMean_mps < hi)
            row=dict(low=lo,high=hi,all_windows=int(((wind>=lo)&(wind<hi)).sum()),
                     fundamental_windows=int(mask.sum()))
            for name in ['EFDD','RDT']:
                vals=dd.loc[mask & dd[name+'_gate'],name]*100
                row[name+'_n']=len(vals)
                for key, quantile in [('q25', .25), ('median', .5), ('q75', .75)]:
                    row[name+'_'+key]=float(vals.quantile(quantile))
            bins.append(row)
        pd.DataFrame(bins).to_csv(out/'wind_bins.csv',index=False)
        # Prespecified descriptive neighborhoods, never used for selection.
        energy_cfg=json.loads((CONFIG/'energy.json').read_text())
        centers=np.array(energy_cfg['centers_hz']); band_low,band_high=energy_cfg['analysis_band_hz']
        halfwidth=energy_cfg['halfwidth_hz']; totalmask=(freq>=band_low-1e-10)&(freq<=band_high+1e-10)
        integrate=np.trapezoid if hasattr(np,'trapezoid') else np.trapz
        den=integrate(PSD[:,totalmask],freq[totalmask],axis=1)
        shares=[]
        for center in centers:
            mask=(np.arange(len(freq))>=int(round((center-halfwidth)*N/fs)))&(np.arange(len(freq))<=int(round((center+halfwidth)*N/fs)))
            shares.append(100*integrate(PSD[:,mask],freq[mask],axis=1)/den)
        np.savez_compressed(out/'energy.npz',centers=centers,shares=np.array(shares).T)
        progress(92, 'Calculating daily coverage')
        calendar_path = data_root / CALENDAR
        coverage_source = None
        if original_root is not None:
            coverage(original_root, out)
            coverage_source = 'original time vectors'
        elif calendar_path.is_file():
            calendar_coverage(calendar_path, out)
            coverage_source = 'public Zenodo calendar index'
        summary=dict(window_count=len(ids),events=len(events),matched=int((events.FE_Mode>0).sum()),
            groups=len(groups),eligible=len(ranked),selected=selected,mass=mass,damping=stats,paired=len(paired),
            damping_windows=len(dd), coverage_source=coverage_source,
            identification_hashes=identification_hashes,full_scan_windows=total,
            scanned_windows_including_partial=len(wmean),sample_count=sample_count,
            monitoring_hours=sample_count/fs/3600,mean_monitoring_wind_mps=mean_wind,
            end_aligned_window_mean_mps=float(wmean[-1]) if remainder else None,
            partial_window_policy='When samples remain, screen one final full-length window aligned to the record end using the same wind threshold')
        write_json(out/'summary.json',summary)
        table1=[]
        for label,lo,hi in [('<2 Hz',0,2),('2–3 Hz',2,3),('≥3 Hz',3,np.inf)]:
            ss=groups[(groups.MedianPeakFreq_Hz>=lo)&(groups.MedianPeakFreq_Hz<hi)]
            table1.append(dict(Band=label,Groups=len(ss),**{c:float(ss[c].median()) for c in ['OccurrenceRate','FE_ModeConsistency','ShapeConsistencyMedian','MedianScore']}))
        pd.DataFrame(table1).to_csv(out/'table1_band_statistics.csv',index=False)
        pd.DataFrame([selected]).to_csv(out/'table2_selected_statistics.csv',index=False)
        for name in (ACC,WIND,CALENDAR):
            if name not in provenance['inputs']: continue
            record=provenance['inputs'][name]; stat=(data_root/name).stat()
            if (stat.st_size,stat.st_mtime_ns)!=(record['bytes'],record['mtime_ns']):
                raise ValueError(f'Monitoring input changed during calculation: {name}')
        provenance['status']='complete'
        provenance['finished_utc']=datetime.now(timezone.utc).isoformat()
        provenance['outputs']={str(p.relative_to(out)):digest(p) for p in out.iterdir() if p.is_file() and p.name!='run.json'}
        write_json(out/'run.json',provenance)
        progress(100, 'Monitoring calculation complete')
        return out
    except Exception as exc:
        provenance['status']='failed'; provenance['error']=str(exc); write_json(out/'run.json',provenance); raise

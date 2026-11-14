"""Compare a completed public-data run with values printed in V19.0.

This independent checker is never imported by the analysis or figures.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def validate(root, require_figures=False):
    run=json.loads((root/'run.json').read_text(encoding='utf-8'))
    s=json.loads((root/'summary.json').read_text(encoding='utf-8'))
    checks=[]
    def check(label, actual, expected, atol=0):
        if atol:
            np.testing.assert_allclose(actual,expected,rtol=0,atol=atol,err_msg=label)
        else:
            np.testing.assert_equal(actual,expected,err_msg=label)
        checks.append(label)
    for name,expected_hash in run['outputs'].items():
        check(f'Analysis output digest: {name}',sha256(root/name),expected_hash)
    check('Samples, total scan, selected windows, peak events, matched events, groups, ranked groups',
          [s[k] for k in ['sample_count','scanned_windows_including_partial','window_count','events','matched','groups','eligible']],
          [67285219,112143,26399,131995,117006,118,106])
    sc=pd.read_csv(root/'window_screening.csv')
    check('Every scanned window appears in the screening ledger',len(sc),112143)
    check('Last end-aligned window: start, flag and high-wind decision',
          [int(sc.StartSample.iloc[-1]),bool(sc.EndAligned.iloc[-1]),bool(sc.Selected.iloc[-1])],
          [67284619,True,False])
    check('End-aligned window mean wind',sc.WindMean_mps.iloc[-1],1.1395,1.e-6)
    check('Full-record wind-speed mean and duration',[s['mean_monitoring_wind_mps'],s['monitoring_hours']],
          [2.2094401553520284,1869.0338611111111],5.e-7)
    e=pd.read_csv(root/'peak_events.csv')
    check('All screened windows contribute exactly five local peaks',e.groupby('WindowID').size().to_numpy(),np.full(26399,5))
    check('Peak ledger window IDs are those screened from wind',e.WindowID.unique(),sc.loc[sc.Selected,'WindowID'].to_numpy())
    g=s['selected']
    check('Group rank, occurrences, events, matched-mode vote counts and pair count',
          [g['OverallRank'],g['OccurrenceCount'],g['PeakEventCount'],*json.loads(g['Votes']).values(),g['ShapePairCount']],
          [1,1245,1261,1241,20,794430])
    check('Fundamental frequency and frequency spread', [g['MedianPeakFreq_Hz'],g['StdPeakFreq_Hz']],[.550,.015331],.0000005)
    check('Equal-peak representative vs FE mode 1 MAC',g['RepresentativeMode1MAC'],.25684,.000005)
    check('Equivalent areal added mass',s['mass']['ma'],1.33334,.000005)
    f=pd.read_csv(root/'frequency_comparison.csv');f=f[~f.Calibration]
    check('17 comparison groups outside calibration',len(f),17)
    check('Median nearest-frequency deviations before and after correction',
          f[['BeforeError','AfterError']].abs().median().tolist(),[3.99,2.01],.005)
    ef,rd=[s['damping'][k] for k in ('EFDD','RDT')]
    check('Damping windows, passes and paired fraction <=1%',
          [s['damping_windows'],ef['passing'],rd['passing'],s['paired'],ef['paired_at_most_one_percent'],rd['paired_at_most_one_percent']],
          [1245,1233,653,648,578,455])
    check('All-pass and paired medians (%)',
          [ef['median'],rd['median'],ef['paired_median'],rd['paired_median']],[.550,.752,.545,.749],.0005)
    w=pd.read_csv(root/'wind_bins.csv')
    check('Wind-bin counts',w[['all_windows','fundamental_windows','EFDD_n','RDT_n']].to_numpy(),
          [[17620,582,574,321],[6449,391,388,195],[1786,181,180,90],[442,73,73,41],[84,16,16,5],[18,2,2,1]])
    check('Wind-bin medians (%)',w[['EFDD_median','RDT_median']].to_numpy(),
          [[.568,.810],[.553,.695],[.490,.687],[.474,.623],[.550,.983],[.288,.592]],.0005)
    d=pd.read_csv(root/'damping_windows.csv').set_index('WindowID')
    base=pd.read_csv(root/'selected_members.csv')
    check('Damping covers every distinct selected-group window exactly once',d.index.to_numpy(),np.sort(base.WindowID.unique()))
    check('Example 30087 RDT (%)',100*d.loc[30087,'RDT'],.84,.005)
    if s['coverage_source']=='original time vectors':
        check('Original time-vector file count',len(pd.read_csv(root/'coverage_files.csv')),1875)
    manifest=root/'figure_manifest.json'
    if require_figures and not manifest.is_file():
        raise FileNotFoundError('No figure manifest; run plot all first')
    if manifest.exists():
        m=json.loads(manifest.read_text(encoding='utf-8'))
        expected=set(range(3,24))
        if not s['coverage_source']:expected.remove(3)
        check('Flowchart and site photograph excluded',sorted(m['excluded']),['1','2'])
        check('Complete English figure set (excluding 1 and 2)', sorted(map(int,m['figures'])),sorted(expected))
        run_hash=sha256(root/'run.json')
        for n,item in m['figures'].items():
            assert item['run_sha256']==run_hash,n
            for name,source in item['sources'].items():
                expected_hash=(run['outputs'][name] if source['origin']=='current calculation' else run['inputs'][name]['sha256'])
                assert source['sha256']==expected_hash,(n,name)
            assert item['sources'], n
            for suffix in ['pdf','png']:
                p=root/item[suffix]
                assert p.is_file() and p.stat().st_size>0 and sha256(p)==item['output_sha256'][suffix],p
        checks.append('Every figure and source fingerprint agrees with its completed run')
        if s['coverage_source']:
            assert m['figures']['3']['coverage_source']==s['coverage_source']
        counts=pd.read_csv(root/'figure09_density.csv').to_numpy().sum()
        check('Figure 9 histogram includes every recorded peak event',counts,len(e))
        if '14' in m['figures']:
            heat=pd.read_csv(root/'figure14_time_counts.csv',index_col=0)
            counts=e[e.PeakGroupID.isin(heat.index)].groupby('PeakGroupID').WindowID.nunique()
            check('Figure 14 time counts sum to distinct occurrences per group',heat.sum(axis=1).to_numpy(),counts.loc[heat.index].to_numpy())
    return checks


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    parser.add_argument('--require-figures',action='store_true')
    args=parser.parse_args()
    for label in validate(args.run,args.require_figures):print(f'PASS: {label}')

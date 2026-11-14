"""All computational figures of V19.0, drawn from this run's computed data.

Figures 1 and 2 are excluded. No manuscript images are read.
"""
from pathlib import Path
import json
import matplotlib
import h5py
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LogNorm
from scipy.signal import get_window

from .analysis import ASSETS, ACC, SENSORS, load_model, write_json, resolve_monitoring, digest, monitoring_window, peak_csd
from .damping import estimate_rdt_window, _leading_csd_curves

BLUE='#285f90'; ORANGE='#bc7132'; GREEN='#21855a'
TITLES={3:'Daily monitoring coverage',
 4:'FE geometry and accelerometers',5:'FE frequencies and spacing',6:'First six signed normal mode shapes',
 7:'Single-window PSD and five strongest peaks',8:'CSD matrix and dominant response vector',9:'All peak-event time-frequency density',
 10:'Uniform grouping and repeatability ranking',11:'FE matching consistency',12:'Selected-group frequency distribution',
 13:'Two single-window matching examples',14:'All low-frequency groups and occurrences',15:'Added mass and fundamental frequency',
 16:'Signed normal displacement of mode 1',17:'Frequency comparison with one added-mass parameter',
 18:'Median spectra by wind conditions',19:'Variance fractions in fixed neighborhoods',20:'Single-window damping estimation',
 21:'Passing damping distributions',22:'Damping estimates across wind bins',23:'Paired damping estimates and empirical CDF'}
STYLE = {'font.family': ['Times New Roman', 'DejaVu Serif'], 'font.size': 11,
         'pdf.fonttype': 42, 'mathtext.fontset': 'stix', 'axes.spines.top': False,
         'axes.spines.right': False}


EXCLUDED = {'1': 'Methods flowchart', '2': 'Site photograph'}

# Exact inputs read by each renderer. No other run outputs are touched.
PLOT_INPUTS = {
    3: ('summary.json','coverage_daily.csv'),
    4: ('full_model_pure.mac',),
    5: ('model_mass.csv',),
    6: ('full_model_pure.mac','modes600.npz','model_mass.csv'),
    7: ('spectra_vectors.npz','peak_events.csv'),
    8: (ACC,'peak_events.csv'),
    9: ('summary.json','peak_events.csv'),
    10: ('groups.csv','ranking.csv'),
    11: ('groups.csv',),
    12: ('groups.csv','peak_events.csv'),
    13: (ACC,'peak_events.csv','FE_modal_sensor_z585.csv'),
    14: ('groups.csv','peak_events.csv','frequency_comparison.csv'),
    15: ('summary.json','model_mass.csv','groups.csv'),
    16: ('full_model_pure.mac','modes600.npz','model_mass.csv'),
    17: ('frequency_comparison.csv','model_mass.csv'),
    18: ('spectra_vectors.npz',),
    19: ('energy.npz','spectra_vectors.npz','groups.csv'),
    20: (ACC,'groups.csv','peak_events.csv','damping_windows.csv'),
    21: ('damping_windows.csv',),
    22: ('wind_bins.csv',),
    23: ('damping_windows.csv',),
}
MODEL_FILES = {'full_model_pure.mac','modes600.npz','FE_modal_sensor_z585.csv'}


def read_run(out):
    run = json.loads((Path(out)/'run.json').read_text(encoding='utf-8'))
    if run.get('status') != 'complete' or run.get('manuscript') != 'V19.0' or run.get('run_schema') != 2:
        raise ValueError('Run a fresh V19.0 monitoring calculation with v0.1.2 or later.')
    return run


def subplots(*args, figsize=None, layout=None):
    fig = Figure(figsize=figsize, layout=layout)
    return fig, fig.subplots(*args)


def build(number, out, data_root=None):
    """Build a native Matplotlib figure from a completed calculation."""
    out = Path(out).resolve()
    if number not in TITLES:
        raise ValueError(f'Figure {number} is not part of this software.')
    run = read_run(out)
    if number==3 and 'coverage_daily.csv' not in run['outputs']:
        raise FileNotFoundError('Figure 3 requires the public calendar CSV or original time vectors.')
    sources={}; files={}
    for name in PLOT_INPUTS[number]:
        if name == ACC:
            record=run['inputs'][name]; expected=record['sha256']
            path=resolve_monitoring(data_root or Path(record['path']).parent)/ACC
            stat=path.stat()
            unchanged=(path==Path(record['path']) and stat.st_mtime_ns==record.get('mtime_ns') and stat.st_size==record['bytes'])
            origin='monitoring acceleration'
        elif name in MODEL_FILES:
            path=ASSETS/name; expected=run['inputs'][name]['sha256']; unchanged=False
            origin='bare FE input'
        else:
            path=out/name; expected=run['outputs'].get(name); unchanged=False
            origin='current calculation'
        if not path.is_file() or expected is None:
            raise ValueError(f'Figure {number}: required input {name} is missing from this calculation.')
        if not unchanged and digest(path)!=expected:
            raise ValueError(f'Figure {number}: {name} differs from the completed calculation. Recalculate from monitoring data.')
        files[name]=path; sources[name]={'origin':origin,'sha256':expected}
    with matplotlib.rc_context(STYLE):
        fig=_build(number,out,run,files)
    fig.airmem_sources=sources
    return fig


def _build(number, out, run, files):
    def csv(name):
        return pd.read_csv(files[name])
    def axes(n=2,size=(10,4)):
        fig,aa=subplots(1,n,figsize=size,layout='constrained'); return fig,np.atleast_1d(aa)
    def example(wid, target=None, peak_id=None):
        with h5py.File(files[ACC]) as h:
            x=monitoring_window(h['analyacc_downsample'],wid,window_samples)
        fft=np.fft.rfft(x*get_window('hann',window_samples),axis=1)
        fr=np.fft.rfftfreq(window_samples,1/fs)
        ss=e[e.WindowID==wid]
        if ss.empty: raise ValueError(f'Example window {wid} absent from this run')
        row=ss[ss.PeakID==peak_id].iloc[0] if peak_id is not None else ss.iloc[np.argmin(abs(ss.MeasFreq_Hz.to_numpy()-target))]
        G,u=peak_csd(fft,fr,row.MeasFreq_Hz,run['config']['csd_halfwidth_hz'])
        return x,row,u,G
    window_samples=run['config']['window_samples']; fs=run['config']['sampling_hz']
    windows_per_hour=3600*fs/window_samples
    if 'summary.json' in files:
        s=json.loads(files['summary.json'].read_text(encoding='utf-8'))
    if 'groups.csv' in files:
        g=csv('groups.csv'); selected=g[g.SelectedFundamental].iloc[0]; f0=float(selected.MedianPeakFreq_Hz)
    if 'peak_events.csv' in files:e=csv('peak_events.csv')
    if 'model_mass.csv' in files:model=csv('model_mass.csv')
    if 'frequency_comparison.csv' in files:comp=csv('frequency_comparison.csv')
    if 'damping_windows.csv' in files:d=csv('damping_windows.csv')
    if 'spectra_vectors.npz' in files:
        with np.load(files['spectra_vectors.npz'],allow_pickle=False) as z:
            if number in (7,18): freq=z['freq']; psd=z['PSD']
            if number==7: ids=z['WindowID']
            if number in (18,19): wind=z['WindMean_mps']
    if number==3:
        cov=csv('coverage_daily.csv'); fig,aa=axes(1,(10,3.5)); a=aa[0]
        a.bar(pd.to_datetime(cov.date),cov.coverage_h,color=BLUE,width=.85); a.axhline(24,ls='--',color='gray')
        a.set(xlabel='Date (2022)',ylabel='Coverage / h',ylim=(0,26),title="Coverage from the public calendar index" if s['coverage_source'] != 'original time vectors' else None)
    elif number in [4,6,16]:
        xyz,tri,seg,area,normal=load_model(); local=xyz-(xyz.min(0)+xyz.max(0))/2
        if number==4:
            fig,aa=subplots(2,1,figsize=(10,6),layout='constrained')
            for i,a in enumerate(aa):
                if i==0:
                    az=np.deg2rad(-65); el=np.deg2rad(25); right=np.array([-np.sin(az),np.cos(az),0]); up=np.array([-np.sin(el)*np.cos(az),-np.sin(el)*np.sin(az),np.cos(el)])
                    xy=np.c_[local@right,local@up]; order=np.argsort((local@np.cross(right,up))[tri].mean(1))
                else: xy=local[:,:2]; order=np.argsort(local[tri,2].mean(1))
                a.add_collection(PolyCollection(xy[tri][order],facecolors='#dbe7ee',edgecolors='#afc2cf',linewidths=.04))
                if i:a.add_collection(LineCollection(xy[seg],colors=BLUE,linewidths=.12))
                a.scatter(*xy[SENSORS].T,c='#f9d547',s=25,edgecolors='black',lw=.4,zorder=4)
                for k,p in enumerate(xy[SENSORS]):a.annotate(str(k+1),p,xytext=(0,7 if k%3 else -12),textcoords='offset points',ha='center',fontsize=9,zorder=5)
                a.autoscale(); a.margins(.06,.18); a.set_aspect('equal'); a.axis('off')
        else:
            with np.load(files['modes600.npz'],allow_pickle=False) as z: modes=z['U']
            fig,aa=subplots(2,3,figsize=(12,6),layout='constrained') if number==6 else subplots(1,1,figsize=(10,4),layout='constrained')
            for i,a in enumerate(np.atleast_1d(aa).flat):
                signed=np.sum(modes[i][tri].mean(1)*normal,axis=1); signed/=max(abs(signed))
                pc=PolyCollection(local[tri,:2],array=signed,cmap='RdBu_r',edgecolors='none',clim=(-1,1));a.add_collection(pc)
                a.autoscale();a.set_aspect('equal');a.axis('off');a.set_title(f'Mode {i+1}: {model.Bare_Hz.iloc[i]:.4f} Hz')
                pd.DataFrame(dict(Element=np.arange(1,len(tri)+1),SignedNormal=signed)).to_csv(out/f'figure{number:02d}_mode{i+1}_data.csv',index=False)
            fig.colorbar(pc,ax=list(np.atleast_1d(aa).flat),shrink=.8,label='Normalized normal displacement')
    elif number==5:
        fig,aa=axes(); f=model.loc[model.Bare_Hz<=5,'Bare_Hz'].to_numpy()
        aa[0].plot(np.arange(1,len(f)+1),f,'.',ms=2,color=BLUE); aa[0].set(xlabel='Mode order',ylabel='Frequency / Hz')
        gaps=np.diff(f); limit=.020
        display=gaps[gaps<=limit]
        aa[1].hist(display,bins=np.linspace(0,limit,41),color=BLUE,edgecolor='white',linewidth=.5)
        for stat,c,label in [(np.median(gaps),ORANGE,'Median'),(np.mean(gaps),BLUE,'Mean')]:
            aa[1].axvline(stat,c=c,ls='--',label=f'{label} {stat:.5f} Hz')
        aa[1].legend(frameon=False,fontsize=9); aa[1].set(xlabel='Adjacent spacing / Hz',ylabel='Count',title=f'(b) {len(gaps):,} gaps; {len(gaps)-len(display)} above plot range')
        aa[0].set_title(f'(a) {len(f):,} modes at or below 5 Hz')
    elif number in [7,8]:
        wid=107431
        if number==7:
            fig,aa=axes(1);a=aa[0];k=int(np.flatnonzero(ids==wid)[0]);a.plot(freq,psd[k],c=BLUE)
            ss=e[e.WindowID==wid];a.vlines(ss.MeasFreq_Hz,0,max(psd[k])*1.04,color=ORANGE,linestyles='dotted'); a.scatter(ss.MeasFreq_Hz,ss.PeakAmplitude,c=ORANGE,s=12,zorder=4)
            a.set(xlabel='Frequency / Hz',ylabel='Mean PSD (record-unit²/Hz)',title=f'Window {wid}')
            pd.DataFrame(dict(Frequency_Hz=freq,PSD=psd[k])).to_csv(out/'figure07_data.csv',index=False)
        else:
            x,row,u,G=example(wid,target=33/60)
            fig,aa=axes(); im=aa[0].imshow(abs(G)/abs(G).max(),origin='lower',cmap='viridis',extent=(.5,12.5,.5,12.5));fig.colorbar(im,ax=aa[0],label='Normalized CSD magnitude')
            aa[1].bar(np.arange(1,13),abs(u)/max(abs(u)),color=BLUE);aa[1].set(xlabel='Sensor',ylabel='Dominant-vector magnitude')
            pd.DataFrame(abs(G)/abs(G).max()).to_csv(out/'figure08_csd_data.csv',index=False)
    elif number==9:
        fig,aa=axes(1);a=aa[0];h,xe,ye=np.histogram2d((e.WindowID-1)/windows_per_hour,e.MeasFreq_Hz,bins=[180,150],range=[[0,s['full_scan_windows']/windows_per_hour],[0,5]])
        im=a.pcolormesh(xe,ye,np.ma.masked_equal(h.T,0),norm=LogNorm(vmin=1,vmax=h.max()),cmap='viridis');fig.colorbar(im,ax=a,label='Peak-event count')
        a.set(xlabel='Concatenated valid time / h',ylabel='Peak frequency / Hz',title=f'{len(e):,} peak events')
        pd.DataFrame(h).to_csv(out/'figure09_density.csv',index=False)
    elif number==10:
        fig,aa=axes();aa[0].scatter(g.MedianPeakFreq_Hz,g.ThreeMetricScore,s=18,color=BLUE);aa[0].scatter(f0,selected.ThreeMetricScore,marker='*',s=150,c=ORANGE,zorder=4)
        aa[0].set(xlabel='Group median frequency / Hz',ylabel='Three-repeatability score',xlim=(0,5))
        top=csv('ranking.csv').head(10).sort_values('ThreeMetricScore');aa[1].barh(np.arange(len(top)),top.ThreeMetricScore,color=[ORANGE if v==selected.PeakGroupID else BLUE for v in top.PeakGroupID]);aa[1].set_yticks(np.arange(len(top)),[f'{x:.3f}' for x in top.MedianPeakFreq_Hz]);aa[1].set(xlabel='Three-repeatability score',ylabel='Group median frequency / Hz')
    elif number==11:
        fig,aa=axes(1);a=aa[0];gg=g[g.ValidMatchedEvents>0];a.scatter(gg.MedianPeakFreq_Hz,gg.FE_ModeConsistency,c=BLUE,s=22);a.scatter(f0,selected.FE_ModeConsistency,marker='*',c=ORANGE,s=150,zorder=4)
        a.set(xlabel='Group median frequency / Hz',ylabel='FE matching consistency',ylim=(0,1.05),xlim=(0,5)); a.annotate(f'{f0:.3f} Hz\n{selected.FE_ModeConsistency:.2%}',(f0,selected.FE_ModeConsistency),xytext=(12,-16),textcoords='offset points')
    elif number==12:
        fig,aa=axes(1);a=aa[0];base=e[e.PeakGroupID==selected.PeakGroupID];counts=base.MeasFreq_Hz.value_counts().sort_index();a.bar(counts.index,counts.values,width=.011,color=BLUE);a.axvline(f0,c=ORANGE,ls='--',label=f'Median {f0:.3f} Hz');a.legend(frameon=False);a.set(xlabel='Peak frequency / Hz',ylabel='Peak-event count',title=f'{len(base)} events; {base.WindowID.nunique()} windows')
    elif number==13:
        fig,aa=subplots(2,2,figsize=(10,6),layout='constrained');fe=csv('FE_modal_sensor_z585.csv').uz.to_numpy().reshape(-1,12)[0].astype(complex)
        for col,(wid,target) in enumerate([(107431, 0.55), (22369, 32/60)]):
            x,row,u,G=example(wid, target);u*=np.exp(-1j*np.angle(np.vdot(fe,u)));u/=max(abs(u));v=fe/max(abs(fe))
            for a,p,q in [(aa[0,col],abs(u),abs(v)),(aa[1,col],u.real,v.real)]:
                a.plot(range(1,13),p,'o-',color=BLUE,label='Measured');a.plot(range(1,13),q,'s--',color=ORANGE,label='FE mode 1');a.set(xlabel='Sensor',ylabel='Normalized component');a.legend(frameon=False)
            aa[0,col].set_title(f'Window {wid}, {row.MeasFreq_Hz:.4f} Hz; MAC={row.MAC:.3f}')
    elif number==14:
        cg=g[g.PeakGroupID.isin(comp.PeakGroupID)].sort_values('MedianPeakFreq_Hz');edges=np.arange(0,1921,60)
        heat=np.array([np.histogram((np.sort(e.loc[e.PeakGroupID==r.PeakGroupID,'WindowID'].unique())-1)/windows_per_hour,edges)[0] for _,r in cg.iterrows()])
        pd.DataFrame(heat,index=cg.PeakGroupID).to_csv(out/'figure14_time_counts.csv')
        fig,aa=axes(size=(11,7));y=np.arange(len(cg));labels=[f'{x:.3f}' for x in cg.MedianPeakFreq_Hz]
        aa[0].barh(y,cg.OccurrenceCount,color=[ORANGE if v==selected.PeakGroupID else BLUE for v in cg.PeakGroupID]);aa[0].set_yticks(y,labels);aa[0].invert_yaxis();aa[0].set(xlabel='Distinct occurrence windows',ylabel='Group frequency / Hz')
        im=aa[1].pcolormesh(edges,np.arange(len(cg)+1)-.5,heat,cmap='Blues');aa[1].set_yticks(y,labels);aa[1].invert_yaxis();aa[1].set(xlabel='Concatenated valid time / h');fig.colorbar(im,ax=aa[1],label='Windows / 60-h block')
    elif number==15:
        fig,aa=axes(1,(9,6));a=aa[0];ma=s['mass']['ma'];mu=s['mass']['mu1'];fb=model.Bare_Hz.iloc[0];calc=lambda x:fb/np.sqrt(1+np.asarray(x)*mu);xx=np.unique(np.r_[np.linspace(0,3,250),np.linspace(3,100,750)])
        for ax in [a,a.inset_axes([.54,.46,.43,.46])]:
            ax.plot(xx,calc(xx),c=BLUE,label='Predicted mode-1 frequency')
            ax.axhline(f0,c=ORANGE,ls='--',label='Selected-group median frequency')
            ax.scatter([ma],[f0],c=ORANGE,zorder=4,label=f'Calibrated mass {ma:.3f} kg/m²')
            for key,c in [('beta_ma',GREEN),('classic_ma','#ac4e36')]:
                label='β-reduced scale reference' if key=='beta_ma' else 'Unreduced classical value'
                ax.scatter([s['mass'][key]],[calc(s['mass'][key])],c=c,marker='s',label=label)
            ax.set(xlabel='Added mass / kg m$^{-2}$',ylabel='Frequency / Hz')
        a.legend(loc='upper left',bbox_to_anchor=(.14,.97),fontsize=8.5,frameon=True,facecolor='white',edgecolor='none',framealpha=1)
        a.set(xlim=(0,100),ylim=(.1,.68));a.child_axes[0].set(xlim=(0,2.5),ylim=(.49,.65));a.child_axes[0].tick_params(labelsize=8)
        pd.DataFrame(dict(AddedMass=xx,Frequency_Hz=calc(xx))).to_csv(out/'figure15_curve.csv',index=False)
    elif number==17:
        fig,aa=axes(size=(11,5))
        for y,v,c in [(0,model.Bare_Hz,ORANGE),(1,model.Corrected_Hz,BLUE),(2,comp.Measured,GREEN)]:aa[0].vlines(v[(v>=.48)&(v<=1.38)],y-.18,y+.18,color=c)
        aa[0].set_title('(a) Frequency locations');aa[0].set_yticks([0,1,2],['Bare FE','Added-mass corrected','Observed group medians']);aa[0].set(xlim=(.48,1.38),xlabel='Frequency / Hz')
        a=aa[1];xx=np.arange(len(comp));a.axhspan(-3,3,color='#eee');a.bar(xx-.18,comp.BeforeError,.36,color=ORANGE,label='Before');a.bar(xx+.18,comp.AfterError,.36,color=BLUE,label='After');a.set_xticks(xx,[f'{x:.3f}' for x in comp.Measured],rotation=65,ha='right');a.set(ylabel='Nearest-frequency deviation / %',xlabel='Measured group frequency / Hz',title='(b) Nearest-frequency deviations');a.legend(frameon=False)
    elif number==18:
        fig,aa=axes(1,(9,6));a=aa[0];ins=a.inset_axes([.44,.65,.53,.3]);masks=[(wind>=3)&(wind<7),wind>=7,wind>=6];meds=[np.median(psd[m],axis=0) for m in masks]
        for ax in [a,ins]:
            for med,c,label,m in zip(meds,[BLUE,ORANGE,GREEN],['3–7 m/s','≥7 m/s','≥6 m/s'],masks):ax.plot(freq,med,c=c,label=f'{label} (n={sum(m)})')
            ax.set(xlabel='Frequency / Hz',ylabel='Median PSD');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        a.set(xlim=(.1,5));a.set_ylim(0,a.get_ylim()[1]*1.8);ins.set(xlim=(.1,1.6));ins.tick_params(labelsize=8);a.legend(loc='upper left',frameon=False,fontsize=9)
        pd.DataFrame(dict(Frequency_Hz=freq,PSD_3to7=meds[0],PSD_ge7=meds[1],PSD_ge6=meds[2])).to_csv(out/'figure18_median_psd.csv',index=False)
    elif number==19:
        with np.load(files['energy.npz'],allow_pickle=False) as en: centers=en['centers']; shares=en['shares']
        fig,aa=axes(size=(11,4.5))
        for a,mask,title in [(aa[0],np.ones(len(wind),bool),'All analysis windows'),(aa[1],wind>=7,'Wind ≥7 m/s')]:
            vals=np.median(shares[mask],axis=0);order=np.argsort(-vals)[:10];a.bar(range(10),vals[order],color=[ORANGE if abs(centers[j]-f0)<.001 else BLUE for j in order]);a.set_xticks(range(10),[f'{centers[j]:.3f}' for j in order],rotation=50,ha='right');a.set(xlabel='Fixed neighborhood center / Hz',ylabel='Median in-band variance / %',title=title)
    elif number==20:
        wid=30087
        member=d.set_index('WindowID').loc[wid]; x,row,u,G=example(wid,peak_id=int(member.PeakID))
        dc=run['damping_config']; fs=run['config']['sampling_hz']; nperseg=x.shape[1]
        rd=estimate_rdt_window(x,f0,u,fs,dc,True,True)
        if '_curve' not in rd:
            raise ValueError(f'Window {wid} did not produce an RDT curve: {rd["FailureReason"]}')
        fr,sv,_=_leading_csd_curves(x,fs,nperseg,nperseg*int(dc['zero_padding_factor']))
        fig,aa=axes();a=aa[0];t=rd['_times'];curve=rd['_curve'];fit=rd['_fit_line'];peaks=rd['_fit_peaks']
        a.plot(t,curve,c=BLUE,label='RD signature');a.plot(t,fit,'--',c=ORANGE,label='Fitted envelope')
        a.plot(t[peaks],curve[peaks],'o',c=ORANGE);a.legend(frameon=False,loc='lower right');a.set(xlabel='Time / s',ylabel='Normalized response',title=f'(a) RDT: {rd["DampingRatio"]*100:.2f}%')
        fraction=dc['local_band_fraction']; mask=(fr>=f0*(1-fraction))&(fr<=f0*(1+fraction))
        aa[1].plot(fr[mask],sv[mask],c=BLUE);aa[1].set(xlabel='Frequency / Hz',ylabel='First singular value',title=f'(b) EFDD: window {wid}')
        pd.DataFrame(dict(Time_s=t,RD=curve,EnvelopeFit=fit)).to_csv(out/'figure20_rd_data.csv',index=False)
        pd.DataFrame(dict(Frequency_Hz=fr[mask],FirstSingularValue=sv[mask])).to_csv(out/'figure20_efdd_spectrum.csv',index=False)
    elif number==21:
        fig,aa=axes()
        for i,(a,name,c) in enumerate(zip(aa,['EFDD','RDT'],[BLUE,ORANGE])):
            v=d.loc[d[name+'_gate'],name].to_numpy()*100
            all_pass=np.r_[d.loc[d.EFDD_gate,'EFDD'], d.loc[d.RDT_gate,'RDT']]*100
            upper=max(6.2,np.ceil(all_pass.max()*5)/5) if len(all_pass) else 6.2
            if len(v):
                q1,md,q3=np.percentile(v,[25,50,75]); a.hist(v,bins=np.arange(0,upper+.1,.2),weights=np.full(len(v),100/len(v)),color=c,edgecolor='white',linewidth=.55)
                a.axvspan(q1,q3,color=c,alpha=.16);a.axvline(md,ls='--',c=c,label=f'Median {md:.3f}%');a.legend(frameon=False)
            a.set(xlabel='Damping ratio / %',ylabel='Passing-window fraction / %',title=f'({chr(97+i)}) {name}, n = {len(v):,}',xlim=(0,upper),ylim=(0,max(40,a.get_ylim()[1])))
    elif number==22:
        bins=csv('wind_bins.csv');fig,(a,b)=subplots(2,1,figsize=(10,7),layout='constrained'); xx=np.arange(len(bins))
        bars=[('Fundamental-group windows','fundamental_windows','#acb4bd'),('EFDD passes','EFDD_n',BLUE),('RDT passes','RDT_n',ORANGE)]
        for sub,(label,col,c) in enumerate(bars):
            bb=b.bar(xx+(sub-1)*.23,bins[col],.23,color=c,label=label)
            b.bar_label(bb,labels=[f'{n:,}' if n else '' for n in bins[col]],padding=2,fontsize=8)
        for name,c,m,shift in [('EFDD',BLUE,'o',-.085),('RDT',ORANGE,'s',.085)]:
            mid=bins[name+'_median']; a.errorbar(xx+shift,mid,yerr=np.maximum(0,np.stack([mid-bins[name+'_q25'],bins[name+'_q75']-mid])),fmt=m,c=c,capsize=3,label=name)
        a.set(title='(a) Median and interquartile range',ylabel='Damping ratio / %'); a.legend(frameon=False,loc='upper right')
        if bins.RDT_n.iloc[-1] == 1:
            a.annotate('n = 1',(xx[-1]+.085,bins.RDT_median.iloc[-1]),xytext=(-30,20),textcoords='offset points',arrowprops={'arrowstyle':'-'})
        b.set(title='(b) Fundamental-group and passing windows',xlabel='Wind-speed bin / m s$^{-1}$',ylabel='Window count',yscale='log',ylim=(.8,max(1,bins.fundamental_windows.max())*2.5));b.legend(frameon=False,loc='upper right')
        for ax in [a,b]:ax.set_xticks(xx,['3–4','4–5','5–6','6–7','7–8','≥8']);ax.set_xlim(-.5,5.5)
    elif number==23:
        paired=d[d.EFDD_gate&d.RDT_gate]
        fig,aa=axes(size=(11,4.5));ev=paired.EFDD.to_numpy()*100;rv=paired.RDT.to_numpy()*100;mx=max(6.2,ev.max()*1.02,rv.max()*1.02) if len(ev) else 6.2
        aa[0].scatter(ev,rv,s=12,color=BLUE,alpha=.53);aa[0].plot([0,mx],[0,mx],'--',c='gray',label='Equal damping');aa[0].set(xlim=(0,mx),ylim=(0,mx),xlabel='EFDD damping ratio / %',ylabel='RDT damping ratio / %',title=f'(a) {len(paired):,} paired passing windows');aa[0].legend(frameon=False,loc='lower right')
        for name,v,c in [('EFDD',ev,BLUE),('RDT',rv,ORANGE)]:
            if not len(v):continue
            vv=np.sort(v);aa[1].step(np.r_[0,vv,mx],np.r_[0,np.arange(1,len(vv)+1)/len(vv)*100,100],where='post',c=c,label=f'{name}: median {np.median(vv):.3f}%')
            md=np.median(vv); aa[1].vlines(md,0,50,color=c,linestyles=':'); aa[1].plot(md,50,'o',c=c,ms=3)
        aa[1].axhline(50,c='gray',ls='--');aa[1].set(xlabel='Damping ratio / %',ylabel='Cumulative windows / %',ylim=(0,103),xlim=(0,mx),title='(b) Empirical cumulative distributions');aa[1].legend(frameon=False,loc='lower right')
        paired.to_csv(out/'figure23_pairs.csv',index=False)
    else: raise ValueError(number)
    for a in fig.axes:
        if a.axison:a.grid(alpha=.16);a.set_axisbelow(True)
    return fig


def export(fig, number, out):
    out=Path(out).resolve(); target=out/'figures'; target.mkdir(exist_ok=True)
    png=target/f'fig{number:02d}.png'; pdf=target/f'fig{number:02d}.pdf'
    fig.savefig(png,dpi=220,facecolor='white'); fig.savefig(pdf,facecolor='white')
    p=out/'figure_manifest.json'
    data=json.loads(p.read_text(encoding='utf-8')) if p.exists() else {
        'manuscript':'V19.0', 'language':'en', 'figures':{}, 'excluded':EXCLUDED}
    item={'title':TITLES[number], 'png':str(png.relative_to(out)), 'pdf':str(pdf.relative_to(out)),
          'plot_sha256':digest(Path(__file__)), 'run_sha256':digest(out/'run.json'),
          'sources':fig.airmem_sources, 'output_sha256':{'png':digest(png),'pdf':digest(pdf)}}
    if number==3:
        item['coverage_source']=json.loads((out/'summary.json').read_text(encoding='utf-8'))['coverage_source']
        item['same_time_input_type_as_manuscript']=item['coverage_source']=='original time vectors'
    data['figures'][str(number)]=item;write_json(p,data)
    return png,pdf

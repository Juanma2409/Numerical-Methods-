from pathlib import Path
import io,json,zipfile
import numpy as np,pandas as pd,soundfile as sf
from scipy.signal import hilbert,find_peaks
import plotly.graph_objects as go
from IPython.display import Audio,display

FEATURES=["time_seconds","rms","rms_dbfs","peak","crest_factor","dominant_frequency_hz","spectral_centroid_hz","spectral_bandwidth_hz","spectral_rolloff_hz","spectral_flatness","spectral_flux","band_0_500_ratio","band_500_1000_ratio","band_1000_2000_ratio","band_2000_4000_ratio","band_4000_8000_ratio","band_8000_12000_ratio","spectral_power"]

def find_results_folder(base="/content"):
    p=list(Path(base).rglob("metadata_final.json"));
    if not p: raise FileNotFoundError("No se encontro metadata_final.json")
    return p[0].parent

def load_project(root):
    root=Path(root); m=json.loads((root/"metadata_final.json").read_text(encoding="utf-8"));
    events=pd.read_csv(root/"events.csv") if (root/"events.csv").exists() else pd.DataFrame();
    feat=np.memmap(root/"features.f32",dtype=np.float32,mode="r",shape=(int(m["feature_frames"]),int(m["feature_count"])))
    spec=np.memmap(root/"spectrogram_clean_magnitude.f32",dtype=np.float32,mode="r",shape=(int(m["feature_frames"]),int(m["full_spectrogram"]["frequency_bins"])))
    return root,m,events,feat,spec

def read_segment(wav,start,stop):
    info=sf.info(str(wav)); a=int(start*info.samplerate); b=int(stop*info.samplerate); x,fs=sf.read(str(wav),start=a,stop=b,always_2d=False,dtype="float32");
    if x.ndim==2: x=x.mean(axis=1)
    return x,fs

def fft_segment(x,fs,nfft=65536):
    y=np.zeros(nfft); n=min(len(x),nfft); s=max(0,(len(x)-n)//2); y[:n]=x[s:s+n]; z=np.fft.rfft(y*np.hanning(nfft)); return np.fft.rfftfreq(nfft,1/fs),np.abs(z)

def maxima(f,m,fmin=100,fmax=12000,n=10,prominence=6):
    q=(f>=fmin)&(f<=fmax); ff=f[q]; mm=m[q]; db=20*np.log10(np.maximum(mm,1e-12)); p,props=find_peaks(db,prominence=prominence);
    if len(p)==0: p=np.argsort(mm)[-n:][::-1]; pr=np.full(len(p),np.nan)
    else: o=np.argsort(mm[p])[::-1][:n]; p=p[o]; pr=props["prominences"][o]
    return pd.DataFrame({"rank":np.arange(1,len(p)+1),"frequency_hz":ff[p],"magnitude":mm[p],"magnitude_db_relative":20*np.log10(np.maximum(mm[p]/max(mm.max(),1e-12),1e-12)),"prominence_db":pr})

def show_segment(root,start,stop,fmin=100,fmax=12000,nfft=65536):
    root=Path(root); x,fs=read_segment(root/"audio_clean.wav",start,stop); display(Audio(x,rate=fs)); t=start+np.arange(len(x))/fs; env=np.abs(hilbert(x.astype(float)));
    fig=go.Figure(); fig.add_scatter(x=t,y=x,mode="lines",name="Señal"); fig.add_scatter(x=t,y=env,mode="lines",name="Envolvente"); fig.add_scatter(x=t,y=-env,mode="lines",name="-Envolvente",showlegend=False); fig.update_layout(title="Forma de onda + envolvente",xaxis_title="Tiempo (s)",yaxis_title="Amplitud"); fig.show()
    f,m=fft_segment(x,fs,nfft); tab=maxima(f,m,fmin,min(fmax,fs/2)); q=(f>=fmin)&(f<=min(fmax,fs/2)); fig=go.Figure(); fig.add_scatter(x=f[q],y=20*np.log10(np.maximum(m[q],1e-12)),mode="lines",name="FFT"); fig.add_scatter(x=tab.frequency_hz,y=20*np.log10(np.maximum(tab.magnitude,1e-12)),mode="markers+text",text=[f"{v:.1f}" for v in tab.frequency_hz],textposition="top center",name="Máximos"); fig.update_layout(title="FFT + máximos representativos",xaxis_title="Frecuencia (Hz)",yaxis_title="Magnitud (dB relativa)"); fig.show(); display(tab); tab.to_csv(root/f"maximos_{start:.2f}_{stop:.2f}.csv",index=False)
    m=json.loads((root/"metadata_final.json").read_text(encoding="utf-8")); spec=np.memmap(root/"spectrogram_clean_magnitude.f32",dtype=np.float32,mode="r",shape=(int(m["feature_frames"]),int(m["full_spectrogram"]["frequency_bins"]))); hop=int(m["hop"]); a=max(0,int(start*fs/hop)); b=min(len(spec),int(stop*fs/hop)+1); z=spec[a:b]; freq=np.arange(z.shape[1])*fs/int(m["nfft"]); q=(freq>=fmin)&(freq<=min(fmax,fs/2)); freq=freq[q]; z=z[:,q]; tt=np.arange(a,b)*hop/fs; zd=20*np.log10(np.maximum(z,1e-12)); zd-=zd.max(); fig=go.Figure(go.Heatmap(x=tt,y=freq,z=zd.T,colorbar_title="dB rel.")); fig.update_layout(title="Sonograma — C++",xaxis_title="Tiempo (s)",yaxis_title="Frecuencia (Hz)"); fig.show(); return {"audio":x,"sample_rate":fs,"maxima":tab,"frequency":f,"magnitude":m}

def show_event(root,events,event_id,**kwargs):
    r=events[events.event_id==event_id].iloc[0]; return show_segment(root,float(r.start_time_s),float(r.end_time_s),**kwargs)

def compare_events(root,events,a,b,fmax=12000):
    ra=events[events.event_id==a].iloc[0]; rb=events[events.event_id==b].iloc[0]; xa,fa=read_segment(Path(root)/"audio_clean.wav",ra.start_time_s,ra.end_time_s); xb,fb=read_segment(Path(root)/"audio_clean.wav",rb.start_time_s,rb.end_time_s); A,MA=fft_segment(xa,fa); B,MB=fft_segment(xb,fb); qa=A<=min(fmax,fa/2); qb=B<=min(fmax,fb/2); fig=go.Figure(); fig.add_scatter(x=A[qa],y=20*np.log10(np.maximum(MA[qa],1e-12)),mode="lines",name=f"Evento {a}"); fig.add_scatter(x=B[qb],y=20*np.log10(np.maximum(MB[qb],1e-12)),mode="lines",name=f"Evento {b}"); fig.update_layout(title=f"Comparación {a} vs {b}",xaxis_title="Frecuencia (Hz)",yaxis_title="Magnitud (dB relativa)"); fig.show(); return fig

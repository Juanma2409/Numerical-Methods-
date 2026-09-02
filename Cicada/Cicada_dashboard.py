
from pathlib import Path
import io, json
import numpy as np
import pandas as pd
import soundfile as sf
from scipy.signal import hilbert, find_peaks
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

FEATURES = ["time_seconds","rms","rms_dbfs","peak","crest_factor","dominant_frequency_hz","spectral_centroid_hz","spectral_bandwidth_hz","spectral_rolloff_hz","spectral_flatness","spectral_flux","band_0_500_ratio","band_500_1000_ratio","band_1000_2000_ratio","band_2000_4000_ratio","band_4000_8000_ratio","band_8000_12000_ratio","spectral_power"]

def read_json(p): return json.loads(Path(p).read_text(encoding="utf-8"))

@st.cache_data
def load_project(root):
    root=Path(root); meta=read_json(root/"metadata_final.json")
    events=pd.read_csv(root/"events.csv") if (root/"events.csv").exists() else pd.DataFrame()
    return root,meta,events

@st.cache_resource
def open_memmaps(root):
    root=Path(root); m=read_json(root/"metadata_final.json")
    nf=int(m["feature_frames"]); nb=int(m["feature_count"]); sb=int(m["full_spectrogram"]["frequency_bins"])
    feat=np.memmap(root/"features.f32",dtype=np.float32,mode="r",shape=(nf,nb))
    spec=np.memmap(root/"spectrogram_clean_magnitude.f32",dtype=np.float32,mode="r",shape=(nf,sb))
    return feat,spec

@st.cache_data
def load_global(root):
    root=Path(root); m=read_json(root/"metadata_final.json"); nt=int(m["global_spectrogram"]["time_bins"]); nf=int(m["global_spectrogram"]["frequency_bins"])
    a=np.fromfile(root/"spectrogram_global.f32",dtype=np.float32).reshape(nt,nf); return a

def audio_slice(path,a,b):
    info=sf.info(str(path)); i=max(0,int(a*info.samplerate)); j=min(info.frames,int(b*info.samplerate))
    x,fs=sf.read(str(path),start=i,stop=j,always_2d=False,dtype="float32")
    if x.ndim==2: x=x.mean(axis=1).astype(np.float32)
    return x,fs

def audio_bytes(x,fs):
    buf=io.BytesIO(); sf.write(buf,x,fs,format="WAV",subtype="PCM_16"); return buf.getvalue()

def fft_segment(x,fs,nfft=65536):
    y=np.zeros(nfft); x=np.asarray(x,dtype=float); n=min(len(x),nfft); start=max(0,(len(x)-n)//2); y[:n]=x[start:start+n]
    s=np.fft.rfft(y*np.hanning(nfft)); return np.fft.rfftfreq(nfft,1/fs),np.abs(s)

def maxima_table(f,m,fmin,fmax,n=10,prom=6):
    mask=(f>=fmin)&(f<=fmax); ff=f[mask]; mm=m[mask]
    if len(ff)<3: return pd.DataFrame()
    db=20*np.log10(np.maximum(mm,1e-12)); pk,pr=find_peaks(db,prominence=prom)
    if len(pk)==0: pk=np.argsort(mm)[-n:][::-1]; prominence=np.full(len(pk),np.nan)
    else:
        order=np.argsort(mm[pk])[::-1][:n]; pk=pk[order]; prominence=pr["prominences"][order]
    ref=max(float(mm.max()),1e-12)
    return pd.DataFrame({"rank":np.arange(1,len(pk)+1),"frequency_hz":ff[pk],"magnitude":mm[pk],"magnitude_db_relative":20*np.log10(np.maximum(mm[pk]/ref,1e-12)),"prominence_db":prominence})

def spectrogram_from_cpp(root,meta,a,b,fmin,fmax,max_t=2500,max_f=700):
    _,spec=open_memmaps(str(root)); fs=float(meta["input_sample_rate"]); hop=int(meta["hop"]); nfft=int(meta["nfft"])
    i=max(0,int(np.floor(a*fs/hop))); j=min(spec.shape[0],int(np.ceil(b*fs/hop))+1); z=spec[i:j]
    freq=np.arange(z.shape[1])*fs/nfft; fm=(freq>=fmin)&(freq<=min(fmax,fs/2)); freq=freq[fm]; z=z[:,fm]; t=np.arange(i,j)*hop/fs
    if len(t)>max_t:
        q=np.linspace(0,len(t)-1,max_t,dtype=int); t=t[q]; z=z[q]
    if len(freq)>max_f:
        q=np.linspace(0,len(freq)-1,max_f,dtype=int); freq=freq[q]; z=z[:,q]
    z=20*np.log10(np.maximum(z,1e-12)); z-=np.nanmax(z)
    return go.Figure(go.Heatmap(x=t,y=freq,z=z.T,colorbar_title="dB rel."))

def spectrogram_3d(root,meta,a,b,fmin,fmax):
    _,spec=open_memmaps(str(root)); fs=float(meta["input_sample_rate"]); hop=int(meta["hop"]); nfft=int(meta["nfft"])
    i=max(0,int(a*fs/hop)); j=min(spec.shape[0],int(np.ceil(b*fs/hop))+1); z=spec[i:j]
    freq=np.arange(z.shape[1])*fs/nfft; fm=(freq>=fmin)&(freq<=min(fmax,fs/2)); freq=freq[fm]; z=z[:,fm]; t=np.arange(i,j)*hop/fs
    if len(t)>180: q=np.linspace(0,len(t)-1,180,dtype=int); t=t[q]; z=z[q]
    if len(freq)>220: q=np.linspace(0,len(freq)-1,220,dtype=int); freq=freq[q]; z=z[:,q]
    z=20*np.log10(np.maximum(z,1e-12)); z-=np.nanmax(z)
    return go.Figure(go.Surface(x=t,y=freq,z=z.T))

def waveform_fig(x,fs,a):
    env=np.abs(hilbert(x.astype(float))); q=np.linspace(0,len(x)-1,min(len(x),20000),dtype=int); t=a+q/fs
    fig=go.Figure(); fig.add_trace(go.Scattergl(x=t,y=x[q],mode="lines",name="Señal")); fig.add_trace(go.Scattergl(x=t,y=env[q],mode="lines",name="Envolvente")); fig.add_trace(go.Scattergl(x=t,y=-env[q],mode="lines",name="-Envolvente",showlegend=False)); fig.update_layout(title="Forma de onda + envolvente",xaxis_title="Tiempo (s)",yaxis_title="Amplitud digital",height=430); return fig

def fft_fig(f,m,fmax,tab):
    q=f<=fmax; fig=go.Figure(go.Scattergl(x=f[q],y=20*np.log10(np.maximum(m[q],1e-12)),mode="lines",name="FFT"))
    if not tab.empty:
        inds=[int(np.argmin(np.abs(f-v))) for v in tab.frequency_hz]; fig.add_trace(go.Scatter(x=f[inds],y=20*np.log10(np.maximum(m[inds],1e-12)),mode="markers+text",text=[f"{v:.1f}" for v in tab.frequency_hz],textposition="top center",name="Máximos"))
    fig.update_layout(title="FFT y máximos representativos",xaxis_title="Frecuencia (Hz)",yaxis_title="Magnitud (dB relativa)",height=430); return fig

def segment_block(root,meta,a,b,fmin,fmax,nfft,prom,peaks):
    audio=root/"audio_clean.wav"; x,fs=audio_slice(audio,a,b); f,m=fft_segment(x,fs,nfft); tab=maxima_table(f,m,fmin,min(fmax,fs/2),peaks,prom); return x,fs,f,m,tab

def main():
    st.set_page_config(page_title="Cicada Analyzer",page_icon="🦗",layout="wide")
    st.title("CICADA ANALYZER")
    st.caption("Interfaz científica: C++ genera los datos; Python los explora, visualiza y exporta.")
    with st.sidebar:
        root_text=st.text_input("Carpeta results", "results")
        fmin=st.number_input("Frec. mínima (Hz)",0.0,24000.0,100.0,50.0)
        fmax=st.number_input("Frec. máxima (Hz)",100.0,24000.0,12000.0,100.0)
        nfft=st.selectbox("NFFT para FFT del segmento",[8192,16384,32768,65536],index=3)
        npeaks=st.slider("Máximos representativos",3,20,10)
        prom=st.slider("Prominencia mínima (dB)",0.0,30.0,6.0,1.0)
    try: root,meta,events=load_project(root_text)
    except Exception as e: st.error(f"No se pudo cargar el proyecto: {e}"); return
    audio=root/"audio_clean.wav"; info=sf.info(str(audio)); dur=info.frames/info.samplerate
    c=st.columns(5); c[0].metric("Duración",f"{dur/60:.2f} min"); c[1].metric("Fs",f"{info.samplerate:,} Hz"); c[2].metric("Frames STFT",f"{meta['feature_frames']:,}"); c[3].metric("Eventos",str(len(events))); c[4].metric("Δf",f"{meta['frequency_resolution_hz']:.2f} Hz")
    tabs=st.tabs(["Resumen","Segmento","Eventos","Comparación A/B","Estadística","Datos"] )
    with tabs[0]:
        g=load_global(root_text); nt,nf=g.shape; tt=np.linspace(0,dur,nt); ff=np.linspace(0,meta["analysis_fmax_hz"],nf); q=ff<=fmax; z=20*np.log10(np.maximum(g[:,q],1e-12)); z-=np.nanmax(z); fig=go.Figure(go.Heatmap(x=tt,y=ff[q],z=z.T,colorbar_title="dB rel.")); fig.update_layout(title="Espectrograma global multirresolución",xaxis_title="Tiempo (s)",yaxis_title="Frecuencia (Hz)",height=560); st.plotly_chart(fig,width="stretch")
        feat,_=open_memmaps(root_text); df=pd.DataFrame(feat,columns=FEATURES); var=st.selectbox("Característica",FEATURES[1:]); tf=go.Figure(go.Scattergl(x=df.time_seconds,y=df[var],mode="lines",name=var)); tf.update_layout(title=f"Evolución temporal — {var}",height=400); st.plotly_chart(tf,width="stretch")
    with tabs[1]:
        s=st.number_input("Inicio (s)",0.0,max(0.0,dur-0.01),0.0,0.1); e=st.number_input("Fin (s)",0.01,dur,min(5.0,dur),0.1)
        if e<=s: st.error("El fin debe ser mayor que el inicio."); return
        x,fs,f,m,tab=segment_block(root,meta,s,e,fmin,fmax,nfft,prom,npeaks); st.audio(audio_bytes(x,fs),format="audio/wav"); st.plotly_chart(waveform_fig(x,fs,s),width="stretch"); st.plotly_chart(spectrogram_from_cpp(root,meta,s,e,fmin,fmax),width="stretch"); st.plotly_chart(fft_fig(f,m,min(fmax,fs/2),tab),width="stretch"); st.dataframe(tab,width="stretch",hide_index=True); st.download_button("Guardar máximos (CSV)",tab.to_csv(index=False).encode(),f"maximos_{s:.2f}_{e:.2f}.csv","text/csv"); st.plotly_chart(spectrogram_3d(root,meta,s,e,fmin,fmax),width="stretch")
    with tabs[2]:
        if events.empty: st.info("No hay eventos detectados.")
        else:
            eid=st.selectbox("Evento",events.event_id.astype(int).tolist()); r=events[events.event_id==eid].iloc[0]; a=float(r.start_time_s); b=float(r.end_time_s); x,fs,f,m,tab=segment_block(root,meta,a,b,fmin,fmax,nfft,prom,npeaks); st.subheader(f"Evento {eid}"); st.dataframe(pd.DataFrame([r]),width="stretch",hide_index=True); st.audio(audio_bytes(x,fs),format="audio/wav"); st.plotly_chart(waveform_fig(x,fs,a),width="stretch"); st.plotly_chart(spectrogram_from_cpp(root,meta,a,b,fmin,fmax),width="stretch"); st.plotly_chart(fft_fig(f,m,min(fmax,fs/2),tab),width="stretch"); st.dataframe(tab,width="stretch",hide_index=True); st.plotly_chart(spectrogram_3d(root,meta,a,b,fmin,fmax),width="stretch")
    with tabs[3]:
        if len(events)<2: st.info("Se necesitan al menos dos eventos.")
        else:
            a_id=st.selectbox("Evento A",events.event_id.astype(int).tolist(),index=0); b_id=st.selectbox("Evento B",events.event_id.astype(int).tolist(),index=1); ar=events[events.event_id==a_id].iloc[0]; br=events[events.event_id==b_id].iloc[0]; xa,fa=audio_slice(audio,float(ar.start_time_s),float(ar.end_time_s)); xb,fb=audio_slice(audio,float(br.start_time_s),float(br.end_time_s)); A,M=fft_segment(xa,fa,nfft); B,N=fft_segment(xb,fb,nfft); q1=A<=min(fmax,fa/2); q2=B<=min(fmax,fb/2); fig=go.Figure(); fig.add_trace(go.Scattergl(x=A[q1],y=20*np.log10(np.maximum(M[q1],1e-12)),name=f"A {a_id}")); fig.add_trace(go.Scattergl(x=B[q2],y=20*np.log10(np.maximum(N[q2],1e-12)),name=f"B {b_id}")); fig.update_layout(title="Espectros superpuestos A/B",xaxis_title="Frecuencia (Hz)",yaxis_title="Magnitud (dB relativa)",height=500); st.plotly_chart(fig,width="stretch"); cmp=pd.DataFrame({"metric":["duration_s","rms_dbfs","peak","dominant_frequency_hz","spectral_centroid_hz","spectral_bandwidth_hz","spectral_rolloff_hz","spectral_flatness","spectral_flux"],"A":[ar.get(k,np.nan) for k in ["duration_s","rms_dbfs","peak","dominant_frequency_hz","spectral_centroid_hz","spectral_bandwidth_hz","spectral_rolloff_hz","spectral_flatness","spectral_flux"]],"B":[br.get(k,np.nan) for k in ["duration_s","rms_dbfs","peak","dominant_frequency_hz","spectral_centroid_hz","spectral_bandwidth_hz","spectral_rolloff_hz","spectral_flatness","spectral_flux"]]}); st.dataframe(cmp,width="stretch",hide_index=True); st.download_button("Guardar comparación A/B",cmp.to_csv(index=False).encode(),f"comparacion_{a_id}_{b_id}.csv","text/csv")
    with tabs[4]:
        if events.empty: st.info("No hay eventos.")
        else:
            cols=[c for c in ["duration_s","rms_dbfs","peak","dominant_frequency_hz","spectral_centroid_hz","spectral_bandwidth_hz","spectral_rolloff_hz","spectral_flatness","spectral_flux","band_4000_8000_ratio","attack_s","decay_s"] if c in events.columns]; v=st.selectbox("Variable",cols); s=pd.to_numeric(events[v],errors="coerce").dropna(); st.write(s.describe().to_frame("value")); h=go.Figure(go.Histogram(x=s)); h.update_layout(title=f"Distribución — {v}",height=400); st.plotly_chart(h,width="stretch"); box=go.Figure(go.Box(y=s,boxpoints="outliers")); box.update_layout(title=f"Boxplot — {v}",height=400); st.plotly_chart(box,width="stretch")
    with tabs[5]:
        st.dataframe(events,width="stretch",height=400,hide_index=True); feat,_=open_memmaps(root_text); rows=st.slider("Filas de features",10,min(5000,len(feat)),min(500,len(feat)),10); st.dataframe(pd.DataFrame(feat[:rows],columns=FEATURES),width="stretch",height=400); st.download_button("Exportar eventos CSV",events.to_csv(index=False).encode(),"events_export.csv","text/csv"); st.download_button("Exportar features CSV",pd.DataFrame(feat,columns=FEATURES).to_csv(index=False).encode(),"features_export.csv","text/csv"); st.json(meta)

if __name__=="__main__": main()

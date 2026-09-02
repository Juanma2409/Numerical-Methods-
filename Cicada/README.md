# CICADA ANALYZER

**Cicada Analyzer** is a scientific bioacoustic analysis platform for long recordings of cicadas (Cicadidae), designed around a two-layer architecture:

- **C++** performs the computationally intensive digital signal processing.
- **Python** provides interactive visualization, exploratory analysis, event comparison, statistics, and export.

The project is designed for academic work in **numerical methods, computational physics, digital signal processing, and bioacoustics**, with special attention to long recordings, limited RAM, reproducibility, and transparent scientific interpretation.

The central design principle is simple:

```text
                   CICADA.WAV
                       │
                       ▼
              ┌─────────────────┐
              │   C++ ENGINE     │
              │                  │
              │ WAV → Mono       │
              │ DC removal       │
              │ STFT / FFT       │
              │ Noise reduction  │
              │ Features         │
              │ Event detection  │
              └────────┬─────────┘
                       │
          ┌────────────┼────────────┐
          ▼            ▼            ▼
      Audio clean   Features    Spectrograms
          │            │            │
          └────────────┼────────────┘
                       ▼
              ┌─────────────────┐
              │ PYTHON LAYER    │
              │                 │
              │ Dashboard       │
              │ Events          │
              │ A/B comparison  │
              │ Statistics      │
              │ 3D analysis     │
              │ Export          │
              └─────────────────┘
```

---

## 1. Project objectives

The goal is not merely to remove noise from an audio file. The project is intended to behave as a small scientific workstation for exploring cicada recordings.

The system can:

1. Read long WAV recordings efficiently.
2. Convert multichannel audio to mono.
3. Estimate and remove DC offset.
4. Calculate STFT/FFT representations.
5. Estimate a noise spectrum from an initial reference interval.
6. Apply Wiener-style spectral attenuation.
7. Reconstruct the cleaned signal using streaming overlap-add.
8. Preserve the original sample count and duration.
9. Calculate temporal and spectral features.
10. Detect acoustic episodes automatically.
11. Characterize individual events.
12. Explore events interactively.
13. Compare two events.
14. Analyze event populations statistically.
15. Inspect time-frequency structure in 2D and 3D.
16. Listen to selected regions instead of loading an entire recording into RAM.
17. Export derived data for reports and further analysis.

The software is deliberately designed so that **observation, signal-processing output, and biological interpretation are kept separate**.

---

# 2. Architecture

The project is divided into two principal layers.

## C++ — scientific processing engine

The C++ program is responsible for the heavy numerical work:

- WAV parsing
- channel conversion
- DC offset estimation/removal
- block-based processing
- STFT
- FFT/IFFT
- noise estimation
- Wiener-style attenuation
- overlap-add reconstruction
- temporal features
- spectral features
- band-energy calculations
- event detection
- event statistics
- reduced/global spectrogram generation
- full-resolution spectrogram storage
- validation of sample count and clipping
- machine-readable metadata

The final engine is `Cicada_final_main.cpp`.

## Python — visualization and interpretation layer

Python is responsible for:

- reading C++ outputs
- interactive exploration
- selecting time ranges
- playing selected audio
- plotting waveforms
- plotting envelopes
- plotting sonograms
- calculating high-resolution FFTs for selected windows
- identifying representative spectral maxima
- browsing detected events
- comparing events A/B
- statistical analysis
- 3D time-frequency visualizations
- tables
- CSV/JSON export
- scientific presentation

The main interactive interface is `Cicada_dashboard.py`.

A lightweight, notebook-oriented interface is provided by:

- `Cicada_Colab.py`
- `Cicada_Colab.ipynb`

---

# 3. Repository structure

A practical final repository can use a structure similar to:

```text
Cicada/
│
├── Cicada.wav
│
├── Cicada_final_main.cpp
├── Cicada.exe
│
├── Cicada_dashboard.py
├── Cicada_Colab.py
├── Cicada_Colab.ipynb
├── requirements_cicada_python.txt
│
├── results_final/
│   ├── audio_clean.wav
│   ├── noise_profile.f32
│   ├── spectrogram_clean_magnitude.f32
│   ├── spectrogram_clean_phase.f32
│   ├── spectrogram_global.f32
│   ├── features.f32
│   ├── events.csv
│   ├── event_features.csv
│   ├── metadata_final.json
│   └── validation_summary.json
│
└── README.md
```

Historical development files from earlier phases may also be kept in the repository, but the recommended user workflow should focus on the final engine and the final Python layer.

---

# 4. The C++ engine

## 4.1 Main source file

The principal C++ source is:

```text
Cicada_final_main.cpp
```

It is written for **C++17** and was designed to be friendly to Windows/MinGW environments without depending on `std::filesystem` or `std::clamp`.

The engine is explicitly designed to avoid loading the complete recording into memory.

Instead, the main signal-processing path is conceptually:

```text
WAV block
   ↓
mono conversion
   ↓
DC correction
   ↓
STFT frame
   ↓
FFT
   ↓
noise model
   ↓
Wiener mask
   ↓
clean spectrum
   ↓
IFFT
   ↓
streaming overlap-add
   ↓
clean WAV
```

At the same time, every STFT frame is used to compute feature vectors and populate the scientific output files.

---

# 5. WAV input

By default, the executable looks for:

```text
Cicada.wav
```

in the same directory as the executable.

This behavior is intentional because it avoids problems caused by relative paths when the program is launched from another working directory.

A different input file can be supplied explicitly with:

```text
Cicada.exe --input otro_audio.wav
```

The reader supports:

- PCM 8-bit
- PCM 16-bit
- PCM 24-bit
- PCM 32-bit
- IEEE float32 WAV

Multichannel signals are converted to mono by averaging the channels.

---

# 6. Default processing parameters

The final engine uses the following defaults:

| Parameter | Default |
|---|---:|
| NFFT | 2048 |
| HOP | 1024 samples |
| Block duration | 20 s |
| Noise reference | 10 s |
| Noise factor | 1.5 |
| Gain floor | 0.08 |
| Minimum frequency | 100 Hz |
| Maximum frequency | 12000 Hz |
| Global time bins | 4096 |
| Global frequency bins | 512 |
| C++ event threshold | 2.5 |
| Minimum event duration | 100 ms |
| Event merge gap | 250 ms |
| Save full STFT phase | yes |

These parameters are configurable from the command line.

---

# 7. Running `Cicada.exe`

## Basic execution

Place:

```text
Cicada.exe
Cicada.wav
```

in the same folder and run:

```powershell
.\Cicada.exe
```

The executable creates:

```text
results_final/
```

automatically.

---

## Show help

```powershell
.\Cicada.exe --help
```

---

## Specify another input

```powershell
.\Cicada.exe --input "otra_grabacion.wav"
```

---

## Change STFT resolution

For greater frequency resolution:

```powershell
.\Cicada.exe --nfft 4096 --hop 1024
```

For the default lightweight configuration:

```powershell
.\Cicada.exe --nfft 2048 --hop 1024
```

The theoretical frequency spacing is approximately:

\[
\Delta f = \frac{f_s}{NFFT}
\]

and the time step between consecutive frames is:

\[
\Delta t = \frac{HOP}{f_s}.
\]

For example, at 48 kHz and NFFT = 2048:

\[
\Delta f = \frac{48000}{2048}
        \approx 23.4375\;\text{Hz}.
\]

---

## Change noise-reduction parameters

Example:

```powershell
.\Cicada.exe --noise 10 --noise-factor 1.5 --floor 0.08
```

Where:

- `--noise` = duration of the initial noise-reference interval
- `--noise-factor` = scaling factor applied to the noise profile
- `--floor` = minimum spectral gain

A higher gain floor generally produces less aggressive attenuation.

---

## Change frequency range

```powershell
.\Cicada.exe --fmin 100 --fmax 12000
```

The upper limit is automatically restricted by the Nyquist frequency.

---

## Change global spectrogram resolution

```powershell
.\Cicada.exe --global-time-bins 4096 --global-freq-bins 512
```

The global representation is deliberately much smaller than the full STFT and is intended for fast visualization of the entire recording.

---

## Disable phase output

If full STFT phase is not required:

```powershell
.\Cicada.exe --no-phase
```

---

# 8. Noise reduction

The current pipeline uses a median spectral estimate from an initial reference interval.

Conceptually:

```text
first N seconds
       ↓
STFT
       ↓
magnitude spectra
       ↓
median over time
       ↓
noise profile
```

For each frequency bin, a Wiener-style gain is then computed from the ratio between the observed magnitude and the estimated noise level.

The implementation keeps a minimum gain floor to avoid completely suppressing spectral components.

### Important assumption

The first `noise_seconds` are assumed to be representative of background noise.

That assumption is not guaranteed to be correct in a biological recording. If cicadas are already active during the first few seconds, the noise model may partially absorb signal that should have been preserved.

Therefore, the noise-reduction stage should be interpreted as an engineering preprocessing step, not as an automatically perfect separation between “noise” and “animal sound”.

---

# 9. Streaming overlap-add reconstruction

One of the important technical issues in the development of this project was reconstruction.

A naïve block-by-block ISTFT can produce discontinuities or very large values at frame boundaries if overlap information is reset between blocks.

The final engine instead maintains continuous overlap-add state.

The reconstruction can be summarized as:

```text
IFFT frame
   ↓
window
   ↓
overlap-add buffer
   ↓
window normalization
   ↓
emit finalized samples
```

This allows the recording to be processed block-by-block without treating every block as an independent signal.

The output is also validated so that the number of samples matches the input and the cleaned signal remains within the expected digital range.

---

# 10. C++ feature extraction

The final engine generates one feature vector per STFT frame.

The feature file contains 18 float32 values per frame:

```text
time_seconds
rms
rms_dbfs
peak
crest_factor
dominant_frequency_hz
spectral_centroid_hz
spectral_bandwidth_hz
spectral_rolloff_hz
spectral_flatness
spectral_flux
band_0_500_ratio
band_500_1000_ratio
band_1000_2000_ratio
band_2000_4000_ratio
band_4000_8000_ratio
band_8000_12000_ratio
spectral_power
```

### RMS

Root mean square amplitude:

\[
RMS = \sqrt{\frac{1}{N}\sum_{n=0}^{N-1}x_n^2}.
\]

### RMS in dBFS

The project uses:

\[
L_{\mathrm{dBFS}}
=
20\log_{10}(RMS).
\]

This is a **digital full-scale reference**.

It is not calibrated sound-pressure level and therefore should not be interpreted as dB SPL.

### Peak

Maximum absolute sample amplitude in the frame.

### Crest factor

\[
CF = \frac{\text{peak}}{RMS}.
\]

It provides information about how impulsive or peak-dominated a signal is.

### Dominant frequency

The frequency bin with maximum magnitude within the analysis range.

Important:

```text
dominant frequency ≠ automatically fundamental frequency
```

The dominant component may be a harmonic or another strong spectral component.

### Spectral centroid

A weighted mean of frequency:

\[
f_c =
\frac{\sum_k f_k |X_k|}
     {\sum_k |X_k|}.
\]

### Spectral bandwidth

A measure of the spread of spectral energy around the centroid.

### Spectral roll-off

The frequency below which a chosen fraction of the spectral power is accumulated.

The default fraction in the engine is 85%.

### Spectral flatness

A measure of how noise-like versus tone-like the spectrum is.

### Spectral flux

Measures changes in the normalized spectrum from one frame to the next.

### Band-energy ratios

The engine stores relative energy in:

```text
0–500 Hz
500–1000 Hz
1000–2000 Hz
2000–4000 Hz
4000–8000 Hz
8000–12000 Hz
```

These values are useful for studying changes in spectral distribution over time.

---

# 11. Event detection

Event detection is one of the most important parts of the project.

The C++ engine contains a robust multi-feature heuristic based on:

- RMS level
- 4–8 kHz relative energy
- spectral flux
- inverse spectral flatness

Robust statistics are calculated using the median and MAD:

\[
\sigma_{\mathrm{robust}}
\approx 1.4826\,MAD.
\]

The resulting normalized components are combined into an event score.

This is intentionally a **signal-processing heuristic**.

It is not:

- species identification
- individual identification
- a biological classifier
- a machine-learning model trained on labeled cicadas

The engine therefore reports **acoustic episodes**, not biological identities.

---

# 12. Why Python has its own event detection

A practical problem appeared during development: a C++ detector can be mathematically valid but still be too restrictive for a particular recording.

For example, a cicada recording may contain meaningful activity outside a single preferred band.

Therefore the Dashboard includes an additional interactive detector that operates on the C++ feature stream.

This has two advantages:

1. The heavy signal processing does not have to be repeated.
2. Detection parameters can be changed interactively.

The Dashboard detector uses several features together rather than relying exclusively on amplitude.

It also has a fallback mechanism that ranks strong acoustic candidates when the main criterion produces no continuous events.

This means the user can still inspect the most active portions of the recording without the interface being forced into an empty event list.

The distinction between:

```text
detected acoustic event
```

and

```text
candidate for manual inspection
```

should always be preserved in scientific interpretation.

---

# 13. C++ output files

After running the executable, `results_final/` contains the following.

## `audio_clean.wav`

Cleaned mono audio reconstructed by the C++ processing pipeline.

This is the main audio file used by the Dashboard.

---

## `noise_profile.f32`

Float32 spectral noise profile.

This is a one-dimensional array with:

```text
NFFT / 2 + 1
```

frequency bins.

---

## `spectrogram_clean_magnitude.f32`

Full-resolution cleaned one-sided STFT magnitude.

Layout:

```text
[frame, frequency_bin]
```

Data type:

```text
float32
```

This is intentionally stored as raw binary because it can become very large.

Python reads it using `numpy.memmap`, so the complete matrix does not need to be copied into RAM.

---

## `spectrogram_clean_phase.f32`

Optional full-resolution STFT phase.

It is not required for normal Dashboard exploration, but it is preserved for future scientific or reconstruction-oriented work.

---

## `spectrogram_global.f32`

Reduced multiresolution spectrogram used for fast full-recording visualization.

Default shape:

```text
4096 × 512
```

This file is much smaller than the full STFT representation.

---

## `features.f32`

Float32 feature matrix.

Layout:

```text
[frame, feature]
```

The corresponding metadata describes the feature order.

Python opens it through memory mapping.

---

## `events.csv`

Events identified by the C++ event detector.

It contains:

- event ID
- start time
- end time
- duration
- amplitude features
- spectral features
- frequency slope
- band ratios
- attack
- decay
- detection score
- inter-event interval

This file can be empty if the C++ event heuristic is too restrictive for the recording.

That does not mean that the recording contains no acoustic events.

---

## `event_features.csv`

Compact mapping between event identifiers and the feature-frame range associated with each event.

---

## `metadata_final.json`

Machine-readable description of the processing.

It includes:

- sample rate
- original channel count
- bit depth
- duration
- DC offset
- NFFT
- HOP
- frequency resolution
- time resolution
- noise-reduction parameters
- feature count
- feature names
- full spectrogram dimensions
- global spectrogram dimensions
- event-detection configuration
- scientific limitations

---

## `validation_summary.json`

Validation information including:

- input sample count
- output sample count
- sample difference
- sample rate
- duration
- input peak
- output peak
- input RMS
- output RMS
- samples outside the digital range
- feature-frame count
- C++ event count
- basic validation status

This file is useful when documenting computational correctness in an academic report.

---

# 14. Python installation

The Python layer requires:

```text
numpy
pandas
scipy
soundfile
plotly
streamlit
```

Install them with:

```powershell
python -m pip install -r requirements_cicada_python.txt
```

Using:

```powershell
python -m pip
```

instead of bare:

```powershell
pip
```

is recommended on Windows because it guarantees that `pip` is associated with the active Python interpreter.

---

# 15. Running the Dashboard

The Dashboard is a Streamlit application.

Do **not** start it with:

```powershell
python Cicada_dashboard.py
```

Use:

```powershell
streamlit run Cicada_dashboard.py
```

or:

```powershell
python -m streamlit run Cicada_dashboard.py
```

The second form is especially convenient when several Python installations or virtual environments exist.

If the virtual environment is called `venv`, an explicit form is:

```powershell
.\venv\Scripts\python.exe -m streamlit run Cicada_dashboard.py
```

---

# 16. Recommended Windows environment setup

From the project directory:

```powershell
py -m venv venv
.\venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements_cicada_python.txt
```

Verify the interpreter:

```powershell
where.exe python
```

The first result should correspond to:

```text
...\Cicada\venv\Scripts\python.exe
```

Then verify pip:

```powershell
python -m pip --version
```

If PowerShell blocks activation temporarily:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

Then activate again.

---

# 17. Dashboard functionality

The Dashboard is organized into several analysis sections.

## Resumen

Provides the global view of the recording.

Typical elements include:

- duration
- sample rate
- number of STFT frames
- frequency resolution
- global sonogram
- event-detection score
- RMS timeline
- temporal feature inspection

The complete recording is displayed using the reduced global spectrogram rather than the full STFT matrix.

---

## Segmento

Allows selection of an arbitrary time interval.

The user can:

- select start time
- select end time
- listen only to that interval
- inspect waveform
- inspect envelope
- inspect sonogram
- inspect FFT
- locate representative maxima
- inspect the segment in 3D
- export peak information

This is particularly useful for fulfilling the basic experimental workflow:

```text
choose segment
      ↓
listen
      ↓
plot waveform
      ↓
plot sonogram
      ↓
calculate FFT
      ↓
identify spectral maxima
      ↓
save results
```

---

## Eventos

Allows the user to browse individual acoustic episodes.

For each event the Dashboard can show:

- start
- end
- duration
- RMS
- peak
- dominant frequency
- spectral centroid
- bandwidth
- roll-off
- flatness
- spectral flux
- frequency slope
- attack
- decay
- detection score
- audio
- waveform
- envelope
- sonogram
- FFT
- 3D representation
- representative maxima

---

# 18. Event detection controls

The Dashboard allows the user to modify detection without rerunning C++.

Typical controls include:

- sensitivity
- threshold
- minimum duration
- event merge gap

Three useful conceptual operating modes are:

```text
Sensitive
Balanced
Strict
```

### Sensitive

Useful when the first objective is to avoid missing candidate events.

### Balanced

Recommended as the starting point.

### Strict

Useful when the objective is to keep only stronger episodes.

The correct values depend on the recording.

There is no universal threshold that can be assumed to be correct for every environment, microphone, species, distance, or recording gain.

---

# 19. Comparison A/B

The Dashboard allows selection of two events:

```text
Event A
Event B
```

and provides:

- superimposed spectra
- individual sonograms
- duration comparison
- RMS comparison
- peak comparison
- dominant-frequency comparison
- centroid comparison
- bandwidth comparison
- flatness comparison
- flux comparison
- attack/decay comparison
- frequency-slope comparison
- audio playback

The purpose is to answer questions such as:

- Are the dominant frequencies similar?
- Is one event longer?
- Does one event contain greater high-frequency energy?
- Does one event show a stronger frequency slope?
- Are the envelopes similar?
- Do the spectra contain similar peaks or harmonic structure?

---

# 20. Statistical analysis

Once several events are available, the Dashboard can calculate descriptive statistics over the event population.

Available analyses include:

- descriptive statistics
- histograms
- boxplots
- correlations
- scatter plots
- correlation matrices

Variables can include:

- duration
- RMS
- peak
- dominant frequency
- centroid
- bandwidth
- roll-off
- flatness
- spectral flux
- band-energy ratios
- attack
- decay

These analyses make it possible to investigate questions such as:

### Frequency stability

Do events cluster around similar dominant frequencies?

### Duration variability

Are some events much longer or shorter than others?

### Frequency-duration relationship

Does event duration correlate with dominant frequency?

### Spectral evolution

Does dominant frequency change during an event?

### Event morphology

Do the envelopes show similar attack/decay patterns?

These are observational questions. A correlation does not automatically establish biological causation.

---

# 21. Time-frequency visualization

The project provides two complementary representations.

## Global view

The global spectrogram is reduced and intended to display the full recording quickly.

Typical axes:

```text
X = time
Y = frequency
color = relative spectral magnitude
```

## Local view

When the user selects a short event or segment, the Dashboard reads a limited time interval from the full-resolution memory-mapped STFT.

This provides substantially better detail without requiring the entire 30-minute recording to be rendered at full resolution.

---

# 22. Color scale

The time-frequency maps use a perceptually stronger multicolor scale rather than an all-blue palette.

Conceptually:

```text
low magnitude
     ↓
dark blue → cyan → white → yellow → orange → red
     ↓
high magnitude
```

The color represents relative spectral magnitude.

It should not be interpreted as an absolute acoustic level in dB SPL.

---

# 23. 3D analysis

The 3D representation is intentionally restricted to selected short windows.

The axes are:

```text
X = time
Y = frequency
Z = magnitude
```

The 3D view is designed for:

- rotating
- zooming
- inspecting time-frequency morphology
- comparing local spectral structures

It is not intended to represent an entire hour-long recording in 3D.

---

# 24. FFT and representative maxima

For selected segments, Python calculates a higher-resolution FFT independently of the global C++ representation.

This is intentional.

The C++ stage is optimized for processing the whole recording.

The Python stage can then perform a detailed FFT only for the small interval the user selected.

Representative maxima are identified with `scipy.signal.find_peaks`.

The Dashboard reports:

- frequency
- magnitude
- relative magnitude in dB
- prominence

These peaks can be used to investigate:

- dominant spectral components
- possible harmonic relationships
- tonal structure
- differences between events

The presence of a peak does not by itself prove that it is the biological fundamental frequency.

---

# 25. Envelopes, attack and decay

The Python layer can calculate an amplitude envelope from the selected waveform using the analytic signal.

Conceptually:

\[
e(t)=|x(t)+j\,\mathcal{H}\{x(t)\}|,
\]

where \(\mathcal{H}\) is the Hilbert transform.

This envelope is useful for examining:

- onset
- attack
- maximum amplitude
- decay
- temporal morphology

The C++ event summary also estimates attack and decay times using relative RMS thresholds.

These values should be interpreted as signal-processing measurements rather than universal biological definitions.

---

# 26. Scientific interpretation and limitations

The project intentionally includes several safeguards.

## dBFS is not dB SPL

The project uses digital signal amplitudes.

Therefore:

```text
dBFS ≠ dB SPL
```

A calibrated microphone, reference signal, and calibration procedure would be required to obtain physically calibrated sound-pressure levels.

---

## Dominant frequency is not automatically F0

A spectral maximum may correspond to:

- the fundamental
- a harmonic
- environmental noise
- resonance
- another strong spectral component

Therefore the project does not automatically declare every dominant frequency to be the fundamental frequency.

---

## Event detection is not species identification

The event detector answers:

```text
Where are there acoustically interesting episodes?
```

It does not answer:

```text
Which species produced this sound?
```

Species identification would require separate biological validation and, ideally, labeled training data and/or expert annotations.

---

## Noise estimation can be imperfect

The first 10 seconds are currently treated as a noise-reference region by default.

If cicadas are already vocalizing there, the noise model may suppress biologically relevant spectral content.

For future work, a user-selectable noise-only interval or adaptive background estimator would be a natural improvement.

---

## Aggressive denoising can remove signal

Noise reduction is always a trade-off.

A stronger attenuation can improve clarity while also removing weak biological components.

Therefore original and cleaned signals should be compared whenever the result is used in a scientific conclusion.

---

# 27. Validation

The C++ engine automatically validates important computational properties.

Among other quantities, it checks:

- input sample count
- output sample count
- duration
- RMS
- peak
- samples outside the valid digital range
- feature-frame count
- event count

The intended condition is:

```text
input samples = output samples
```

and:

```text
no output clipping
```

The project should be considered computationally healthy only when these checks pass.

---

# 28. Memory strategy

Long recordings can be very large.

For example, a roughly 30-minute recording at 48 kHz contains tens of millions of samples.

The project therefore avoids this pattern:

```text
load entire WAV
        ↓
large STFT matrix
        ↓
large Python arrays
```

Instead, the intended architecture is:

```text
long WAV
  ↓
small processing block
  ↓
STFT frame
  ↓
write result
  ↓
next block
```

For the Python layer:

```text
large .f32 file
       ↓
numpy.memmap
       ↓
read only selected region
       ↓
visualize
```

This separation is fundamental to the project's scalability.

---

# 29. Reproducible workflow

A recommended complete workflow is:

## Step 1 — prepare the project

```text
Cicada/
    Cicada.wav
    Cicada.exe
    Cicada_dashboard.py
    requirements_cicada_python.txt
```

## Step 2 — execute C++

```powershell
.\Cicada.exe
```

Wait for the engine to complete.

## Step 3 — verify the results

Check:

```text
results_final/
```

and inspect:

```text
metadata_final.json
validation_summary.json
```

## Step 4 — create/activate Python environment

```powershell
py -m venv venv
.\venv\Scripts\Activate.ps1
python -m pip install -r requirements_cicada_python.txt
```

## Step 5 — start the Dashboard

```powershell
python -m streamlit run Cicada_dashboard.py
```

## Step 6 — explore

Recommended order:

```text
Resumen
   ↓
Evento / Segmento
   ↓
audio + waveform
   ↓
sonogram
   ↓
FFT
   ↓
maxima
   ↓
3D
   ↓
A/B
   ↓
Statistics
```

## Step 7 — export

Save:

- event CSV
- maxima CSV
- statistical tables
- selected event WAV
- figures/HTML when required

---

# 30. Google Colab

For environments where a desktop Streamlit application is inconvenient, the project includes a notebook-oriented path:

```text
Cicada_Colab.py
Cicada_Colab.ipynb
```

The notebook provides functions for:

- locating the results directory
- loading metadata
- opening memory-mapped arrays
- reading selected audio segments
- plotting segments
- calculating FFTs
- finding maxima
- displaying events
- comparing events

A typical conceptual Colab workflow is:

```text
upload project/results
        ↓
load C++ outputs
        ↓
select segment
        ↓
listen
        ↓
waveform
        ↓
sonogram
        ↓
FFT
        ↓
maxima
        ↓
event comparison
```

The Colab layer is intended mainly for reproducible analysis and demonstration rather than serving as a replacement for the full local Dashboard.

---

# 31. Basic troubleshooting

## `pip` launcher error

If Windows reports an error similar to:

```text
Fatal error in launcher:
Unable to create process...
```

the virtual environment is often broken or points to another Python installation.

Recommended fix:

```powershell
deactivate
Remove-Item -Recurse -Force .\venv
py -m venv venv
.\venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements_cicada_python.txt
```

Then verify:

```powershell
where.exe python
python -m pip --version
```

---

## Streamlit warnings about `ScriptRunContext`

Do not start the application with:

```powershell
python Cicada_dashboard.py
```

Use:

```powershell
streamlit run Cicada_dashboard.py
```

or:

```powershell
python -m streamlit run Cicada_dashboard.py
```

---

## Dashboard cannot find results

Make sure that:

```text
metadata_final.json
```

exists in the selected results directory.

The Dashboard normally expects:

```text
results_final/
```

---

## No events appear

First verify that:

```text
features.f32
metadata_final.json
```

belong to the same C++ run.

Then adjust the event-detection controls.

Important: the C++ `events.csv` and the Dashboard's interactive event detector are separate mechanisms. A zero-event C++ output does not necessarily mean that Python cannot identify useful candidate episodes.

---

## Dashboard becomes slow

Avoid displaying extremely long, full-resolution windows.

Use:

```text
short segment → detailed analysis
```

rather than:

```text
30-minute full-resolution 3D visualization
```

The architecture is designed around this principle.

---

# 32. Current scientific scope

The current implementation is strongest for:

- exploratory bioacoustic analysis
- long-recording inspection
- signal-processing coursework
- numerical methods demonstrations
- FFT/STFT experiments
- feature extraction
- event detection research prototypes
- comparative event analysis

It should not yet be described as:

- a validated species classifier
- an individual-identification system
- an ecological monitoring product
- a calibrated acoustic measurement instrument

Those capabilities would require additional validation and datasets.

---

# 33. Possible future improvements

The current architecture makes several extensions possible without changing the basic C++/Python separation.

Potential future developments include:

### Better event detection

- adaptive background estimation
- user-defined noise-only interval
- hysteresis thresholds
- spectral morphology
- band-pass activity models
- annotation-assisted validation
- precision/recall evaluation against manually labeled events

### Fundamental-frequency analysis

A future Python module could estimate F0 using:

- harmonic spacing
- autocorrelation
- cepstral methods
- harmonic product spectrum
- probabilistic pitch tracking

This should remain separate from the simple dominant-frequency measure.

### Automatic event reports

A future reporting system could generate:

```text
recording summary
       ↓
event table
       ↓
selected representative events
       ↓
figures
       ↓
statistics
       ↓
PDF/HTML report
```

### Multiple recordings

The project can eventually be extended to:

```text
recording
individual
date
location
temperature
time
metadata
```

which would allow population-level comparisons.

### Clustering

Events could eventually be grouped by acoustic similarity using:

- PCA
- hierarchical clustering
- k-means
- density-based clustering

Any cluster interpretation would still require biological validation.

---

# 34. Academic reproducibility

For an academic submission, it is recommended to preserve:

1. The original WAV.
2. The exact C++ source.
3. The exact command used to execute the C++ engine.
4. `metadata_final.json`.
5. `validation_summary.json`.
6. The Python source.
7. The Python dependency file.
8. The event CSV.
9. Representative exported figures.
10. A short description of the parameter values used.

This allows another person to reproduce the processing pipeline rather than only seeing the final graphs.

---

# 35. Example command log

A reproducible experiment can be documented like this:

```powershell
.\Cicada.exe --nfft 2048 --hop 1024 --noise 10 --noise-factor 1.5 --floor 0.08 --fmin 100 --fmax 12000
```

Then:

```powershell
python -m streamlit run Cicada_dashboard.py
```

The exact parameters used should be retained with the experiment.

---

# 36. Design philosophy

The project follows several principles.

### Do the expensive processing once

C++ performs the complete recording analysis.

### Do not load unnecessary data

Large arrays remain in binary storage and are accessed selectively.

### Keep raw and derived information distinguishable

The original recording should never be confused with processed data.

### Prefer inspectable heuristics

Event detection parameters can be changed and visualized.

### Do not overinterpret a single metric

A dominant frequency, RMS value, or event score is not a complete biological explanation.

### Make every result auditable

Metadata, validation information, and exported tables should accompany scientific conclusions.

---

# 37. Citation and attribution

When using this project in an academic context, cite the software repository and describe:

- the C++ preprocessing pipeline
- the STFT configuration
- the noise-reduction method
- the feature definitions
- the event-detection criterion
- the Python analysis layer
- the parameters used for the reported experiment

Do not report generated event labels as validated species annotations unless independent biological evidence supports that conclusion.

---

# 38. Summary

The complete workflow can be reduced to:

```text
                  ORIGINAL RECORDING
                          │
                          ▼
                 ┌─────────────────┐
                 │    C++ ENGINE   │
                 │                 │
                 │ WAV             │
                 │ Mono            │
                 │ DC removal      │
                 │ STFT / FFT      │
                 │ Denoising       │
                 │ Features        │
                 │ Events          │
                 │ Validation      │
                 └────────┬────────┘
                          │
             ┌────────────┼────────────┐
             │            │            │
             ▼            ▼            ▼
        clean WAV     features      spectra
             │            │            │
             └────────────┼────────────┘
                          ▼
                 ┌─────────────────┐
                 │ PYTHON DASHBOARD│
                 │                 │
                 │ Explore         │
                 │ Listen          │
                 │ Sonogram        │
                 │ FFT             │
                 │ Maxima          │
                 │ Events          │
                 │ Compare A/B     │
                 │ Statistics      │
                 │ 3D              │
                 │ Export          │
                 └─────────────────┘
```

The result is a modular, resource-conscious, scientifically interpretable workflow for long cicada recordings, combining efficient C++ numerical processing with an interactive Python analysis environment.

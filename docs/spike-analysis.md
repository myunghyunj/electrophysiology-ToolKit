# Spike Analysis Algorithm Documentation

This document describes the spike (tic event) detection algorithm implemented in this toolkit. The algorithm is used in two modules:

- **Module 1** (`eegSpectogram.m`) - Full detection with LF validation, SNR, and spectrogram visualization
- **Module 2** (`tFUS_EventAnalyzer.m:detectPeaksFromFile()`) - Lightweight detection with z-score output for temporal distribution analysis

---

## Algorithm Overview

The detection pipeline is a **two-pass** system:

```
Raw EEG (1 kHz)
  --> Bandpass filter (1-100 Hz)
  --> Moving RMS (median-based, 200 ms window)
  --> Adaptive threshold (3.5x median RMS)
  --> Pass 1: Duration filter (>= 50 ms)
  --> Pass 2: Refractory merge (50 ms inter-peak)
  --> [Module 1 only] LF power validation (0.5-4 Hz coincidence)
  --> Output: event list with times, durations, TP/FP labels, SNR
```

---

## Tunable Parameters

All parameters are defined in the `params` struct at the top of each script.

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `win_sec` | 0.2 | s | RMS sliding window length |
| `k_thresh` | 3.5 | dimensionless | Threshold multiplier (baseline = median RMS) |
| `min_dur_ms` | 50 | ms | Minimum event duration to accept |
| `refractory_ms` | 50 | ms | Minimum inter-peak interval; closer peaks are merged |
| `LF_band` | [0.5, 4] | Hz | Low-frequency band for true-positive validation |
| `LF_win_sec` | 0.25 | s | Window for LF power envelope calculation |
| `k_LF_thresh` | 3 | dimensionless | LF power threshold multiplier (baseline = median LF power) |
| `tp_window_ms` | 50 | ms | Coincidence window between EEG peak and LF burst |

---

## Step-by-Step Algorithm

### 1. Preprocessing: Bandpass Filter

**Source**: `eegSpectogram.m:116-118`, `tFUS_EventAnalyzer.m:231-232`

```matlab
[b, a] = butter(4, [1 100]/(Fs/2), 'bandpass');
eeg_filtered = filtfilt(b, a, eeg);
```

- **Filter**: 4th-order Butterworth bandpass (1-100 Hz)
- **Application**: Zero-phase filtering via `filtfilt` (no phase distortion)
- **Nyquist**: 500 Hz (Fs = 1000 Hz)
- **Purpose**: Removes DC offset, slow drift (<1 Hz), and high-frequency noise (>100 Hz) while preserving the tic waveform morphology

### 2. Moving RMS Calculation (Median-Based)

**Source**: `eegSpectogram.m:128-133`, `tFUS_EventAnalyzer.m:235-240`

```matlab
win_samples = round(0.2 * 1000);  % 200 samples
eeg_sq = eeg_filtered.^2;
rms_values = sqrt(movmedian(eeg_sq, win_samples, 'omitnan'));
```

- **Window**: 200 ms (200 samples at 1 kHz), forced to odd length for `movmedian` symmetry
- **Method**: `movmedian` of the squared signal, then square root
- **Why median instead of mean?**: The median is robust to outliers. A single large spike does not inflate the local RMS estimate as severely, yielding a more stable envelope that better represents the local "typical" amplitude
- **`'omitnan'`**: Ensures robustness at signal edges

### 3. Adaptive Threshold

**Source**: `eegSpectogram.m:136-137`, `tFUS_EventAnalyzer.m:243-244`

```matlab
baseline_rms = median(rms_values, 'omitnan');
threshold = 3.5 * baseline_rms;
```

- **Baseline**: Global median of the entire RMS trace
- **Threshold**: `k_thresh` (default 3.5) times the baseline
- **Adaptive**: Each recording gets its own threshold, accounting for different noise floors across animals, sessions, or electrode impedances
- **Rationale**: The 3.5x multiplier was empirically tuned for Tourette-like tic detection in preclinical rodent models. This is deliberately above the common 3x standard to reduce false positives in noisy recordings

### 4. Pass 1: Event Extraction with Duration Filter

**Source**: `eegSpectogram.m:148-180`, `tFUS_EventAnalyzer.m:250-291`

1. **Find supra-threshold points**: All indices where `rms > threshold`
2. **Form contiguous runs**: Group consecutive supra-threshold samples into runs (gap > 1 sample starts a new run)
3. **Duration filter**: Discard runs shorter than `min_dur_ms` (50 ms = 50 samples)
4. **Peak extraction**: Within each surviving run, find the sample with the maximum RMS value

**Why 50 ms minimum duration?** This filters out transient noise spikes and electrical artifacts that are too brief to represent genuine tic events. Neurophysiological tic-related EEG discharges typically last >= 50 ms.

### 5. Pass 2: Refractory Period Merge

**Source**: `eegSpectogram.m:182-211`, `tFUS_EventAnalyzer.m:293-320`

This is the core refractory period implementation:

```
For each pair of consecutive candidate events (sorted by peak time):
  inter_peak_interval = peak_time[i] - peak_time[i-1]

  If inter_peak_interval < refractory_ms (50 ms):
    MERGE the two events:
      - Extend boundaries to encompass both events
      - Keep the peak with the HIGHER RMS value
      - Discard the lower peak
  Else:
    Accept as a distinct, separate event
```

#### What the refractory period does

The **50 ms refractory period** enforces a minimum separation between detected events. When two candidate events have peaks closer than 50 ms apart, they are merged into a single event. This serves multiple purposes:

1. **Prevents double-counting**: A single tic can produce multiple RMS excursions above threshold (e.g., the initial discharge and afterwave). The refractory merge ensures this counts as one event.

2. **Biological motivation**: The 50 ms value is inspired by neuronal refractory periods. After a tic-related discharge, the neural circuit requires a recovery period before it can generate another independent event. Events arriving within this window are treated as part of the same discharge complex.

3. **Peak selection**: When merging, the algorithm retains the peak with the higher RMS value. This ensures the reported peak time corresponds to the strongest part of the discharge.

#### Implementation differences

**Module 1** (`eegSpectogram.m:182-211`) uses a **forward-scanning** merge:
- Sorts candidates by peak time
- Iterates forward, comparing each new event to the last accepted event
- Extends the merged event's end boundary to encompass both

**Module 2** (`tFUS_EventAnalyzer.m:293-320`) uses a **backward-scanning** merge:
- Starts from the last event and works backward
- Deletes the weaker of two too-close peaks from the arrays
- Uses `i = i - 2` to recheck after deletion (avoids cascading merge errors)

Both produce equivalent results for well-separated events.

### 6. Low-Frequency Power Validation (Module 1 Only)

**Source**: `eegSpectogram.m:229-300`

This step classifies each detected event as a **True Positive (TP)** or **False Positive (FP)**:

```
1. Bandpass filter the signal to 0.5-10 Hz (2nd-order Butterworth)
   - With 500-sample symmetric padding to reduce edge artifacts
2. Compute instantaneous LF power: movmean(eeg_LF^2, 250 ms window)
3. Set LF threshold: 3x median(LF_power)
4. Create binary burst mask: LF_power > LF_threshold
5. For each detected EEG event:
   - Check if ANY sample in the LF burst mask is TRUE
     within +/- 50 ms of the event's peak time
   - If yes: label as True Positive (TP)
   - If no:  label as False Positive (FP)
```

**Rationale**: Genuine tic events in EEG are accompanied by a coincident low-frequency power burst (0.5-4 Hz delta band). This reflects the slow cortical potential shift associated with motor tics. Events that exceed the RMS threshold but lack this LF signature are likely artifacts (e.g., movement artifacts, electrode noise).

### 7. SNR Calculation (Module 1 Only)

**Source**: `eegSpectogram.m:306-413`

For each event:
```matlab
SNR_dB = 20 * log10(peak_amplitude / noise_std)
```

- **Peak amplitude**: `abs(eeg_filtered(peak_sample))` - absolute voltage at the event peak
- **Noise std**: Standard deviation of a 1-second baseline window, taken from immediately before the event start (or after the event end if insufficient pre-event signal)
- **Units**: dB (voltage ratio, hence 20*log10)

### 8. Z-Score Calculation (Module 2 Only)

**Source**: `tFUS_EventAnalyzer.m:322-340`

Per-file z-scoring of peak intensities:
```matlab
peak_zscores = (peak_intensities - mean(peak_intensities)) / std(peak_intensities)
```

- Normalizes peak RMS values relative to the distribution within each file
- Enables cross-file comparison of relative event magnitude
- Used downstream for temporal distribution color-coding (plotTemporalPSDmap)

---

## Output Summary

| Output | Module 1 | Module 2 |
|--------|----------|----------|
| Event start/peak/end times | Yes | Peak times only |
| Event duration (ms) | Yes | No |
| Peak RMS value | Yes | Yes (as intensity) |
| TP/FP classification | Yes (LF validation) | No |
| SNR (dB) | Yes | No |
| Z-scores | No | Yes |
| Event rate (events/min) | Yes | Derived from temporal binning |

---

## Spectrogram Parameters (Module 1)

Used for visualization alongside detected events:

| Parameter | Value |
|-----------|-------|
| Window | Hamming, 256 ms |
| FFT size | 256 samples |
| Overlap | 90% |
| Max frequency (EEG) | 50 Hz |
| Max frequency (EMG) | 150 Hz |
| PSD floor | 1e-12 (uV)^2/Hz |
| Log scale | log10(ASD) where ASD = sqrt(PSD) |

---

## References

- Detection parameters tuned for preclinical Tourette-like tic models (Kim, Kang, Jeong et al., 2025 - paper under review)
- Intan RHD file reader based on [Intan Technologies MATLAB file readers](https://www.intantech.com/downloads.html)

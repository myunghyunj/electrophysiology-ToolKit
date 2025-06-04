# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This MATLAB project analyzes EEG peak temporal distributions from tFUS (transcranial focused ultrasound) experiments. It processes EEG recordings to detect spike events, compares STIM vs SHAM conditions, and generates temporal Power Spectral Density (PSD) maps with z-score overlays.

## Architecture

### Core Components

1. **tFUS_EventAnalyzer_v3.m** - Main analysis script
   - Processes multiple .mat files containing EEG recordings
   - Detects peaks using RMS-based spike detection via `detectPeaksFromFile()` helper function
   - Groups peaks by condition (STIM vs SHAM)
   - Key parameters: win_sec=0.2s, k_thresh=3.5, min_dur_ms=50ms, refractory_ms=50ms
   - Generates z-score distribution bar plots

2. **plotTemporalPSDmap.m** - Visualization function
   - Creates 2-panel temporal PSD maps with z-score spectrum overlay
   - Generates bar plots with KDE smoothing curves
   - Calculates Standard Error of Mean (SEM) across files
   - Exports results to Excel and saves figures (.png and .fig)
   - Supports Adobe .ase color palettes

### Data Processing Pipeline

1. **Input**: User selects STIM and SHAM folders containing .mat files
   - Expected format: `raw_data` field in each .mat file
   - First row contains EEG data sampled at 1 kHz

2. **Peak Detection** (in `detectPeaksFromFile`):
   - Bandpass filter EEG data (1-100 Hz)
   - Calculate moving RMS with adaptive threshold
   - Identify supra-threshold runs meeting duration criteria
   - Apply refractory period merging
   - Calculate z-scores relative to each file's peak distribution

3. **Visualization**:
   - Temporal binning (user-defined, default 5 min)
   - Z-score spectrum overlay with customizable increments
   - Optional Mean+SEM error bars and KDE smoothing

4. **Output**: All results saved to `analysis/output/`

### Key Modifications (per README)

- Z-score colorbar uses 2 sigma intervals (plotTemporalPSDmap.m:469)
- ASE color palette is inverted when loaded (plotTemporalPSDmap.m:265)
- Y-axis ticks dynamically adjust based on data range (plotTemporalPSDmap.m:293-298)

## Commands

### Running Analysis
```matlab
% In MATLAB command window:
tFUS_EventAnalyzer_v3
```

### User Inputs During Execution
1. **Folder Selection**:
   - Select folder containing STIM .mat files
   - Select folder containing SHAM .mat files

2. **Analysis Parameters** (via dialog):
   - Bin size in minutes (default: 5)
   
3. **Visualization Parameters** (via dialog in plotTemporalPSDmap):
   - Z-score increment for spectrum overlay (default: 0.5)
   - Display Mean+SEM error bars (yes/no, default: yes)
   - Display KDE smoothing curves (yes/no, default: yes)

### Expected Data Format
- Input files: `.mat` files with `raw_data` field
- `raw_data`: Matrix where first row is EEG data (1 kHz sampling rate)

### Output Files
All outputs are saved to: `analysis/output/`
- `Temporal_PSD_RawData.xlsx` - Raw binned data for both conditions
- `Temporal_PSDmap.fig/.png` - Main temporal PSD visualization
- `Zscore_BarPlot_Distribution.fig/.png` - Z-score distribution comparison
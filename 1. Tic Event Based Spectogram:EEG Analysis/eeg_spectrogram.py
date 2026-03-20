from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, filtfilt, medfilt, spectrogram


DEFAULT_FS = 1000


@dataclass
class AnalysisParams:
    win_sec: float = 0.2
    k_thresh: float = 3.5
    min_dur_ms: float = 50
    refractory_ms: float = 50
    lf_band: tuple[float, float] = (0.5, 10.0)
    lf_win_sec: float = 0.25
    k_lf_thresh: float = 3.0
    tp_window_ms: float = 50


@dataclass
class EEGEvent:
    start_time: float
    peak_time: float
    end_time: float
    duration_ms: float
    peak_rms_value: float
    is_tp: bool
    snr_db: float


@dataclass
class AnalysisResult:
    file: Path
    events: list[EEGEvent]
    threshold_eeg: float
    baseline_eeg_rms: float
    lf_threshold: float
    lf_baseline: float
    has_emg: bool


def load_raw_data(mat_file: Path) -> np.ndarray:
    try:
        data = loadmat(mat_file, squeeze_me=True, struct_as_record=False)
        for key in ('raw_data', 'data'):
            if key in data:
                return ensure_2d_signal(np.asarray(data[key], dtype=float))
        for key, value in data.items():
            if key.startswith('__'):
                continue
            arr = np.asarray(value)
            if arr.ndim >= 1 and np.issubdtype(arr.dtype, np.number):
                return ensure_2d_signal(arr.astype(float))
    except NotImplementedError:
        pass
    except ValueError as exc:
        if 'Unknown mat file type' not in str(exc):
            raise

    with h5py.File(mat_file, 'r') as f:
        for key in ('raw_data', 'data'):
            if key in f:
                return ensure_2d_signal(np.array(f[key], dtype=float).T)
        for key in f.keys():
            obj = f[key]
            if isinstance(obj, h5py.Dataset):
                return ensure_2d_signal(np.array(obj, dtype=float).T)
    raise ValueError(f'Could not find EEG array in {mat_file}')


def ensure_2d_signal(arr: np.ndarray) -> np.ndarray:
    arr = np.squeeze(arr)
    if arr.ndim == 1:
        return arr[np.newaxis, :]
    if arr.shape[0] > arr.shape[1]:
        return arr.T
    return arr


def moving_rms(signal: np.ndarray, win_samples: int) -> np.ndarray:
    if win_samples % 2 == 0:
        win_samples += 1
    return np.sqrt(np.maximum(medfilt(signal ** 2, kernel_size=win_samples), 0.0))


def moving_mean_power(signal: np.ndarray, win_samples: int) -> np.ndarray:
    return uniform_filter1d(signal ** 2, size=max(1, win_samples), mode='nearest')


def bandpass_filter(signal: np.ndarray, fs: int, low: float, high: float, order: int = 4) -> np.ndarray:
    nyq = fs / 2.0
    b, a = butter(order, [low / nyq, high / nyq], btype='bandpass')
    return filtfilt(b, a, signal)


def lf_bandpass_with_padding(signal: np.ndarray, fs: int, low: float, high: float, order: int = 2, pad_samples: int = 500) -> np.ndarray:
    nyq = fs / 2.0
    b, a = butter(order, [low / nyq, high / nyq], btype='bandpass')
    if signal.size <= 2 * pad_samples:
        return filtfilt(b, a, signal)
    padded = np.concatenate((signal[:pad_samples][::-1], signal, signal[-pad_samples:][::-1]))
    filtered = filtfilt(b, a, padded)
    return filtered[pad_samples:-pad_samples]


def find_candidate_events(rms_values: np.ndarray, threshold: float, min_dur_samples: int) -> list[dict[str, float]]:
    supra = np.flatnonzero(rms_values > threshold)
    if supra.size == 0:
        return []
    runs = np.split(supra, np.where(np.diff(supra) > 1)[0] + 1)
    candidates: list[dict[str, float]] = []
    for run in runs:
        if run.size < min_dur_samples:
            continue
        peak_idx = int(run[np.argmax(rms_values[run])])
        candidates.append({'start_idx': int(run[0]), 'peak_idx': peak_idx, 'end_idx': int(run[-1]), 'peak_rms_value': float(rms_values[peak_idx])})
    return candidates


def merge_events(candidates: list[dict[str, float]], refractory_samples: int) -> list[dict[str, float]]:
    if not candidates:
        return []
    ordered = sorted(candidates, key=lambda item: item['peak_idx'])
    merged = [ordered[0].copy()]
    for current in ordered[1:]:
        last = merged[-1]
        if current['peak_idx'] - last['peak_idx'] < refractory_samples:
            last['end_idx'] = max(last['end_idx'], current['end_idx'])
            if current['peak_rms_value'] > last['peak_rms_value']:
                last['peak_idx'] = current['peak_idx']
                last['peak_rms_value'] = current['peak_rms_value']
        else:
            merged.append(current.copy())
    return merged


def compute_snr_db(signal: np.ndarray, start_idx: int, peak_idx: int, end_idx: int, fs: int) -> float:
    baseline_samples = fs
    if start_idx > baseline_samples:
        noise_idx = np.arange(start_idx - baseline_samples, start_idx)
    elif end_idx + baseline_samples < signal.size:
        noise_idx = np.arange(end_idx + 1, end_idx + baseline_samples + 1)
    else:
        noise_idx = np.arange(max(0, peak_idx - baseline_samples // 2), min(signal.size, peak_idx + baseline_samples // 2))
    noise = signal[noise_idx]
    noise_std = float(np.std(noise)) if noise.size else float(np.std(signal))
    noise_std = max(noise_std, np.finfo(float).eps)
    peak_amp = abs(float(signal[min(max(peak_idx, 0), signal.size - 1)]))
    return 20.0 * math.log10(max(peak_amp, np.finfo(float).eps) / noise_std)


def compute_log_spectrogram(signal: np.ndarray, fs: int, nfft: int, noverlap: int, f_cut: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    freqs, times, psd = spectrogram(signal.astype(float), fs=fs, window='hamming', nperseg=nfft, noverlap=noverlap, nfft=nfft, scaling='density', mode='psd')
    keep = freqs <= f_cut
    psd = np.maximum(psd[keep], np.finfo(float).tiny)
    asd = np.sqrt(psd)
    return freqs[keep], times, np.log10(asd)


def analyze_file(mat_file: Path, output_dir: Path, fs: int, params: AnalysisParams, save_detection_envelope: bool = False, save_eps: bool = False) -> AnalysisResult:
    raw = load_raw_data(mat_file)
    eeg = raw[0].astype(float)
    emg = raw[1].astype(float) if raw.shape[0] >= 2 and raw.shape[1] == raw[1].size else None
    has_emg = emg is not None and np.any(np.isfinite(emg)) and np.any(emg != 0)

    eeg_plot = bandpass_filter(eeg, fs, 1, 100)
    emg_plot = bandpass_filter(emg, fs, 1, 100) if has_emg else None
    rms = moving_rms(eeg_plot, round(params.win_sec * fs))
    baseline_eeg_rms = float(np.median(rms))
    threshold_eeg = float(params.k_thresh * baseline_eeg_rms)

    candidates = find_candidate_events(rms, threshold_eeg, round(params.min_dur_ms / 1000 * fs))
    merged = merge_events(candidates, round(params.refractory_ms / 1000 * fs))

    eeg_lf = lf_bandpass_with_padding(eeg_plot, fs, params.lf_band[0], params.lf_band[1], order=2, pad_samples=500)
    lf_power = moving_mean_power(eeg_lf, round(params.lf_win_sec * fs))
    finite = lf_power[np.isfinite(lf_power)]
    lf_baseline = float(np.median(finite)) if finite.size else float('nan')
    lf_threshold = float(params.k_lf_thresh * lf_baseline) if np.isfinite(lf_baseline) else float('nan')
    lf_mask = lf_power > lf_threshold if np.isfinite(lf_threshold) else np.zeros_like(lf_power, dtype=bool)
    tp_samples = round(params.tp_window_ms / 1000 * fs)

    events: list[EEGEvent] = []
    for event in merged:
        peak_idx = event['peak_idx']
        start_idx = event['start_idx']
        end_idx = event['end_idx']
        lo = max(0, peak_idx - tp_samples)
        hi = min(lf_mask.size, peak_idx + tp_samples + 1)
        is_tp = bool(np.any(lf_mask[lo:hi]))
        events.append(EEGEvent(start_time=start_idx / fs, peak_time=peak_idx / fs, end_time=end_idx / fs, duration_ms=(end_idx - start_idx + 1) / fs * 1000.0, peak_rms_value=float(event['peak_rms_value']), is_tp=is_tp, snr_db=compute_snr_db(eeg_plot, start_idx, peak_idx, end_idx, fs)))

    save_main_figure(mat_file, output_dir, eeg, eeg_plot, emg, emg_plot, events, fs, save_eps=save_eps)
    if save_detection_envelope:
        save_detection_figure(mat_file, output_dir, rms, threshold_eeg, fs, save_eps=save_eps)
    save_event_table(mat_file, output_dir, events)

    return AnalysisResult(file=mat_file, events=events, threshold_eeg=threshold_eeg, baseline_eeg_rms=baseline_eeg_rms, lf_threshold=lf_threshold, lf_baseline=lf_baseline, has_emg=has_emg)


def save_event_table(mat_file: Path, output_dir: Path, events: list[EEGEvent]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f'{mat_file.stem}_events.csv'
    with csv_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(EEGEvent(0, 0, 0, 0, 0, False, 0)).keys()))
        writer.writeheader()
        for event in events:
            writer.writerow(asdict(event))


def save_main_figure(mat_file: Path, output_dir: Path, eeg: np.ndarray, eeg_plot: np.ndarray, emg: np.ndarray | None, emg_plot: np.ndarray | None, events: list[EEGEvent], fs: int, save_eps: bool = False) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    t = np.arange(eeg.size) / fs
    fig, axes = plt.subplots(2, 2, figsize=(16, 8), constrained_layout=True)
    fig.suptitle(f'EEG Event Analysis - {mat_file.name}')

    ax = axes[0, 0]
    ax.plot(t, eeg_plot, color='black', linewidth=0.8)
    if events:
        all_peaks = np.array([event.peak_time for event in events])
        ax.scatter(all_peaks, np.zeros_like(all_peaks), c='red', s=24, label='EEG event')
        tp_peaks = np.array([event.peak_time for event in events if event.is_tp])
        if tp_peaks.size:
            ax.scatter(tp_peaks, np.zeros_like(tp_peaks), facecolors='none', edgecolors='blue', s=80, label='LF validated')
        ax.legend(loc='upper right', framealpha=1.0, fancybox=False)
    ax.set_title('Filtered EEG (1-100 Hz)')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude (µV)')
    ax.grid(True)

    ax = axes[0, 1]
    if emg_plot is not None:
        ax.plot(t, emg_plot, color='gray', linewidth=0.8)
        ax.set_title('Filtered EMG (1-100 Hz)')
        ax.set_ylabel('Amplitude (µV)')
    else:
        ax.text(0.5, 0.5, 'No EMG data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Filtered EMG (No Data)')
    ax.set_xlabel('Time (s)')
    ax.grid(True)

    ax = axes[1, 0]
    freqs, times, logspec = compute_log_spectrogram(eeg, fs, nfft=256, noverlap=230, f_cut=50)
    mesh = ax.pcolormesh(times, freqs, logspec, shading='auto', cmap='jet')
    mesh.set_rasterized(True)
    fig.colorbar(mesh, ax=ax, label='log10(µV/√Hz)')
    ax.set_title('EEG Spectrogram (Unfiltered)')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Frequency (Hz)')

    ax = axes[1, 1]
    if emg is not None and emg.size >= 256:
        emg_freqs, emg_times, emg_logspec = compute_log_spectrogram(emg, fs, nfft=256, noverlap=230, f_cut=150)
        mesh_emg = ax.pcolormesh(emg_times, emg_freqs, emg_logspec, shading='auto', cmap='jet')
        mesh_emg.set_rasterized(True)
        fig.colorbar(mesh_emg, ax=ax, label='log10(µV/√Hz)')
        ax.set_title('EMG Spectrogram (Unfiltered)')
        ax.set_ylabel('Frequency (Hz)')
    else:
        ax.text(0.5, 0.5, 'No EMG data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('EMG Spectrogram (No Data)')
    ax.set_xlabel('Time (s)')

    fig.savefig(output_dir / f'{mat_file.stem}_analysis.png', dpi=200)
    if save_eps:
        fig.savefig(output_dir / f'{mat_file.stem}_analysis.eps', format='eps')
    plt.close(fig)


def save_detection_figure(mat_file: Path, output_dir: Path, rms: np.ndarray, threshold_eeg: float, fs: int, save_eps: bool = False) -> None:
    t = np.arange(rms.size) / fs
    fig, ax = plt.subplots(figsize=(12, 4), constrained_layout=True)
    ax.plot(t, rms, color='purple', linewidth=0.8, label='RMS envelope')
    ax.axhline(threshold_eeg, color='red', linestyle='--', label='RMS threshold')
    ax.set_title('Detection Envelope')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('RMS (µV)')
    ax.grid(True)
    ax.legend(loc='upper right', framealpha=1.0, fancybox=False)
    fig.savefig(output_dir / f'{mat_file.stem}_detection_envelope.png', dpi=200)
    if save_eps:
        fig.savefig(output_dir / f'{mat_file.stem}_detection_envelope.eps', format='eps')
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Python port of eegSpectogram.m')
    parser.add_argument('inputs', nargs='+', help='MAT files or directories containing MAT files')
    parser.add_argument('--output-dir', default='analyzed_figures_py', help='Directory for figures and event tables')
    parser.add_argument('--fs', type=int, default=DEFAULT_FS, help='Sampling frequency in Hz')
    parser.add_argument('--save-detection-envelope', action='store_true', help='Save detection-envelope figure in addition to MATLAB-parity 4-panel figure')
    parser.add_argument('--save-eps', action='store_true', help='Also export EPS vector outputs (slower on long recordings)')
    return parser.parse_args()


def iter_mat_files(inputs: Iterable[str]) -> list[Path]:
    files: list[Path] = []
    for raw in inputs:
        path = Path(raw)
        if path.is_dir():
            files.extend(sorted(path.glob('*.mat')))
        elif path.suffix.lower() == '.mat':
            files.append(path)
    deduped: list[Path] = []
    seen = set()
    for file in files:
        resolved = file.resolve()
        if resolved not in seen:
            seen.add(resolved)
            deduped.append(file)
    return deduped


def main() -> None:
    args = parse_args()
    params = AnalysisParams()
    files = iter_mat_files(args.inputs)
    if not files:
        raise SystemExit('No .mat files found.')
    output_dir = Path(args.output_dir)
    results = [analyze_file(file, output_dir, args.fs, params, save_detection_envelope=args.save_detection_envelope, save_eps=args.save_eps) for file in files]
    total_events = sum(len(result.events) for result in results)
    total_tp = sum(sum(event.is_tp for event in result.events) for result in results)
    print(f'Processed {len(results)} file(s).')
    print(f'Total events: {total_events}')
    print(f'True positives: {total_tp}')
    for result in results:
        print(f'- {result.file.name}: {len(result.events)} events, RMS thr={result.threshold_eeg:.2f}, LF thr={result.lf_threshold:.2f}')


if __name__ == '__main__':
    main()

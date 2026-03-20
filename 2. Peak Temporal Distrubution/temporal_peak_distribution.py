from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.signal import butter, filtfilt, medfilt
from scipy.stats import gaussian_kde

try:
    from openpyxl import Workbook
except Exception:
    Workbook = None

DEFAULT_FS = 1000


@dataclass
class DetectionParams:
    win_sec: float = 0.2
    k_thresh: float = 3.5
    min_dur_ms: float = 50
    refractory_ms: float = 50


def load_raw_data(mat_file: Path) -> np.ndarray:
    try:
        data = loadmat(mat_file, squeeze_me=True, struct_as_record=False)
        for key in ('raw_data', 'data'):
            if key in data:
                return ensure_2d(np.asarray(data[key], dtype=float))
    except NotImplementedError:
        pass
    with h5py.File(mat_file, 'r') as f:
        for key in ('raw_data', 'data'):
            if key in f:
                return ensure_2d(np.array(f[key], dtype=float).T)
    raise ValueError(f'Could not find raw_data in {mat_file}')


def ensure_2d(arr: np.ndarray) -> np.ndarray:
    arr = np.squeeze(arr)
    if arr.ndim == 1:
        return arr[np.newaxis, :]
    if arr.shape[0] > arr.shape[1]:
        return arr.T
    return arr


def bandpass(signal: np.ndarray, fs: int, low: float, high: float, order: int = 4) -> np.ndarray:
    b, a = butter(order, [low / (fs / 2), high / (fs / 2)], btype='bandpass')
    return filtfilt(b, a, signal)


def moving_rms(signal: np.ndarray, win_samples: int) -> np.ndarray:
    if win_samples % 2 == 0:
        win_samples += 1
    return np.sqrt(np.maximum(medfilt(signal ** 2, kernel_size=win_samples), 0.0))


def detect_peaks_from_file(mat_file: Path, fs: int, params: DetectionParams) -> tuple[np.ndarray, np.ndarray]:
    raw = load_raw_data(mat_file)
    eeg = raw[0].astype(float)
    filtered = bandpass(eeg, fs, 1, 100)
    rms = moving_rms(filtered, round(params.win_sec * fs))
    threshold = params.k_thresh * float(np.median(rms))
    supra = np.flatnonzero(rms > threshold)
    if supra.size == 0:
        return np.array([]), np.array([])

    min_dur = round(params.min_dur_ms / 1000 * fs)
    refractory = round(params.refractory_ms / 1000 * fs)
    runs = np.split(supra, np.where(np.diff(supra) > 1)[0] + 1)

    peak_times: list[float] = []
    peak_vals: list[float] = []
    for run in runs:
        if run.size < min_dur:
            continue
        peak_idx = int(run[np.argmax(rms[run])])
        peak_times.append(peak_idx / fs)
        peak_vals.append(float(rms[peak_idx]))

    if not peak_times:
        return np.array([]), np.array([])

    order = np.argsort(peak_times)
    peak_times = list(np.asarray(peak_times)[order])
    peak_vals = list(np.asarray(peak_vals)[order])

    merged_times = [peak_times[0]]
    merged_vals = [peak_vals[0]]
    for time_val, peak_val in zip(peak_times[1:], peak_vals[1:]):
        if (time_val - merged_times[-1]) * fs < refractory:
            if peak_val > merged_vals[-1]:
                merged_times[-1] = time_val
                merged_vals[-1] = peak_val
        else:
            merged_times.append(time_val)
            merged_vals.append(peak_val)

    vals = np.asarray(merged_vals)
    if vals.size <= 1:
        zscores = np.zeros_like(vals)
    else:
        sample_std = np.std(vals, ddof=1)
        zscores = np.zeros_like(vals) if sample_std == 0 else (vals - np.mean(vals)) / sample_std
    return np.asarray(merged_times), zscores


def collect_group(files: list[Path], fs: int, params: DetectionParams) -> tuple[np.ndarray, np.ndarray, list[str]]:
    all_times: list[float] = []
    all_z: list[float] = []
    all_names: list[str] = []
    for file in files:
        times, zscores = detect_peaks_from_file(file, fs, params)
        if times.size:
            all_times.extend(times.tolist())
            all_z.extend(zscores.tolist())
            all_names.extend([file.name] * len(times))
    if not all_times:
        return np.array([]), np.array([]), []
    order = np.argsort(all_times)
    return np.asarray(all_times)[order], np.asarray(all_z)[order], list(np.asarray(all_names)[order])


def stacked_counts(peak_times: np.ndarray, zscores: np.ndarray, filenames: list[str], rec_dur_min: int, bin_min: float, z_increment: float):
    edges = np.arange(0, rec_dur_min + bin_min, bin_min)
    centers = edges[:-1] + bin_min / 2
    z_min, z_max = -4.0, 4.0
    z_edges = np.concatenate(([-np.inf], np.arange(z_min, z_max + z_increment, z_increment), [np.inf]))
    num_z_bins = len(z_edges) - 1
    unique_files = sorted(set(filenames))
    per_file = np.zeros((len(unique_files), len(edges) - 1, num_z_bins), dtype=float)
    counts_per_file = np.zeros((len(unique_files), len(edges) - 1), dtype=float)

    if peak_times.size:
        binned = np.digitize(peak_times / 60.0, edges, right=False) - 1
        z_binned = np.digitize(zscores, z_edges, right=False) - 1
        filename_array = np.asarray(filenames)
        for file_idx, name in enumerate(unique_files):
            mask = filename_array == name
            for time_bin, z_bin in zip(binned[mask], z_binned[mask]):
                if 0 <= time_bin < len(edges) - 1 and 0 <= z_bin < num_z_bins:
                    per_file[file_idx, time_bin, z_bin] += 1
                    counts_per_file[file_idx, time_bin] += 1

    mean_counts = per_file.mean(axis=0) if unique_files else np.zeros((len(edges) - 1, num_z_bins))
    sem = counts_per_file.std(axis=0, ddof=1) / np.sqrt(len(unique_files)) if len(unique_files) > 1 else np.zeros(len(edges) - 1)
    mean_total = counts_per_file.mean(axis=0) if unique_files else np.zeros(len(edges) - 1)
    return centers, mean_counts, mean_total, sem, unique_files, z_edges


def kde_curve(peak_times: np.ndarray, rec_dur_min: int, bin_min: float) -> tuple[np.ndarray, np.ndarray]:
    x = np.linspace(0, rec_dur_min, 200)
    if peak_times.size < 2:
        return x, np.zeros_like(x)
    data_min = peak_times / 60.0
    sample_std = np.std(data_min, ddof=1)
    if not np.isfinite(sample_std) or sample_std == 0:
        return x, np.zeros_like(x)
    try:
        kde = gaussian_kde(data_min, bw_method=bin_min / sample_std)
        return x, kde(x) * peak_times.size
    except Exception:
        return x, np.zeros_like(x)


def save_temporal_psd_map(stim, sham, out_dir: Path, rec_dur_min: int, bin_min: float, z_increment: float, show_sem: bool, show_kde: bool) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    stim_times, stim_z, stim_names = stim
    sham_times, sham_z, sham_names = sham
    centers, stim_stack, stim_mean, stim_sem, stim_files, z_edges = stacked_counts(stim_times, stim_z, stim_names, rec_dur_min, bin_min, z_increment)
    _, sham_stack, sham_mean, sham_sem, sham_files, _ = stacked_counts(sham_times, sham_z, sham_names, rec_dur_min, bin_min, z_increment)
    num_z_bins = stim_stack.shape[1] if stim_stack.ndim == 2 else sham_stack.shape[1]
    cmap = plt.get_cmap('jet', num_z_bins)
    stim_x, stim_kde = kde_curve(stim_times, rec_dur_min, bin_min)
    sham_x, sham_kde = kde_curve(sham_times, rec_dur_min, bin_min)
    stim_totals = np.sum(stim_stack, axis=1) if stim_stack.size else np.array([0.0])
    sham_totals = np.sum(sham_stack, axis=1) if sham_stack.size else np.array([0.0])
    max_count = max(float(np.max(stim_totals)), float(np.max(sham_totals)), float(np.max(stim_mean + stim_sem)) if stim_mean.size else 0.0, float(np.max(sham_mean + sham_sem)) if sham_mean.size else 0.0, 1.0)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True, constrained_layout=True)
    for ax, stack, mean_vals, sem_vals, kde_x, kde_y, title in [(axes[0], sham_stack, sham_mean, sham_sem, sham_x, sham_kde, 'SHAM'), (axes[1], stim_stack, stim_mean, stim_sem, stim_x, stim_kde, 'STIM')]:
        bottom = np.zeros(len(centers))
        for idx in range(stack.shape[1]):
            ax.bar(centers, stack[:, idx], bottom=bottom, width=bin_min * 0.9, color=cmap(idx), edgecolor='black', linewidth=0.4)
            bottom += stack[:, idx]
        if show_sem and mean_vals.size:
            ax.errorbar(centers, mean_vals, yerr=sem_vals, fmt='k.', linewidth=1.2, capsize=3)
        ax.set_title(title)
        ax.set_xlabel('Time [min]')
        ax.set_xlim(0, rec_dur_min)
        ax.set_ylim(0, max_count * 1.1)
        ax.set_ylabel('Average number of peaks / file')
        ax.grid(True)
        if show_kde:
            twin = ax.twinx()
            twin.plot(kde_x, kde_y, color='red', linewidth=2)
            twin.set_ylabel('Peak density (KDE scaled)', color='red')
            twin.tick_params(axis='y', colors='red')
            twin.set_ylim(0, max(np.max(stim_kde), np.max(sham_kde), 1) * 1.1)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-4, vmax=4))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.03, pad=0.03)
    cbar.set_label('z-score range (σ)')
    cbar.set_ticks(np.arange(-4, 5, 2))
    fig.suptitle(f'Temporal distribution of detected peaks (bin = {bin_min:g} min)')
    fig.savefig(out_dir / 'Temporal_PSDmap.png', dpi=200)
    fig.savefig(out_dir / 'Temporal_PSDmap.eps', format='eps')
    plt.close(fig)
    return {'stim_peak_names': stim_names, 'sham_peak_names': sham_names, 'stim_times': stim_times, 'sham_times': sham_times, 'stim_z': stim_z, 'sham_z': sham_z, 'stim_file_count': len(stim_files), 'sham_file_count': len(sham_files)}


def save_zscore_distribution_plot(stim_z: np.ndarray, sham_z: np.ndarray, num_stim_files: int, num_sham_files: int, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    z_edges = np.arange(-4, 4.1, 0.1)
    z_centers = z_edges[:-1] + np.diff(z_edges) / 2
    counts_stim = np.histogram(stim_z, bins=z_edges)[0] / max(num_stim_files, 1)
    counts_sham = np.histogram(sham_z, bins=z_edges)[0] / max(num_sham_files, 1)
    ax.barh(z_centers, -counts_sham, color=(0.8, 0.3, 0.3), edgecolor='none', label='SHAM')
    ax.barh(z_centers, counts_stim, color=(0.2, 0.6, 0.8), edgecolor='none', label='STIM')
    ax.axvline(0, color='black', linewidth=1.5)
    ax.set_xlabel('Number of data points / file')
    ax.set_ylabel('Z-score (σ)')
    ax.set_ylim(-4, 4)
    ax.set_title('Z-score Distribution Bar Plot')
    ax.grid(True, axis='x')
    ax.legend(loc='best', framealpha=1.0, fancybox=False)
    fig.savefig(out_dir / 'Zscore_BarPlot_Distribution.png', dpi=200)
    fig.savefig(out_dir / 'Zscore_BarPlot_Distribution.eps', format='eps')
    plt.close(fig)


def save_raw_summary(out_dir: Path, stim_names: list[str], sham_names: list[str], stim_times: np.ndarray, sham_times: np.ndarray, stim_z: np.ndarray, sham_z: np.ndarray) -> None:
    rows = []
    for group, peak_names, times, zscores in [('STIM', stim_names, stim_times, stim_z), ('SHAM', sham_names, sham_times, sham_z)]:
        for file_name, peak_time, zscore in zip(peak_names, times, zscores):
            rows.append({'Group': group, 'Filename': file_name, 'PeakTime_sec': float(peak_time), 'Zscore': float(zscore)})
    csv_path = out_dir / 'temporal_peak_raw_data.csv'
    with csv_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Group', 'Filename', 'PeakTime_sec', 'Zscore'])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    if Workbook is not None:
        wb = Workbook()
        ws = wb.active
        ws.append(['Group', 'Filename', 'PeakTime_sec', 'Zscore'])
        for row in rows:
            ws.append([row['Group'], row['Filename'], row['PeakTime_sec'], row['Zscore']])
        wb.save(out_dir / 'temporal_peak_raw_data.xlsx')


def find_mat_files(inputs: Iterable[str]) -> list[Path]:
    files: list[Path] = []
    for item in inputs:
        path = Path(item)
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Python port of tFUS_EventAnalyzer.m + plotTemporalPSDmap.m')
    parser.add_argument('--stim', nargs='+', required=True, help='STIM .mat files or directories')
    parser.add_argument('--sham', nargs='+', required=True, help='SHAM .mat files or directories')
    parser.add_argument('--output-dir', default='analysis/output_py', help='Output directory')
    parser.add_argument('--fs', type=int, default=DEFAULT_FS, help='Sampling frequency in Hz')
    parser.add_argument('--recording-minutes', type=int, default=30, help='Recording duration in minutes')
    parser.add_argument('--bin-min', type=float, default=5.0, help='Time bin size in minutes')
    parser.add_argument('--z-increment', type=float, default=0.5, help='Z-score increment for stack bands')
    parser.add_argument('--no-sem', action='store_true', help='Disable SEM overlay')
    parser.add_argument('--no-kde', action='store_true', help='Disable KDE overlay')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = DetectionParams()
    stim_files = find_mat_files(args.stim)
    sham_files = find_mat_files(args.sham)
    if not stim_files or not sham_files:
        raise SystemExit('Both STIM and SHAM inputs must include at least one .mat file.')

    stim = collect_group(stim_files, args.fs, params)
    sham = collect_group(sham_files, args.fs, params)
    out_dir = Path(args.output_dir)
    summary = save_temporal_psd_map(stim, sham, out_dir, rec_dur_min=args.recording_minutes, bin_min=args.bin_min, z_increment=args.z_increment, show_sem=not args.no_sem, show_kde=not args.no_kde)
    save_zscore_distribution_plot(summary['stim_z'], summary['sham_z'], summary['stim_file_count'], summary['sham_file_count'], out_dir)
    save_raw_summary(out_dir, summary['stim_peak_names'], summary['sham_peak_names'], summary['stim_times'], summary['sham_times'], summary['stim_z'], summary['sham_z'])
    print(f'STIM files: {len(stim_files)}, peaks: {stim[0].size}')
    print(f'SHAM files: {len(sham_files)}, peaks: {sham[0].size}')
    print(f'Outputs written to {out_dir}')


if __name__ == '__main__':
    main()

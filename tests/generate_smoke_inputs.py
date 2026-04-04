from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.io import savemat


def make_dlc_csv(path: Path, frames: int = 90) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    t = np.arange(frames, dtype=float)
    x = 100.0 + 0.5 * t
    y = 50.0 + 0.25 * t
    likelihood = np.ones(frames, dtype=float) * 0.99
    likelihood[10] = 0.5
    likelihood[40] = 0.2

    with path.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['scorer', 'DLC_resnet50', 'DLC_resnet50', 'DLC_resnet50'])
        writer.writerow(['bodyparts', 'head', 'head', 'head'])
        writer.writerow(['coords', 'x', 'y', 'likelihood'])
        for idx in range(frames):
            writer.writerow([int(idx), f'{x[idx]:.3f}', f'{y[idx]:.3f}', f'{likelihood[idx]:.3f}'])


def make_eeg_mat(path: Path, samples: int = 2000, fs: int = 1000) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    t = np.arange(samples, dtype=float) / fs
    signal = 800 * np.sin(2 * np.pi * 8 * t) + 100 * np.sin(2 * np.pi * 40 * t)
    savemat(path, {'raw_data': signal.astype(np.float64)})


def main() -> None:
    parser = argparse.ArgumentParser(description='Generate tiny deterministic smoke-test inputs')
    parser.add_argument('--dlc-csv', type=Path, required=True)
    parser.add_argument('--eeg-mat', type=Path, required=True)
    args = parser.parse_args()

    make_dlc_csv(args.dlc_csv)
    make_eeg_mat(args.eeg_mat)


if __name__ == '__main__':
    main()

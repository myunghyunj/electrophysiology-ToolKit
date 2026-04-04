from __future__ import annotations

import csv
from pathlib import Path


def write_dlc_csv(path: Path, x_start: float, y_start: float, dx: float, dy: float, frames: int = 120) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["scorer", "synthetic_model", "synthetic_model", "synthetic_model"])
        writer.writerow(["bodyparts", "head", "head", "head"])
        writer.writerow(["coords", "x", "y", "likelihood"])
        for frame in range(frames):
            x = x_start + dx * frame
            y = y_start + dy * frame
            likelihood = 0.97 if frame % 11 else 0.85
            writer.writerow([frame, f"{x:.3f}", f"{y:.3f}", f"{likelihood:.3f}"])


def main() -> None:
    out_dir = Path("/tmp/locomotion_synthetic")
    write_dlc_csv(out_dir / "mouse_SHAM_head.csv", x_start=100.0, y_start=50.0, dx=0.75, dy=0.20)
    write_dlc_csv(out_dir / "mouse_STIM_head.csv", x_start=95.0, y_start=60.0, dx=0.60, dy=0.35)
    print(out_dir)


if __name__ == "__main__":
    main()

from __future__ import annotations

import glob
import shlex
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent


def ask(prompt: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{prompt}{suffix}: ").strip()
    return value or (default or "")


def parse_dragdrop_list(raw: str) -> list[str]:
    if not raw:
        return []
    return shlex.split(raw)


def run_script(script: Path, args: list[str]) -> None:
    cmd = [sys.executable, str(script), *args]
    print("\nRunning:", " ".join(shlex.quote(part) for part in cmd))
    subprocess.run(cmd, check=True)


def task_1() -> None:
    raw = ask("Drop one or more .mat files or directories for Task 1 (tic event analysis)")
    inputs = parse_dragdrop_list(raw)
    out_dir = ask("Output directory", "analyzed_figures_py")
    run_script(REPO_ROOT / "1. Tic Event Based Spectogram:EEG Analysis" / "eeg_spectrogram.py", [*inputs, "--output-dir", out_dir])


def task_2() -> None:
    stim = parse_dragdrop_list(ask("Drop STIM .mat file(s) or directory"))
    sham = parse_dragdrop_list(ask("Drop SHAM .mat file(s) or directory"))
    out_dir = ask("Output directory", "analysis/output_py")
    bin_min = ask("Bin size in minutes", "5")
    z_inc = ask("Z-score increment", "0.5")
    run_script(
        REPO_ROOT / "2. Peak Temporal Distrubution" / "temporal_peak_distribution.py",
        ["--stim", *stim, "--sham", *sham, "--output-dir", out_dir, "--bin-min", bin_min, "--z-increment", z_inc],
    )


def task_3() -> None:
    mode = ask("Task 3 mode: [1] adaptive scale, [2] fixed ±5000 µV", "1")
    raw = ask("Drop one .mat file, or enter a glob pattern such as './*.mat'")
    output = ask("Output video file (leave blank for automatic naming)", "")
    window = ask("Visible window size in seconds", "6")
    fps = ask("Video FPS", "30")
    script = "final.py" if mode != "2" else "final_fixed5000uv.py"

    paths = parse_dragdrop_list(raw)
    raw_target = raw.strip()
    parsed_glob = len(paths) == 1 and glob.has_magic(paths[0])
    raw_glob = glob.has_magic(raw_target)
    raw_file = Path(raw_target)

    if parsed_glob or raw_glob:
        target = paths[0] if parsed_glob else raw_target
        if output:
            raise SystemExit("Custom output filename is only supported for single-file mode.")
        args = ["--all", "--pattern", target, "--window", window, "--fps", fps]
    else:
        if raw_file.suffix.lower() == ".mat" and raw_file.exists():
            target = raw_target
        else:
            if len(paths) != 1:
                raise SystemExit('Task 3 expects exactly one .mat file or one glob pattern')
            target = paths[0]
        if Path(target).suffix.lower() != ".mat":
            raise SystemExit("Task 3 expects a .mat file or a glob pattern matching .mat files.")
        args = [target, "--window", window, "--fps", fps]
        if output:
            args += ["--output", output]

    run_script(REPO_ROOT / "3. EEG Monitor (1600 dpi, 30 fps for Fs=1000Hz)" / script, args)


def task_4() -> None:
    mode = ask("Task 4 mode: [1] DLC preprocess, [2] binned locomotion analysis", "1")
    script = REPO_ROOT / "4. Locomotion Analysis Tools" / "locomotion_analysis.py"
    if mode == "1":
        raw = ask("Drop DLC CSV file(s) or directory containing raw CSV files")
        inputs = parse_dragdrop_list(raw)
        out_dir = ask("Output directory", "4. Locomotion Analysis Tools/python_output")
        run_script(script, ["preprocess", *inputs, "--output-dir", out_dir])
    else:
        raw = ask("Drop tracking CSV/XLSX file(s) or directory")
        inputs = parse_dragdrop_list(raw)
        out_dir = ask("Output directory", "4. Locomotion Analysis Tools/python_output")
        bin_size = ask("Bin size in minutes", "5")
        mm_per_pix = ask("mm per pixel", "0.533")
        run_script(script, ["analyze", *inputs, "--output-dir", out_dir, "--bin-size-min", bin_size, "--mm-per-pix", mm_per_pix])


def main() -> None:
    print("Electrophysiology Toolkit Python CLI")
    print("  1) Tic Event / Spectrogram")
    print("  2) Temporal Peak Distribution")
    print("  3) EEG Video Monitor")
    print("  4) Locomotion Analysis")
    choice = ask("Choose task number", "1")
    if choice == "1":
        task_1()
    elif choice == "2":
        task_2()
    elif choice == "3":
        task_3()
    elif choice == "4":
        task_4()
    else:
        raise SystemExit(f"Unsupported task number: {choice}")


if __name__ == "__main__":
    main()

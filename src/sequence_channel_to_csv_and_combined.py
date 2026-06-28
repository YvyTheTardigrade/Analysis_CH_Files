"""Sequence/channel .ch extractor v20 (English edition).

This script scans Agilent .D folders, reads .ch files with the v12 parser,
exports one signal CSV and one header CSV per target, builds combined and
summary CSV files, and can optionally export plots.

Version 19 selects target folders whose names begin with a 3-digit number and contain -Pip_, and labels each target from the immediately preceding numbered folder name without restricting the label to amino acids.
"""

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import sys
import re
import struct
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


def get_matplotlib_pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ModuleNotFoundError as e:
        raise RuntimeError(
            "matplotlib is required for plot export. Install it with: python -m pip install matplotlib"
        ) from e


def open_file_with_default_app(path: Path) -> None:
    """Open a file with the operating-system default application."""
    # used for plot visualisation
    try:
        if sys.platform == "darwin":
            subprocess.Popen(["open", str(path)])
        elif os.name == "nt":
            os.startfile(str(path))  # type: ignore[attr-defined]
        else:
            subprocess.Popen(["xdg-open", str(path)])
    except Exception as exc:
        print(f"[WARNING] Could not open summary plot automatically: {exc}")

"""
AMINO_ACIDS: Set[str] = {
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly",
    "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser",
    "Thr", "Trp", "Tyr", "Val",
}
"""

PeakRecord = Dict[str, object]


# ---------------------------------------------------------------------
# Big-endian readers
# ---------------------------------------------------------------------
def be_i16(b: bytes, off: int) -> int:
    return int.from_bytes(b[off:off + 2], byteorder="big", signed=True)


def be_u16(b: bytes, off: int) -> int:
    return int.from_bytes(b[off:off + 2], byteorder="big", signed=False)


def be_u32(b: bytes, off: int) -> int:
    return int.from_bytes(b[off:off + 4], byteorder="big", signed=False)


def be_f64(b: bytes, off: int) -> float:
    return struct.unpack(">d", b[off:off + 8])[0] # ">" = big-endian  , "d" = double precision 


def be_i32_from_words(hi_word: int, lo_word: int) -> int:
    value = (hi_word << 16) | lo_word
    if value >= 0x80000000:
        value -= 0x100000000
    return value


# ---------------------------------------------------------------------
# Header reader derived from ch_to_csv_two_csvs.py
# ---------------------------------------------------------------------
def read_pascal_utf16_like(data: bytes, offset: int) -> str:
    if offset >= len(data):
        return ""

    n = data[offset]
    if n <= 0:
        return ""

    start = offset + 1
    end = start + 2 * n
    if end > len(data):
        return ""

    raw = data[start:end:2]
    try:
        return raw.decode("latin1").strip()
    except Exception:
        return "".join(chr(x) for x in raw).strip()



def extract_header(data: bytes) -> Dict[str, Any]:
    header: Dict[str, Any] = {}

    header["file_type_number"] = read_pascal_utf16_like(data, 0x0146)
    header["file_type_name"] = read_pascal_utf16_like(data, 0x015B)
    header["notebook_name"] = read_pascal_utf16_like(data, 0x035A)
    header["parent_directory"] = read_pascal_utf16_like(data, 0x0758)
    header["date"] = read_pascal_utf16_like(data, 0x0957)
    header["unknown_09BC"] = read_pascal_utf16_like(data, 0x09BC)
    header["unknown_09E5"] = read_pascal_utf16_like(data, 0x09E5)
    header["method"] = read_pascal_utf16_like(data, 0x0A0E)
    header["instrument"] = read_pascal_utf16_like(data, 0x0C11)
    header["version_string"] = read_pascal_utf16_like(data, 0x0E11)
    header["unknown_0EDA"] = read_pascal_utf16_like(data, 0x0EDA)
    header["units"] = read_pascal_utf16_like(data, 0x104C)
    header["signal_name"] = read_pascal_utf16_like(data, 0x1075)

    header["first_time_ms"] = be_u32(data, 0x011A)
    header["last_time_ms"] = be_u32(data, 0x011E)
    header["scaling_factor"] = be_f64(data, 0x127C)
    header["header_size_bytes"] = 0x1800

    return header


# ---------------------------------------------------------------------
# Signal decoding
# ---------------------------------------------------------------------
def decode_signal_from_offset(data: bytes, start: int) -> List[int]:
    values: List[int] = []
    current = 0
    i = start
    n = len(data)

    while i + 1 < n:
        word = be_i16(data, i)

        if word == -32768:
            if i + 5 >= n:
                break
            hi = be_u16(data, i + 2)
            lo = be_u16(data, i + 4)
            current = be_i32_from_words(hi, lo)
            values.append(current)
            i += 6
        else:
            current += word
            values.append(current)
            i += 2

    return values



def align_signal(signal: List[int], drop: int, target_points: Optional[int]) -> List[int]:
    if drop < 0:
        raise ValueError("drop must be >= 0")

    trimmed = signal[drop:] if len(signal) > drop else []
    if not trimmed:
        raise RuntimeError("Decoded signal is empty after alignment.")

    if target_points is None:
        return trimmed

    if len(trimmed) >= target_points:
        return trimmed[:target_points]

    out = trimmed[:]
    while len(out) < target_points:
        if len(out) >= 2:
            extra = out[-1] + (out[-1] - out[-2])
        else:
            extra = out[-1]
        out.append(extra)
    return out


def deduce_output_point_count(
    first_time_ms: Optional[int],
    last_time_ms: Optional[int],
    decoded_point_count: int,
) -> Optional[int]:
    if (
        first_time_ms is None
        or last_time_ms is None
        or decoded_point_count <= 0
    ):
        return None

    if decoded_point_count == 1:
        return 1

    delta_ms = (last_time_ms - first_time_ms) / (decoded_point_count - 1)
    if delta_ms == 0:
        return decoded_point_count

    # User-requested logic based on:
    # (last_time_ms - first_time_ms) / (time_ms(2) - time_ms(1))
    deduced = int(round((last_time_ms - first_time_ms) / delta_ms)) + 1
    return max(deduced, 1)


def add_header_timing_details(header: Dict[str, Any], decoded_point_count: int) -> None:
    first_time_ms = header.get("first_time_ms")
    last_time_ms = header.get("last_time_ms")
    header["decoded_point_count"] = decoded_point_count
    deduced_count = deduce_output_point_count(first_time_ms, last_time_ms, decoded_point_count)
    header["deduced_output_point_count"] = deduced_count if deduced_count is not None else ""
    if (
        first_time_ms is not None
        and last_time_ms is not None
        and decoded_point_count > 1
    ):
        delta_ms = (last_time_ms - first_time_ms) / (decoded_point_count - 1)
        header["second_time_ms"] = first_time_ms + delta_ms
        header["time_step_ms"] = delta_ms
    else:
        header["second_time_ms"] = ""
        header["time_step_ms"] = ""



def build_time_axis(
    n_points: int,
    total_time_s: Optional[float],
    time_start_s: float,
    first_time_ms: Optional[int] = None,
    last_time_ms: Optional[int] = None,
    use_header_time: bool = True,
) -> List[float]:
    if use_header_time and first_time_ms is not None and last_time_ms is not None and n_points > 0:
        if n_points == 1:
            return [first_time_ms / 1000.0]
        dt = (last_time_ms - first_time_ms) / (n_points - 1)
        return [(first_time_ms + i * dt) / 1000.0 for i in range(n_points)]

    if total_time_s is None:
        total_time_s = 2.0

    if n_points < 2:
        return [time_start_s]

    dt = total_time_s / (n_points - 1)
    return [time_start_s + i * dt for i in range(n_points)]



def calculate_integral(times: List[float], signal: List[int]) -> float:
    if len(times) != len(signal):
        raise ValueError("times and signal must have the same length")

    n = len(signal)
    if n < 2:
        return 0.0

    dt = times[1] - times[0]
    if n == 2:
        return dt * (signal[0] + signal[1]) / 2.0

    intervals = n - 1
    if intervals % 2 == 0:
        odd_sum = sum(signal[i] for i in range(1, n - 1, 2))
        even_sum = sum(signal[i] for i in range(2, n - 1, 2))
        return (dt / 3.0) * (signal[0] + signal[-1] + 4.0 * odd_sum + 2.0 * even_sum)

    usable_n = n - 1
    odd_sum = sum(signal[i] for i in range(1, usable_n - 1, 2))
    even_sum = sum(signal[i] for i in range(2, usable_n - 1, 2))
    simpson_area = (dt / 3.0) * (signal[0] + signal[usable_n - 1] + 4.0 * odd_sum + 2.0 * even_sum)
    trapezoid_tail = dt * (signal[usable_n - 1] + signal[usable_n]) / 2.0
    return simpson_area + trapezoid_tail


# ---------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------
def write_header_csv(out_csv: Path, header: Dict[str, Any]) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["header_key", "header_value"])
        for key, value in header.items():
            writer.writerow([key, value])



def write_signal_csv(
    out_csv: Path,
    times: List[float],
    signal_raw: List[int],
    signal_scaled: List[float],
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["point", "time_s", "signal_raw", "signal_scaled"])
        for i, (t, raw_y, scaled_y) in enumerate(zip(times, signal_raw, signal_scaled)):
            writer.writerow([i, f"{t:.9f}", raw_y, f"{scaled_y:.12g}"])



def read_ch_data(
    ch_file: Path,
    offset: Optional[int],
    target_points: Optional[int],
    run_time_s: Optional[float],
    time_start_s: float,
    use_header_time: bool,
) -> Tuple[Dict[str, Any], List[float], List[int], List[float]]:
    if not ch_file.exists():
        raise FileNotFoundError(f"Missing input file: {ch_file}")

    data = ch_file.read_bytes()
    header = extract_header(data)

    signal_offset = offset if offset is not None else int(header.get("header_size_bytes", 0x1800))

    raw_signal_all = decode_signal_from_offset(data, signal_offset)
    if not raw_signal_all:
        raise RuntimeError(f"Could not decode signal from {ch_file.name}")

    add_header_timing_details(header, decoded_point_count=len(raw_signal_all))

    effective_target_points = target_points
    if effective_target_points is None:
        effective_target_points = deduce_output_point_count(
            header.get("first_time_ms"),
            header.get("last_time_ms"),
            len(raw_signal_all),
        )

    signal_raw = align_signal(raw_signal_all, drop=0, target_points=effective_target_points)
    scaling = float(header.get("scaling_factor", 1.0))
    signal_scaled = [x * scaling for x in signal_raw]

    times = build_time_axis(
        n_points=len(signal_raw),
        total_time_s=run_time_s,
        time_start_s=time_start_s,
        first_time_ms=header.get("first_time_ms"),
        last_time_ms=header.get("last_time_ms"),
        use_header_time=use_header_time,
    )

    return header, times, signal_raw, signal_scaled


# ---------------------------------------------------------------------
# Folder helpers
# ---------------------------------------------------------------------
def extract_folder_order(folder_name: str) -> Optional[int]:
    """Return the numeric order only for folders beginning with exactly 3 digits and a dash.

    Examples accepted: 001-P1-F1-Pip_1.D, 017-P1-B1-Ile.D
    Examples rejected: 1-P1-F1-Pip_1.D, ABC-P1-F1-Pip_1.D
    """
    m = re.match(r"^(\d{3})-", folder_name)
    return int(m.group(1)) if m else None


def extract_previous_data_name(folder_name: str) -> Optional[str]:
    """Extract the data name from a numbered .D folder without restricting the name.

    The name is taken from the text after the final dash and before .D.
    This intentionally allows names such as Ile, Lys2, Dawson link, t, and LMt.
    """
    m = re.match(r"^\d{3}-.+-(.+)\.D$", folder_name)
    if not m:
        return None
    name = m.group(1).strip()
    return name or None


def is_pip_folder(folder_name: str) -> bool:
    """Target folders must start with 3 digits and contain -Pip_."""
    return re.match(r"^\d{3}-", folder_name) is not None and "-Pip_" in folder_name and folder_name.endswith(".D")


def collect_ordered_subfolders(parent_dir: Path) -> List[Path]:
    subdirs = [p for p in parent_dir.iterdir() if p.is_dir()]
    parsed: List[Tuple[int, Path]] = []
    for p in subdirs:
        order = extract_folder_order(p.name)
        if order is not None:
            parsed.append((order, p))
    parsed.sort(key=lambda x: x[0])
    return [p for _, p in parsed]


def build_pip_targets(parent_dir: Path) -> List[Tuple[Path, str, str]]:
    """Build targets from Pip folders and label them with the preceding folder's data name.

    For example, if 018-P1-F1-Pip_1.D follows 017-P1-B1-Ile.D, the
    target label name is Ile. The target counter still increases by one for
    each processed Pip folder: 1_Ile, 2_..., etc.
    """
    ordered_folders = collect_ordered_subfolders(parent_dir)
    previous_data_name: Optional[str] = None
    targets: List[Tuple[Path, str, str]] = []

    for folder in ordered_folders:
        if is_pip_folder(folder.name):
            if previous_data_name is None:
                print(f"[SKIP] {folder.name}: no preceding numbered data folder")
                continue
            folder_number = folder.name.split("-", 1)[0]
            targets.append((folder, folder_number, previous_data_name))
            continue

        data_name = extract_previous_data_name(folder.name)
        if data_name is not None:
            previous_data_name = data_name

    return targets


def sanitize_channel_stem(channel_filename: str) -> str:
    return Path(channel_filename).stem



def write_combined_csv(combined_csv: Path, combined_rows: List[Tuple[str, int, float, int, float, str, str]]) -> None:
    combined_csv.parent.mkdir(parents=True, exist_ok=True)
    with combined_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "sequence_label",
            "global_point_index",
            "time_s",
            "signal_raw",
            "signal_scaled",
            "source_folder",
            "integral",
        ])
        for label, point_index, t, raw_y, scaled_y, source_folder, integral in combined_rows:
            writer.writerow([label, point_index, f"{t:.9f}", raw_y, f"{scaled_y:.12g}", source_folder, integral])



def write_summary_csv(summary_csv: Path, summary_rows: List[Tuple[str, str, str, str, str, str, str, str, str]]) -> None:
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "sequence_label",
            "source_folder",
            "amino_acid",
            "folder_number",
            "channel_file_name",
            "signal_csv_path",
            "header_csv_path",
            "integral_raw",
            "integral_normalized",
        ])
        for row in summary_rows:
            writer.writerow(row)



def _label_position_x(point_indices: List[int], max_idx: int) -> int:
    if len(point_indices) < 2:
        return point_indices[max_idx]
    span = point_indices[-1] - point_indices[0]
    shift = max(1, int(round(span * 0.03)))
    return point_indices[max_idx] - shift



def _draw_peak_panel(ax1, peak: PeakRecord, show_ylabel_left: bool = True, show_ylabel_right: bool = True) -> None:
    point_indices = peak["point_indices"]
    signal = peak["signal"]
    normalized_line = peak["normalized_line"]
    peak_label = peak["label"]
    normalized_integral = peak["normalized_integral"]
    total_integral = peak["total_integral"]

    ax1.plot(point_indices, signal, color="black", linewidth=1.0)
    ax1.set_xlabel("global point index")
    if show_ylabel_left:
        ax1.set_ylabel("signal", color="black")
    ax1.tick_params(axis="y", labelcolor="black")

    ax2 = ax1.twinx()
    ax2.plot(point_indices, normalized_line, color="red", linewidth=1.2)

    max_idx = max(range(len(signal)), key=lambda i: signal[i])
    label_x = _label_position_x(point_indices, max_idx)
    ax2.plot([point_indices[max_idx]], [normalized_integral], marker="o", color="red", markersize=4)
    ax1.annotate(
        peak_label,
        xy=(label_x, signal[max_idx]),
        xytext=(0, 6),
        textcoords="offset points",
        fontsize=8,
        color="black",
        rotation=90,
        va="bottom",
        ha="center",
    )

    if show_ylabel_right:
        ax2.set_ylabel("normalized integral", color="red")
    ax2.tick_params(axis="y", labelcolor="red")
    ax2.set_ylim(0.0, 1.05 if normalized_integral <= 1.0 else max(1.05, normalized_integral * 1.05))
    ax1.set_title(f"{peak_label} | integral = {total_integral:.9f} | normalized = {normalized_integral:.6f}", fontsize=9)



def export_individual_subplot_pages(plots_dir: Path, channel_stem: str, peak_records: List[PeakRecord], rows: int = 3, cols: int = 2) -> List[Path]:
    plt = get_matplotlib_pyplot()
    plots_dir.mkdir(parents=True, exist_ok=True)
    page_paths: List[Path] = []
    per_page = rows * cols

    for page_idx in range(math.ceil(len(peak_records) / per_page)):
        start = page_idx * per_page
        end = min(start + per_page, len(peak_records))
        page_peaks = peak_records[start:end]

        fig, axes = plt.subplots(rows, cols, figsize=(14, 12), dpi=150)
        axes_flat = list(axes.flat)

        for ax, peak in zip(axes_flat, page_peaks):
            _draw_peak_panel(ax, peak)

        for ax in axes_flat[len(page_peaks):]:
            ax.axis("off")

        fig.suptitle(f"Individual peaks ({channel_stem}) - page {page_idx + 1}", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        out_png = plots_dir / f"individual_plots_page_{page_idx + 1:02d}.png"
        fig.savefig(out_png, format="png", bbox_inches="tight")
        plt.close(fig)
        page_paths.append(out_png)

    return page_paths



def export_combined_plot(out_png: Path, combined_rows: List[Tuple[str, int, float, int, float, str, str]], peak_records: List[PeakRecord]) -> None:
    plt = get_matplotlib_pyplot()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    point_indices = [row[1] for row in combined_rows]
    signal = [row[3] for row in combined_rows]

    normalized_series: List[float] = []
    for peak in peak_records:
        normalized_series.extend(peak["normalized_line"])

    fig, ax1 = plt.subplots(figsize=(15, 6), dpi=150)
    ax1.plot(point_indices, signal, color="black", linewidth=0.8)
    ax1.set_xlabel("global point index")
    ax1.set_ylabel("signal", color="black")
    ax1.tick_params(axis="y", labelcolor="black")

    ax2 = ax1.twinx()
    ax2.plot(point_indices, normalized_series, color="red", linewidth=1.0)

    for peak in peak_records:
        peak_points = peak["point_indices"]
        peak_signal = peak["signal"]
        peak_label = peak["label"]
        normalized_integral = peak["normalized_integral"]
        max_idx = max(range(len(peak_signal)), key=lambda i: peak_signal[i])
        label_x = _label_position_x(peak_points, max_idx)
        ax2.plot([peak_points[max_idx]], [normalized_integral], marker="o", color="red", markersize=4)
        ax1.annotate(
            peak_label,
            xy=(label_x, peak_signal[max_idx]),
            xytext=(0, 6),
            textcoords="offset points",
            fontsize=8,
            color="black",
            rotation=90,
            va="bottom",
            ha="center",
        )

    ax2.set_ylabel("normalized integral", color="red")
    ax2.tick_params(axis="y", labelcolor="red")
    ax2.set_ylim(0.0, 1.05)
    ax1.set_title("Combined raw signal and normalized integrals")
    fig.tight_layout()
    fig.savefig(out_png, format="png", bbox_inches="tight")
    plt.close(fig)



def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Scan a parent folder of Agilent .D folders, assign each Pip_x.D folder "
            "to the immediately preceding numbered data folder, convert a chosen channel "
            "file using the v12 .ch reader, write one signal CSV and one header CSV per target, "
            "one combined CSV, one summary CSV, and optionally export PNG plots. "
            "Labels use a sequential integer starting at 1 followed by the preceding data name. "
            "The red plot line uses normalized peak integrals. Version 12.1 uses "
            "fully harmonized English output names."
        )
    )
    parser.add_argument("parent_dir", type=Path, help='Parent folder, e.g. "58mer 2026-02-18 09-40-30"')
    parser.add_argument("channel_file", type=str, help='Channel file inside each Pip folder, e.g. "DAD1D.ch"')
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory where output CSV files will be stored")
    parser.add_argument("--combined-name", type=str, default="combined.csv", help='Filename for the combined CSV in output-dir')
    parser.add_argument("--summary-name", type=str, default="summary.csv", help='Filename for the summary CSV in output-dir')
    parser.add_argument("--offset", type=int, default=None, help="Signal start offset in bytes. Default: header_size_bytes from the .ch file header")
    parser.add_argument("--points", type=int, default=None, help="Optional manual output point count. Default: deduced automatically from the .ch header timing")
    parser.add_argument("--run-time", type=float, default=None, help="Fallback total run time in seconds when header time is not used")
    parser.add_argument("--time-start", type=float, default=0.0, help="Fallback start time in seconds when header time is not used")
    parser.add_argument("--no-header-time", action="store_true", help="Use --run-time/--time-start instead of first_time_ms/last_time_ms from the .ch header (header times are converted from ms to s)")
    parser.add_argument("--export-csv", action="store_true", help="Export CSV files into output-dir/csv")
    parser.add_argument("--export-summary-plot", action="store_true", help="Export the summary plot to output-dir/plots/summary/combined_plot.png")
    parser.add_argument("--export-detail-plots", action="store_true", help="Export all detail plots to output-dir/plots/detail")
    parser.add_argument("--show-summary-plot", action="store_true", help="Open the summary plot at the end of the run using the system default image viewer")
    args = parser.parse_args()

    parent_dir = args.parent_dir
    output_dir = args.output_dir
    channel_file = args.channel_file

    if not parent_dir.exists() or not parent_dir.is_dir():
        raise RuntimeError(f"Parent directory not found: {parent_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = output_dir / "csv"
    summary_plots_dir = output_dir / "plots" / "summary"
    detail_plots_dir = output_dir / "plots" / "detail"
    if args.export_csv:
        csv_dir.mkdir(parents=True, exist_ok=True)
    if args.export_summary_plot:
        summary_plots_dir.mkdir(parents=True, exist_ok=True)
    if args.export_detail_plots:
        detail_plots_dir.mkdir(parents=True, exist_ok=True)

    targets = build_pip_targets(parent_dir)
    if not targets:
        raise RuntimeError("No Pip_x.D folders with a preceding numbered data folder were found.")

    channel_stem = sanitize_channel_stem(channel_file)
    combined_csv = csv_dir / args.combined_name
    summary_csv = csv_dir / args.summary_name
    combined_plot_path = summary_plots_dir / "combined_plot.png"

    print(f"Parent directory: {parent_dir}")
    print(f"Channel file:     {channel_file}")
    print(f"Output directory: {output_dir}")
    if args.export_csv:
        print(f"CSV directory:    {csv_dir}")
        print(f"Combined CSV:     {combined_csv}")
        print(f"Summary CSV:      {summary_csv}")
    else:
        print("CSV export:       disabled")
    print(f"Targets found:    {len(targets)}")
    print(f"Header-based time axis: {'disabled' if args.no_header_time else 'enabled'}")
    if args.export_summary_plot:
        print(f"Summary plot:     {combined_plot_path}")
    if args.show_summary_plot:
        print("Show summary plot: enabled")
    if args.export_detail_plots:
        print(f"Detail plots dir: {detail_plots_dir}")
    print()

    ok = 0
    fail = 0
    global_point_index = 1
    combined_rows: List[Tuple[str, int, float, int, float, str, str]] = []
    peak_records: List[PeakRecord] = []

    for file_counter, (pip_folder, folder_number, data_name) in enumerate(targets, start=1):
        ch_file = pip_folder / channel_file
        label = f"{file_counter}_{data_name}"
        signal_csv = csv_dir / f"{label}_{channel_stem}_signal.csv"
        header_csv = csv_dir / f"{label}_{channel_stem}_header.csv"
        source_folder = pip_folder.name

        print(f"[PROCESSING] {pip_folder.name}")
        print(f"             Label: {label}")
        print(f"             Input .ch: {ch_file}")
        print(f"             Signal CSV: {signal_csv}")
        print(f"             Header CSV: {header_csv}")

        try:
            header, times, signal_raw, signal_scaled = read_ch_data(
                ch_file=ch_file,
                offset=args.offset,
                target_points=args.points,
                run_time_s=args.run_time,
                time_start_s=args.time_start,
                use_header_time=not args.no_header_time,
            )
            point_indices = list(range(global_point_index, global_point_index + len(signal_raw)))
            total_integral = calculate_integral(times=times, signal=signal_raw)

            if args.export_csv:
                write_signal_csv(out_csv=signal_csv, times=times, signal_raw=signal_raw, signal_scaled=signal_scaled)
                write_header_csv(out_csv=header_csv, header=header)

            peak_records.append({
                "label": label,
                "source_folder": source_folder,
                "amino_acid": data_name,
                "folder_number": folder_number,
                "channel_file": channel_file,
                "signal_csv": str(signal_csv),
                "header_csv": str(header_csv),
                "times": times,
                "signal": signal_raw,
                "signal_scaled": signal_scaled,
                "point_indices": point_indices,
                "total_integral": total_integral,
            })

            for point_index, (t, raw_y, scaled_y) in zip(point_indices, zip(times, signal_raw, signal_scaled)):
                is_last_point = point_index == point_indices[-1]
                integral_value = f"{total_integral:.9f}" if is_last_point else ""
                combined_rows.append((label, point_index, t, raw_y, scaled_y, source_folder, integral_value))

            global_point_index += len(signal_raw)
            print("             Status: OK\n")
            ok += 1

        except Exception as e:
            print(f"             Status: FAILED - {e}\n")
            fail += 1

    if peak_records:
        max_integral = max(float(peak["total_integral"]) for peak in peak_records)
        if max_integral == 0:
            max_integral = 1.0
        for peak in peak_records:
            normalized_integral = float(peak["total_integral"]) / max_integral
            peak["normalized_integral"] = normalized_integral
            peak["normalized_line"] = [normalized_integral] * len(peak["point_indices"])

    summary_rows: List[Tuple[str, str, str, str, str, str, str, str, str]] = []
    for peak in peak_records:
        summary_rows.append((
            str(peak["label"]),
            str(peak["source_folder"]),
            str(peak["amino_acid"]),
            str(peak["folder_number"]),
            str(peak["channel_file"]),
            str(peak["signal_csv"]),
            str(peak["header_csv"]),
            f"{float(peak['total_integral']):.9f}",
            f"{float(peak['normalized_integral']):.9f}",
        ))

    if args.export_csv:
        write_combined_csv(combined_csv=combined_csv, combined_rows=combined_rows)
        write_summary_csv(summary_csv=summary_csv, summary_rows=summary_rows)

    combined_plot = None
    page_paths: List[Path] = []
    if args.export_detail_plots and peak_records:
        page_paths = export_individual_subplot_pages(plots_dir=detail_plots_dir, channel_stem=channel_stem, peak_records=peak_records, rows=3, cols=2)
    if args.export_summary_plot and peak_records:
        combined_plot = combined_plot_path
        export_combined_plot(out_png=combined_plot, combined_rows=combined_rows, peak_records=peak_records)

    print("Run summary")
    print(f"  Converted successfully: {ok}")
    print(f"  Failed:                 {fail}")
    if args.export_csv:
        print(f"  CSV directory:          {csv_dir}")
        print(f"  Combined CSV:           {combined_csv}")
        print(f"  Summary CSV:            {summary_csv}")
    if args.export_detail_plots:
        print(f"  Detail plots directory: {detail_plots_dir}")
        for page_path in page_paths:
            print(f"  Individual plot page:   {page_path}")
    if args.export_summary_plot and combined_plot is not None:
        print(f"  Summary plot:           {combined_plot}")
        if args.show_summary_plot:
            print("  Opening summary plot...")
            open_file_with_default_app(combined_plot)
    elif args.show_summary_plot:
        print("  Summary plot was not opened because --export-summary-plot was not enabled or no peaks were processed.")

    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

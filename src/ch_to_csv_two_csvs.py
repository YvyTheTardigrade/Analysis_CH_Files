#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import struct
from pathlib import Path
from typing import Any, Dict, List


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
    return struct.unpack(">d", b[off:off + 8])[0]


def be_i32_from_words(hi_word: int, lo_word: int) -> int:
    value = (hi_word << 16) | lo_word
    if value >= 0x80000000:
        value -= 0x100000000
    return value


# ---------------------------------------------------------------------
# Header string reader
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


def align_signal(signal: List[int], drop: int, target_points: int | None) -> List[int]:
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


def build_time_axis(
    n_points: int,
    total_time_min: float | None,
    time_start_min: float = 0.0,
    first_time_ms: int | None = None,
    last_time_ms: int | None = None,
) -> List[float]:
    if first_time_ms is not None and last_time_ms is not None and n_points > 0:
        if n_points == 1:
            return [first_time_ms / 60000.0]
        dt = (last_time_ms - first_time_ms) / (n_points - 1)
        return [(first_time_ms + i * dt) / 60000.0 for i in range(n_points)]

    if total_time_min is None:
        total_time_min = 2.0

    if n_points < 2:
        return [time_start_min]

    dt = total_time_min / (n_points - 1)
    return [time_start_min + i * dt for i in range(n_points)]


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
        writer.writerow(["point", "time_min", "signal_raw", "signal_scaled"])
        for i, (t, raw_y, scaled_y) in enumerate(zip(times, signal_raw, signal_scaled)):
            writer.writerow([i, f"{t:.9f}", raw_y, f"{scaled_y:.12g}"])


# ---------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------
def convert_ch_to_csvs(
    ch_file: Path,
    signal_csv: Path,
    header_csv: Path,
    offset: int = 6172,
    drop: int = 10,
    points: int | None = 601,
    run_time: float | None = None,
    time_start: float = 0.0,
    use_header_time: bool = True,
) -> None:
    if not ch_file.exists():
        raise FileNotFoundError(f"Missing input file: {ch_file}")

    data = ch_file.read_bytes()

    header = extract_header(data)

    raw_signal_all = decode_signal_from_offset(data, offset)
    if not raw_signal_all:
        raise RuntimeError(f"Could not decode signal from {ch_file.name}")

    signal_raw = align_signal(raw_signal_all, drop=drop, target_points=points)

    scaling = float(header.get("scaling_factor", 1.0))
    signal_scaled = [x * scaling for x in signal_raw]

    if use_header_time:
        times = build_time_axis(
            n_points=len(signal_raw),
            total_time_min=run_time,
            time_start_min=time_start,
            first_time_ms=header.get("first_time_ms"),
            last_time_ms=header.get("last_time_ms"),
        )
    else:
        times = build_time_axis(
            n_points=len(signal_raw),
            total_time_min=run_time if run_time is not None else 2.0,
            time_start_min=time_start,
        )

    write_header_csv(header_csv, header)
    write_signal_csv(signal_csv, times, signal_raw, signal_scaled)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert Agilent .ch file to separate signal CSV and header CSV."
    )
    parser.add_argument("input_ch", type=Path, help="Input .ch file")
    parser.add_argument("--signal-csv", type=Path, help="Output signal CSV")
    parser.add_argument("--header-csv", type=Path, help="Output header CSV")
    parser.add_argument("--offset", type=int, default=6172, help="Signal start offset in bytes")
    parser.add_argument("--drop", type=int, default=10, help="Number of initial decoded points to drop")
    parser.add_argument("--points", type=int, default=601, help="Number of output points")
    parser.add_argument("--run-time", type=float, default=None, help="Run time in minutes if not using header time")
    parser.add_argument("--time-start", type=float, default=0.0, help="Start time in minutes")
    parser.add_argument(
        "--no-header-time",
        action="store_true",
        help="Use --run-time instead of first_time_ms/last_time_ms from header",
    )
    args = parser.parse_args()

    input_ch = args.input_ch
    stem = input_ch.with_suffix("")

    signal_csv = args.signal_csv or Path(f"{stem}_signal.csv")
    header_csv = args.header_csv or Path(f"{stem}_header.csv")

    convert_ch_to_csvs(
        ch_file=input_ch,
        signal_csv=signal_csv,
        header_csv=header_csv,
        offset=args.offset,
        drop=args.drop,
        points=args.points,
        run_time=args.run_time,
        time_start=args.time_start,
        use_header_time=not args.no_header_time,
    )

    print(f"Saved signal CSV: {signal_csv}")
    print(f"Saved header CSV: {header_csv}")


if __name__ == "__main__":
    main()# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


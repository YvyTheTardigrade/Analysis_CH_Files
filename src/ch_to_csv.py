#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List


def be_i16(b: bytes, off: int) -> int:
    return int.from_bytes(b[off:off + 2], byteorder="big", signed=True)


def be_u16(b: bytes, off: int) -> int:
    return int.from_bytes(b[off:off + 2], byteorder="big", signed=False)


def be_i32_from_words(hi_word: int, lo_word: int) -> int:
    value = (hi_word << 16) | lo_word
    if value >= 0x80000000:
        value -= 0x100000000
    return value


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


def build_time_axis(n_points: int, total_time_min: float, time_start_min: float) -> List[float]:
    if n_points < 2:
        return [time_start_min]
    dt = total_time_min / (n_points - 1)
    return [time_start_min + i * dt for i in range(n_points)]


def convert_ch_to_csv(
    ch_file: Path,
    out_csv: Path,
    offset: int = 6172,
    drop: int = 10,
    points: int | None = 601,
    run_time: float = 2.0,
    time_start: float = 0.0,
) -> None:
    if not ch_file.exists():
        raise FileNotFoundError(f"Missing input file: {ch_file}")

    data = ch_file.read_bytes()
    raw_signal = decode_signal_from_offset(data, offset)
    if not raw_signal:
        raise RuntimeError(f"Could not decode signal from {ch_file.name}")

    signal = align_signal(raw_signal, drop=drop, target_points=points)
    times = build_time_axis(len(signal), total_time_min=run_time, time_start_min=time_start)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["point", "time_min", "signal"])
        for i, (t, y) in enumerate(zip(times, signal)):
            writer.writerow([i, f"{t:.9f}", y])


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert one .ch file to CSV.")
    parser.add_argument("input_ch", type=Path, help="Input .ch file")
    parser.add_argument("output_csv", type=Path, nargs="?", help="Output CSV file (default: same name as input)")
    parser.add_argument("--offset", type=int, default=6172, help="Signal start offset in bytes")
    parser.add_argument("--drop", type=int, default=10, help="Number of initial decoded points to drop")
    parser.add_argument("--points", type=int, default=601, help="Number of output points")
    parser.add_argument("--run-time", type=float, default=2.0, help="Total run time in minutes")
    parser.add_argument("--time-start", type=float, default=0.0, help="Start time in minutes")
    args = parser.parse_args()

    output_csv = args.output_csv or args.input_ch.with_suffix(".csv")
    convert_ch_to_csv(
        ch_file=args.input_ch,
        out_csv=output_csv,
        offset=args.offset,
        drop=args.drop,
        points=args.points,
        run_time=args.run_time,
        time_start=args.time_start,
    )
    print(f"Saved CSV: {output_csv}")


if __name__ == "__main__":
    main()
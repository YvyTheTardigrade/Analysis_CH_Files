from binary_helpers import *
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple



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



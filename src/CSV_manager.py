from typing import Any, Dict, List, Optional, Set, Tuple
import csv 
from pathlib import Path

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




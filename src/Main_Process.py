import CH_file_parser as ch
import plot_helpers as plot
from plot_helpers import PeakRecord
import CSV_manager as csv
import integrator 
import folder_helpers as folder

import argparse
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


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

    targets = folder.build_pip_targets(parent_dir)
    if not targets:
        raise RuntimeError("No Pip_x.D folders with a preceding numbered data folder were found.")

    channel_stem = folder.sanitize_channel_stem(channel_file)
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
            header, times, signal_raw, signal_scaled = ch.read_ch_data(
                ch_file=ch_file,
                offset=args.offset,
                target_points=args.points,
                run_time_s=args.run_time,
                time_start_s=args.time_start,
                use_header_time=not args.no_header_time,
            )
            point_indices = list(range(global_point_index, global_point_index + len(signal_raw)))
            total_integral = integrator.calculate_integral(times=times, signal=signal_raw)

            if args.export_csv:
                csv.write_signal_csv(out_csv=signal_csv, times=times, signal_raw=signal_raw, signal_scaled=signal_scaled)
                csv.write_header_csv(out_csv=header_csv, header=header)

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
        csv.write_combined_csv(combined_csv=combined_csv, combined_rows=combined_rows)
        csv.write_summary_csv(summary_csv=summary_csv, summary_rows=summary_rows)

    combined_plot = None
    page_paths: List[Path] = []
    if args.export_detail_plots and peak_records:
        page_paths = plot.export_individual_subplot_pages(plots_dir=detail_plots_dir, channel_stem=channel_stem, peak_records=peak_records, rows=3, cols=2)
    if args.export_summary_plot and peak_records:
        combined_plot = combined_plot_path
        plot.export_combined_plot(out_png=combined_plot, combined_rows=combined_rows, peak_records=peak_records)

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
            plot.open_file_with_default_app(combined_plot)
    elif args.show_summary_plot:
        print("  Summary plot was not opened because --export-summary-plot was not enabled or no peaks were processed.")

    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

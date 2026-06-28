import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

PeakRecord = Dict[str, object]

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




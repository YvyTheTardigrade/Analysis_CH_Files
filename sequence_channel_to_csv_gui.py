"""Tkinter GUI for sequence_channel_to_csv_and_combined_v20.py."""

#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
import threading
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CLI_SCRIPT = SCRIPT_DIR / "sequence_channel_to_csv_and_combined_v20.py"


def be_u32(b: bytes, off: int) -> int:
    return int.from_bytes(b[off:off + 4], byteorder="big", signed=False)


def be_i16(b: bytes, off: int) -> int:
    return int.from_bytes(b[off:off + 2], byteorder="big", signed=True)


def extract_folder_order(name: str) -> int | None:
    """Return the numeric order only for folders beginning with exactly 3 digits and a dash."""
    match = __import__("re").match(r"^(\d{3})-", name)
    return int(match.group(1)) if match else None


def extract_previous_data_name(name: str) -> str | None:
    """Extract the data name after the final dash and before .D, without amino-acid restrictions."""
    match = __import__("re").match(r"^\d{3}-.+-(.+)\.D$", name)
    if not match:
        return None
    data_name = match.group(1).strip()
    return data_name or None


def is_pip_folder(name: str) -> bool:
    """Target folders must start with 3 digits and contain -Pip_."""
    return bool(__import__("re").match(r"^\d{3}-", name)) and "-Pip_" in name and name.endswith(".D")


def build_pip_targets(parent_dir: Path) -> list[tuple[Path, str, str]]:
    folders: list[tuple[int, Path]] = []
    for p in parent_dir.iterdir():
        if p.is_dir():
            order = extract_folder_order(p.name)
            if order is not None:
                folders.append((order, p))
    folders.sort(key=lambda x: x[0])

    targets: list[tuple[Path, str, str]] = []
    previous_data_name: str | None = None
    for _order, folder in folders:
        if is_pip_folder(folder.name):
            if previous_data_name:
                folder_number = folder.name.split("-", 1)[0]
                targets.append((folder, folder_number, previous_data_name))
            continue

        data_name = extract_previous_data_name(folder.name)
        if data_name:
            previous_data_name = data_name
    return targets


def decode_signal_from_offset(data: bytes, start: int) -> list[int]:
    values: list[int] = []
    cursor = start
    current = 0

    while cursor + 2 <= len(data):
        delta = be_i16(data, cursor)
        cursor += 2

        if delta == -32768:
            if cursor + 4 > len(data):
                break
            hi_word = int.from_bytes(data[cursor:cursor + 2], "big", signed=False)
            lo_word = int.from_bytes(data[cursor + 2:cursor + 4], "big", signed=False)
            cursor += 4
            value = (hi_word << 16) | lo_word
            if value >= 0x80000000:
                value -= 0x100000000
            current = value
        else:
            current += delta
        values.append(current)

    return values


def extract_header(data: bytes) -> dict[str, int]:
    return {
        "first_time_ms": be_u32(data, 0x011A),
        "last_time_ms": be_u32(data, 0x011E),
        "header_size_bytes": 0x1800,
    }


def deduce_output_point_count(first_time_ms: int | None, last_time_ms: int | None, decoded_point_count: int) -> int | None:
    if first_time_ms is None or last_time_ms is None or decoded_point_count <= 0:
        return None
    if decoded_point_count == 1:
        return 1
    delta_ms = (last_time_ms - first_time_ms) / (decoded_point_count - 1)
    if delta_ms == 0:
        return decoded_point_count
    return max(int(round((last_time_ms - first_time_ms) / delta_ms)) + 1, 1)


class App:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("Sequence Channel to CSV - GUI v20 (English)")
        self.root.geometry("980x760")

        self.process: subprocess.Popen[str] | None = None

        self.parent_dir_var = tk.StringVar()
        self.channel_file_var = tk.StringVar(value="DAD1D.ch")
        self.output_dir_var = tk.StringVar()
        self.preview_file_var = tk.StringVar()
        self.combined_name_var = tk.StringVar(value="combined.csv")
        self.summary_name_var = tk.StringVar(value="summary.csv")
        self.points_var = tk.StringVar(value="Auto from .ch header")
        self.header_info_var = tk.StringVar(value="Select a parent directory to preview header-derived values.")
        self.export_csv_var = tk.BooleanVar(value=False)
        self.export_summary_plot_var = tk.BooleanVar(value=True)
        self.export_detail_plots_var = tk.BooleanVar(value=False)
        self.show_summary_plot_var = tk.BooleanVar(value=False)
        self.cli_script_var = tk.StringVar(value=str(DEFAULT_CLI_SCRIPT))

        self._build_ui()
        self.root.protocol("WM_DELETE_WINDOW", self.close_window)
        self._bind_updates()
        self.refresh_preview_file_list()
        self._set_header_preview_text(self.header_info_var.get())

    def _build_ui(self) -> None:
        container = ttk.Frame(self.root, padding=12)
        container.pack(fill="both", expand=True)
        container.columnconfigure(1, weight=1)

        title = ttk.Label(
            container,
            text="Sequence Channel to CSV - GUI v20 (English)",
            font=("TkDefaultFont", 14, "bold"),
        )
        title.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 6))

        subtitle = ttk.Label(
            container,
            text=(
                "v20 selects target folders that start with three digits and contain -Pip_, "
                "then labels each Pip folder from the immediately preceding numbered data folder. "
                "The label name is not restricted to amino acids, so names like Dawson link, t, LMt, and Lys2 are accepted."
            ),
            wraplength=920,
            justify="left",
        )
        subtitle.grid(row=1, column=0, columnspan=3, sticky="w", pady=(0, 10))

        row = 2
        row = self._add_file_row(container, row, "CLI script", self.cli_script_var, self._browse_cli_script, file_mode=True)
        row = self._add_file_row(container, row, "Parent input directory", self.parent_dir_var, self._browse_parent_dir)
        row = self._add_labeled_entry(container, row, "Channel filename", self.channel_file_var)

        ttk.Label(container, text="Selected .ch preview file").grid(row=row, column=0, sticky="w", pady=4, padx=(0, 8))
        preview_row = ttk.Frame(container)
        preview_row.grid(row=row, column=1, columnspan=2, sticky="ew", pady=4)
        preview_row.columnconfigure(0, weight=1)
        self.preview_file_combo = ttk.Combobox(preview_row, textvariable=self.preview_file_var, state="normal")
        self.preview_file_combo.grid(row=0, column=0, sticky="ew")
        self.preview_file_combo.bind("<<ComboboxSelected>>", lambda _event: self.refresh_header_preview())
        ttk.Button(preview_row, text="Browse file...", command=self._browse_preview_file).grid(row=0, column=1, padx=(8, 0))
        ttk.Button(preview_row, text="Refresh list", command=self.refresh_preview_file_list).grid(row=0, column=2, padx=(8, 0))
        row += 1

        ttk.Label(container, text="Header preview").grid(row=row, column=0, sticky="nw", pady=4, padx=(0, 8))
        self.header_preview_box = ScrolledText(container, height=3, wrap="word")
        self.header_preview_box.grid(row=row, column=1, columnspan=2, sticky="ew", pady=4)
        self.header_preview_box.configure(state="disabled")
        row += 1

        row = self._add_file_row(container, row, "Output directory", self.output_dir_var, self._browse_output_dir)
        row = self._add_labeled_entry(container, row, "Combined CSV file name", self.combined_name_var)
        row = self._add_labeled_entry(container, row, "Summary CSV file name", self.summary_name_var)

        ttk.Label(container, text="Output point count").grid(row=row, column=0, sticky="w", pady=4, padx=(0, 8))
        ttk.Entry(container, textvariable=self.points_var, state="readonly").grid(row=row, column=1, sticky="ew", pady=4)
        ttk.Button(container, text="Refresh from header", command=self.refresh_header_preview).grid(row=row, column=2, sticky="ew", pady=4, padx=(8, 0))
        row += 1

        options_frame = ttk.LabelFrame(container, text="Outputs")
        options_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(8, 8))
        options_frame.columnconfigure(0, weight=1)
        ttk.Checkbutton(options_frame, text="CSV files", variable=self.export_csv_var, command=self.update_command_preview).grid(row=0, column=0, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(options_frame, text="Summary plot (plots/summary/combined_plot.png)", variable=self.export_summary_plot_var, command=self.update_command_preview).grid(row=1, column=0, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(options_frame, text="All detail plots (plots/detail)", variable=self.export_detail_plots_var, command=self.update_command_preview).grid(row=2, column=0, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(options_frame, text="Show summary plot at end", variable=self.show_summary_plot_var, command=self.update_command_preview).grid(row=3, column=0, sticky="w", padx=8, pady=4)
        row += 1

        button_row = ttk.Frame(container)
        button_row.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(10, 10))
        button_row.columnconfigure(0, weight=0)
        button_row.columnconfigure(1, weight=0)
        button_row.columnconfigure(2, weight=0)
        button_row.columnconfigure(3, weight=1)

        self.run_button = ttk.Button(button_row, text="Run", command=self.run_program)
        self.run_button.grid(row=0, column=0, padx=(0, 8))

        self.stop_button = ttk.Button(button_row, text="Stop", command=self.stop_program, state="disabled")
        self.stop_button.grid(row=0, column=1, padx=(0, 8))

        self.close_button = ttk.Button(button_row, text="Close", command=self.close_window)
        self.close_button.grid(row=0, column=2, padx=(0, 8))

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(button_row, textvariable=self.status_var).grid(row=0, column=3, sticky="w")
        row += 1

        ttk.Label(container, text="Command preview").grid(row=row, column=0, columnspan=3, sticky="w")
        row += 1
        self.command_preview = ScrolledText(container, height=5, wrap="word")
        self.command_preview.grid(row=row, column=0, columnspan=3, sticky="nsew", pady=(0, 10))
        self.command_preview.configure(state="disabled")
        row += 1

        ttk.Label(container, text="Program output").grid(row=row, column=0, columnspan=3, sticky="w")
        row += 1
        self.output_box = ScrolledText(container, height=18, wrap="word")
        self.output_box.grid(row=row, column=0, columnspan=3, sticky="nsew")
        row += 1

        container.rowconfigure(row - 1, weight=1)
        self.update_command_preview()

    def _bind_updates(self) -> None:
        for variable in (
            self.cli_script_var,
            self.parent_dir_var,
            self.channel_file_var,
            self.output_dir_var,
            self.combined_name_var,
            self.summary_name_var,
        ):
            variable.trace_add("write", lambda *_: self._on_inputs_changed())
        self.export_csv_var.trace_add("write", lambda *_: self.update_command_preview())
        self.export_summary_plot_var.trace_add("write", lambda *_: self.update_command_preview())
        self.export_detail_plots_var.trace_add("write", lambda *_: self.update_command_preview())
        self.show_summary_plot_var.trace_add("write", lambda *_: self.update_command_preview())
        self.preview_file_var.trace_add("write", lambda *_: self.update_command_preview())

    def _on_inputs_changed(self) -> None:
        parent_text = self.parent_dir_var.get().strip()
        if parent_text and not self.output_dir_var.get().strip():
            self.output_dir_var.set(str(Path(parent_text) / "Analysis"))
            return
        self.refresh_preview_file_list()
        self.update_command_preview()

    def _add_labeled_entry(self, parent: ttk.Frame, row: int, label: str, variable: tk.StringVar) -> int:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4, padx=(0, 8))
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, columnspan=2, sticky="ew", pady=4)
        return row + 1

    def _add_file_row(self, parent: ttk.Frame, row: int, label: str, variable: tk.StringVar, browse_command, file_mode: bool = False) -> int:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4, padx=(0, 8))
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, sticky="ew", pady=4)
        ttk.Button(parent, text="Browse...", command=browse_command).grid(row=row, column=2, sticky="ew", pady=4, padx=(8, 0))
        return row + 1

    def _browse_cli_script(self) -> None:
        path = filedialog.askopenfilename(
            title="Select the CLI Python script",
            filetypes=[("Python files", "*.py"), ("All files", "*.*")],
            initialdir=str(SCRIPT_DIR),
        )
        if path:
            self.cli_script_var.set(path)

    def _browse_parent_dir(self) -> None:
        path = filedialog.askdirectory(title="Select the parent input directory")
        if path:
            self.parent_dir_var.set(path)
            self.output_dir_var.set(str(Path(path) / "Analysis"))
            self.refresh_header_preview()

    def _browse_output_dir(self) -> None:
        path = filedialog.askdirectory(title="Select the output directory")
        if path:
            self.output_dir_var.set(path)

    def _browse_preview_file(self) -> None:
        initial_dir = self.parent_dir_var.get().strip() or str(SCRIPT_DIR)
        initial_file = self.channel_file_var.get().strip()
        path = filedialog.askopenfilename(
            title="Select the .ch file to preview",
            initialdir=initial_dir,
            initialfile=initial_file,
            filetypes=[("Channel files", "*.ch"), ("All files", "*.*")],
        )
        if path:
            values = list(self.preview_file_combo.cget("values"))
            if path not in values:
                values.append(path)
                self.preview_file_combo["values"] = values
            self.preview_file_var.set(path)
            self.refresh_header_preview()

    def _list_preview_ch_files(self) -> list[Path]:
        parent_text = self.parent_dir_var.get().strip()
        channel_name = self.channel_file_var.get().strip()
        if not parent_text or not channel_name:
            return []
        parent_dir = Path(parent_text)
        if not parent_dir.is_dir():
            return []

        matches: list[Path] = []
        for pip_folder, _folder_number, _data_name in build_pip_targets(parent_dir):
            candidate = pip_folder / channel_name
            if candidate.is_file():
                matches.append(candidate)
        return matches

    def refresh_preview_file_list(self) -> None:
        matches = self._list_preview_ch_files()
        display_values = [str(path.relative_to(path.parents[1])) if len(path.parents) >= 2 else path.name for path in matches]

        current_value = self.preview_file_var.get().strip()
        if current_value and Path(current_value).is_file() and current_value not in display_values:
            display_values.append(current_value)

        self.preview_file_combo["values"] = display_values

        if display_values:
            if current_value not in display_values:
                self.preview_file_var.set(display_values[0])
        else:
            self.preview_file_var.set(current_value if Path(current_value).is_file() else "")

    def _find_preview_ch_file(self) -> Path | None:
        selected_display = self.preview_file_var.get().strip()
        if selected_display:
            selected_path = Path(selected_display)
            if selected_path.is_file():
                return selected_path

        matches = self._list_preview_ch_files()
        if not matches:
            return None

        if not selected_display:
            return matches[0]

        for path in matches:
            display = str(path.relative_to(path.parents[1])) if len(path.parents) >= 2 else path.name
            if display == selected_display:
                return path
        return matches[0]

    def refresh_header_preview(self) -> None:
        ch_file = self._find_preview_ch_file()
        if ch_file is None:
            self.points_var.set("Auto from .ch header")
            self.header_info_var.set("No preview .ch file found yet. Select a valid parent directory containing numbered -Pip_ .D folders.")
            self._set_header_preview_text(self.header_info_var.get())
            self.update_command_preview()
            return

        try:
            data = ch_file.read_bytes()
            header = extract_header(data)
            decoded = decode_signal_from_offset(data, int(header["header_size_bytes"]))
            decoded_count = len(decoded)
            first_time_ms = header.get("first_time_ms")
            last_time_ms = header.get("last_time_ms")
            point_count = deduce_output_point_count(first_time_ms, last_time_ms, decoded_count)

            if decoded_count > 1 and first_time_ms is not None and last_time_ms is not None:
                delta_ms = (last_time_ms - first_time_ms) / (decoded_count - 1)
                second_time_ms = first_time_ms + delta_ms
                info = (
                    f"Preview file: {ch_file.name} | first_time_ms={first_time_ms} | "
                    f"second_time_ms={second_time_ms:.6f} | last_time_ms={last_time_ms} | "
                    f"header_size_bytes={header['header_size_bytes']} | decoded points={decoded_count}"
                )
            else:
                info = (
                    f"Preview file: {ch_file.name} | first_time_ms={first_time_ms} | "
                    f"last_time_ms={last_time_ms} | header_size_bytes={header['header_size_bytes']} | "
                    f"decoded points={decoded_count}"
                )

            self.points_var.set(str(point_count) if point_count is not None else str(decoded_count))
            self.header_info_var.set(info)
            self._set_header_preview_text(info)
        except Exception as exc:
            self.points_var.set("Auto from .ch header")
            self.header_info_var.set(f"Could not read header preview: {exc}")
            self._set_header_preview_text(self.header_info_var.get())

        self.update_command_preview()

    def _build_command(self) -> list[str]:
        cmd = [
            sys.executable,
            self.cli_script_var.get().strip(),
            self.parent_dir_var.get().strip(),
            self.channel_file_var.get().strip(),
            "--output-dir",
            self.output_dir_var.get().strip(),
            "--combined-name",
            self.combined_name_var.get().strip(),
            "--summary-name",
            self.summary_name_var.get().strip(),
        ]

        points_text = self.points_var.get().strip()
        if points_text and points_text.isdigit():
            cmd.extend(["--points", points_text])

        if self.export_csv_var.get():
            cmd.append("--export-csv")
        if self.export_summary_plot_var.get():
            cmd.append("--export-summary-plot")
        if self.export_detail_plots_var.get():
            cmd.append("--export-detail-plots")
        if self.show_summary_plot_var.get():
            cmd.append("--show-summary-plot")

        return cmd

    def _set_header_preview_text(self, text: str) -> None:
        self.header_preview_box.configure(state="normal")
        self.header_preview_box.delete("1.0", tk.END)
        self.header_preview_box.insert(tk.END, text)
        self.header_preview_box.configure(state="disabled")

    def update_command_preview(self) -> None:
        cmd = self._build_command()
        preview = " ".join(self._quote_if_needed(part) for part in cmd)
        self.command_preview.configure(state="normal")
        self.command_preview.delete("1.0", tk.END)
        self.command_preview.insert(tk.END, preview)
        self.command_preview.configure(state="disabled")

    @staticmethod
    def _quote_if_needed(value: str) -> str:
        if not value:
            return '""'
        if any(ch.isspace() for ch in value) or any(ch in value for ch in '"\''):
            return f'"{value}"'
        return value

    def validate_inputs(self) -> bool:
        cli_script = Path(self.cli_script_var.get().strip())
        if not cli_script.is_file():
            messagebox.showerror("Missing CLI script", "Please select a valid CLI Python script.")
            return False

        parent_dir = Path(self.parent_dir_var.get().strip())
        if not parent_dir.is_dir():
            messagebox.showerror("Missing parent directory", "Please select a valid parent input directory.")
            return False

        if not self.channel_file_var.get().strip():
            messagebox.showerror("Missing channel filename", "Please enter the channel filename.")
            return False

        output_dir_text = self.output_dir_var.get().strip()
        if not output_dir_text:
            messagebox.showerror("Missing output directory", "Please select an output directory.")
            return False

        points_text = self.points_var.get().strip()
        if points_text and points_text != "Auto from .ch header" and not points_text.isdigit():
            messagebox.showerror("Invalid point count", "Output point count must be numeric if shown.")
            return False

        return True

    def run_program(self) -> None:
        if not self.validate_inputs():
            return
        if self.process is not None:
            messagebox.showinfo("Already running", "The program is already running.")
            return

        self.refresh_header_preview()
        cmd = self._build_command()
        self.output_box.delete("1.0", tk.END)
        self._append_output("Running command:\n")
        self._append_output(" ".join(self._quote_if_needed(part) for part in cmd) + "\n\n")

        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        self.status_var.set("Running...")

        thread = threading.Thread(target=self._run_subprocess, args=(cmd,), daemon=True)
        thread.start()

    def _run_subprocess(self, cmd: list[str]) -> None:
        try:
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            assert self.process.stdout is not None
            for line in self.process.stdout:
                self.root.after(0, self._append_output, line)
            return_code = self.process.wait()
            self.root.after(0, self._finish_run, return_code)
        except Exception as exc:
            self.root.after(0, self._append_output, f"\nERROR: {exc}\n")
            self.root.after(0, self._finish_run, -1)

    def stop_program(self) -> None:
        if self.process is not None:
            self.process.terminate()
            self._append_output("\nTermination requested by the user.\n")
            self.status_var.set("Stopping...")
            self.stop_button.configure(state="disabled")

    def close_window(self) -> None:
        if self.process is not None:
            try:
                self.process.terminate()
            except Exception:
                pass
            self.process = None
        self.root.destroy()

    def _finish_run(self, return_code: int) -> None:
        self.process = None
        self.run_button.configure(state="normal")
        self.stop_button.configure(state="disabled")
        if return_code == 0:
            self.status_var.set("Finished successfully")
            self._append_output("\nFinished successfully.\n")
        elif return_code == -15:
            self.status_var.set("Stopped")
            self._append_output("\nProgram stopped.\n")
        else:
            self.status_var.set(f"Finished with exit code {return_code}")
            self._append_output(f"\nProgram finished with exit code {return_code}.\n")

    def _append_output(self, text: str) -> None:
        self.output_box.insert(tk.END, text)
        self.output_box.see(tk.END)


def main() -> int:
    root = tk.Tk()
    app = App(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

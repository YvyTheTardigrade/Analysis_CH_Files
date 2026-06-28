from typing import Any, Dict, List, Optional, Set, Tuple
import re
from pathlib import Path




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



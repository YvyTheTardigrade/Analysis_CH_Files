# PepSynth Analyzer

A tool for analyzing peptide synthesis sequences from Agilent `.ch` channel files.

---

## Requirements

All Python dependencies are listed in `requirements.txt`.
It is recommended to install [`make`](https://www.gnu.org/software/make/) to use the predefined shortcuts in `Makefile`.

---

## Setup

### Using pipenv

```bash
pipenv install -r requirements.txt
pipenv shell
```

### Using venv

```bash
python -m venv .venv
source .venv/bin/activate        # Linux / macOS
.venv\Scripts\activate           # Windows
pip install -r requirements.txt
```

---

## Run

| Method | Command |
|---|---|
| With make | `make run` |
| Without make | `PYTHONPATH=. python src/Main.py` |

---

## Test

Tests are written with [`pytest`](https://pytest.org).

| Method | Command |
|---|---|
| With make | `make test` |
| Without make | `PYTHONPATH=. pytest` |

### Adding new tests

Add a file named `test_*.py` in the `tests/` folder. Each test function must follow this pattern:

```python
def test_my_feature():
    result = my_function(input)
    assert result == expected_output
```

---

## Build

Compiles the project into a standalone executable using [PyInstaller](https://pyinstaller.org).

| Method | Command |
|---|---|
| With make | `make build` |
| Without make | See below |

```bash
pyinstaller --onefile --windowed \
    --name "PepSynth_Analyzer" \
    --paths=src \
    --distpath "./build/executable/" \
    --workpath "./build/tmp/" \
    --specpath "./build/spec/" \
    src/Main.py
```

The resulting executable is placed in `build/executable/`.

### Running the executable

| Method | Command |
|---|---|
| With make | `make run_build` |
| Without make | `./build/executable/PepSynth_Analyzer` |

You can also double-click the file in your file explorer.

> [!IMPORTANT]
> The executable is platform-specific. A build made on Linux will only run on Linux, and similarly for Windows and macOS. Build on the target platform.

---

## Deploy

Send the file located in `build/executable/` — that is all that is needed.
No Python installation is required on the target machine.

"""
write_ch.py
-----------
Génère un fichier Agilent .ch (format "other" : UV/CAD/ELSD) à partir d'une
série temporelle Python.

Référence de format :
  https://rainbow-api.readthedocs.io/en/latest/agilent/ch_other.html

Utilisation rapide
------------------
    from write_ch import write_ch
    import numpy as np

    times   = np.linspace(0, 10, 1000)          # minutes
    values  = np.sin(2 * np.pi * times) * 500   # mAU

    write_ch(
        path        = "DAD1A.ch",
        times_min   = times,
        values      = values,
        units       = "mAU",
        signal      = "DAD1A, Sig=215.0,16.0 Ref=off",
        method      = "MY_METHOD.M",
        instrument  = "My ChemStation",
        date        = "28-Jun-26, 12:00:00",
        notebook    = "MY_NOTEBOOK",
    )
"""

import struct
import datetime
from typing import Optional, Sequence
import numpy as np


# ---------------------------------------------------------------------------
# Constantes de format
# ---------------------------------------------------------------------------

HEADER_SIZE   = 0x1800          # 6144 octets — taille fixe du header
SEGMENT_LABEL = 0x10            # octet fixe qui commence chaque segment (= 16)
MAX_PER_SEG   = 255             # max de valeurs par segment (tient sur 1 octet)
DELTA_FLAG    = b'\x80\x00'     # marqueur "valeur absolue" (–0x8000 en signé)


# ---------------------------------------------------------------------------
# Helpers d'encodage du format "string with null separators"
# ---------------------------------------------------------------------------

def _encode_string(s: str) -> bytes:
    """
    Encode une chaîne selon le format Agilent :
      1 octet longueur  +  chaque caractère suivi d'un octet nul.
    Ex : "hello" → b'\\x05h\\x00e\\x00l\\x00l\\x00o\\x00'
    """
    encoded = b""
    for ch in s:
        encoded += ch.encode("latin-1") + b"\x00"
    return bytes([len(s)]) + encoded


def _write_field(buf: bytearray, offset: int, s: str, max_bytes: int = 512) -> None:
    """
    Écrit une chaîne encodée dans `buf` à `offset`, tronquée si nécessaire.
    """
    raw = _encode_string(s)
    raw = raw[:max_bytes]
    buf[offset: offset + len(raw)] = raw


# ---------------------------------------------------------------------------
# Encodage du corps de données (delta encoding)
# ---------------------------------------------------------------------------

def _encode_body(raw_integers: Sequence[int]) -> bytes:
    """
    Encode la liste d'entiers `raw_integers` en segments Agilent.

    Chaque valeur est soit :
      - absolue  : 0x8000 (2 octets) + valeur int32 big-endian (4 octets) → 6 oct.
      - delta    : delta int16 big-endian par rapport à la valeur précédente → 2 oct.

    La première valeur d'un segment (ou tout saut qui ne tient pas en int16)
    est stockée de façon absolue.
    """
    body = bytearray()
    i = 0
    n = len(raw_integers)

    while i < n:
        # Déterminer la taille du segment (≤ MAX_PER_SEG valeurs)
        seg_size = min(MAX_PER_SEG, n - i)
        seg_values = raw_integers[i: i + seg_size]

        # En-tête du segment
        seg_bytes = bytearray()
        seg_bytes.append(SEGMENT_LABEL)          # 0x10
        seg_bytes.append(seg_size)               # nb de valeurs

        current_abs = None  # dernière valeur absolue mémorisée

        for val in seg_values:
            if current_abs is None:
                # Toujours absolue pour la première valeur du segment
                seg_bytes += DELTA_FLAG
                seg_bytes += struct.pack(">i", val)
                current_abs = val
            else:
                delta = val - current_abs
                if -32768 <= delta <= 32767 and delta != -32768:
                    # Delta tient dans un int16 (et pas confondu avec le flag)
                    seg_bytes += struct.pack(">h", delta)
                    current_abs += delta
                else:
                    # Valeur absolue
                    seg_bytes += DELTA_FLAG
                    seg_bytes += struct.pack(">i", val)
                    current_abs = val

        body += seg_bytes
        i += seg_size

    body += b"\x00\x00"   # 2 octets nuls de fin de fichier
    return bytes(body)


# ---------------------------------------------------------------------------
# Fonction principale
# ---------------------------------------------------------------------------

def write_ch(
    path: str,
    times_min: Sequence[float],
    values: Sequence[float],
    *,
    units: str = "mAU",
    signal: str = "DAD1A, Sig=215.0,16.0 Ref=off",
    method: str = "METHOD.M",
    instrument: str = "Agilent ChemStation",
    date: Optional[str] = None,
    notebook: str = "",
    parent_dir: str = "",
    scaling_factor: Optional[float] = None,
) -> None:
    """
    Écrit un fichier Agilent .ch (format UV/CAD/ELSD).

    Paramètres
    ----------
    path          : chemin du fichier de sortie (ex. "DAD1A.ch")
    times_min     : tableau de temps en **minutes**
    values        : tableau des valeurs du signal (même longueur que times_min)
    units         : unités (ex. "mAU", "pA", "counts")
    signal        : description du signal (ex. "DAD1A, Sig=215.0,16.0 Ref=off")
    method        : nom de la méthode (ex. "MY_METHOD.M")
    instrument    : nom de l'instrument
    date          : date/heure au format "DD-Mon-YY, HH:MM:SS" ; auto si None
    notebook      : nom du notebook
    parent_dir    : répertoire parent (métadonnée interne)
    scaling_factor: si None, calculé automatiquement pour maximiser la précision
    """
    times_min = np.asarray(times_min, dtype=float)
    values    = np.asarray(values,    dtype=float)

    if len(times_min) != len(values):
        raise ValueError("times_min et values doivent avoir la même longueur.")
    if len(times_min) == 0:
        raise ValueError("La série temporelle est vide.")

    # --- Temps de rétention en ms (entiers non-signés, big-endian) -----------
    t_start_ms = int(round(times_min[0]  * 60_000))
    t_end_ms   = int(round(times_min[-1] * 60_000))

    # --- Facteur d'échelle ---------------------------------------------------
    # Les valeurs stockées sont des int32 ; on choisit le facteur pour que
    # la valeur maximale (en absolu) occupe ~80 % de la plage int32.
    if scaling_factor is None:
        max_abs = np.max(np.abs(values))
        if max_abs == 0:
            scaling_factor = 1.0
        else:
            scaling_factor = max_abs / (0.8 * 2**31)

    raw_integers = np.round(values / scaling_factor).astype(np.int64)

    # Vérification des bornes int32
    if raw_integers.max() > 2**31 - 1 or raw_integers.min() < -(2**31):
        raise OverflowError(
            "Les valeurs dépassent la plage int32 avec ce scaling_factor. "
            "Essayez un scaling_factor plus grand."
        )
    raw_integers = raw_integers.astype(np.int32).tolist()

    # --- Date ----------------------------------------------------------------
    if date is None:
        date = datetime.datetime.now().strftime("%d-%b-%y, %H:%M:%S")

    # =========================================================================
    # Construction du header (6144 octets, tout à 0x00 par défaut)
    # =========================================================================
    header = bytearray(HEADER_SIZE)

    # Champ 0x146 : type de fichier (nombre) → "130"
    _write_field(header, 0x146, "130")

    # Champ 0x15B : type de fichier (nom) → "LC DATA FILE"
    _write_field(header, 0x15B, "LC DATA FILE")

    # Champ 0x35A : nom du notebook
    _write_field(header, 0x35A, notebook, max_bytes=256)

    # Champ 0x758 : répertoire parent
    _write_field(header, 0x758, parent_dir, max_bytes=256)

    # Champ 0x957 : date
    _write_field(header, 0x957, date, max_bytes=128)

    # Champ 0xA0E : méthode
    _write_field(header, 0xA0E, method, max_bytes=256)

    # Champ 0xC11 : instrument
    _write_field(header, 0xC11, instrument, max_bytes=256)

    # Champ 0x104C : unités
    _write_field(header, 0x104C, units, max_bytes=64)

    # Champ 0x1075 : signal
    _write_field(header, 0x1075, signal, max_bytes=256)

    # Champ 0x11A : premier temps de rétention (ms) — uint32 big-endian
    struct.pack_into(">I", header, 0x11A, t_start_ms & 0xFFFFFFFF)

    # Champ 0x11E : dernier temps de rétention (ms) — uint32 big-endian
    struct.pack_into(">I", header, 0x11E, t_end_ms & 0xFFFFFFFF)

    # Champ 0x127C : facteur d'échelle — double big-endian
    struct.pack_into(">d", header, 0x127C, scaling_factor)

    # =========================================================================
    # Corps de données
    # =========================================================================
    body = _encode_body(raw_integers)

    # =========================================================================
    # Écriture du fichier
    # =========================================================================
    with open(path, "wb") as f:
        f.write(header)
        f.write(body)

    print(
        f"✅  Fichier écrit : {path}\n"
        f"   Points        : {len(values)}\n"
        f"   Temps         : {times_min[0]:.4f} – {times_min[-1]:.4f} min\n"
        f"   Scaling factor: {scaling_factor:.6e}\n"
        f"   Taille totale : {HEADER_SIZE + len(body)} octets"
    )

    return locals()


# ---------------------------------------------------------------------------
# Exemple d'utilisation (lancez ce fichier directement pour tester)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import numpy as np

    # --- Série temporelle de démonstration -----------------------------------
    # Chromatogramme synthétique : 3 pics gaussiens sur 12 minutes
    t = np.linspace(0, 12, 2400)   # 2400 points, 1 point toutes les 0.3 s

    def gaussian(t, mu, sigma, amplitude):
        return amplitude * np.exp(-0.5 * ((t - mu) / sigma) ** 2)

    signal = (
        gaussian(t, 3.0, 0.15, 850)   # pic 1 à 3 min,  hauteur 850 mAU
      + gaussian(t, 5.5, 0.20, 1200)  # pic 2 à 5.5 min, hauteur 1200 mAU
      + gaussian(t, 9.0, 0.10, 400)   # pic 3 à 9 min,  hauteur 400 mAU
      + np.random.normal(0, 2, len(t)) # bruit de fond
    )

    write_ch(
        path        = "DAD1A_demo.ch",
        times_min   = t,
        values      = signal,
        units       = "mAU",
        signal      = "DAD1A, Sig=215.0,16.0 Ref=off",
        method      = "DEMO_METHOD.M",
        instrument  = "Agilent 1260 Infinity II",
        date        = "28-Jun-26, 12:00:00",
        notebook    = "DEMO_NOTEBOOK",
        parent_dir  = "demo",
    )

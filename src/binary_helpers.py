
import struct




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
    return struct.unpack(">d", b[off:off + 8])[0] # ">" = big-endian  , "d" = double precision 


def be_i32_from_words(hi_word: int, lo_word: int) -> int:
    value = (hi_word << 16) | lo_word
    if value >= 0x80000000:
        value -= 0x100000000
    return value





# ---------------------------------------------------------------------
# Header reader derived from ch_to_csv_two_csvs.py
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




import numpy as np

from owl_lip.params import INT_NBYTES, MATRIX_ENTRY_TYPE
from owl_lip.typing import GroupElementT, MatrixEntryT, MatrixT, SetElementT


def encode_int_fixed(x: MatrixEntryT, nbytes: int = INT_NBYTES) -> bytes:
    return int(x).to_bytes(nbytes, byteorder="big", signed=False)


def decode_int_fixed(b: bytes) -> MatrixEntryT:
    return MATRIX_ENTRY_TYPE(int.from_bytes(b, byteorder="big", signed=False))


def encode_matrix(b: MatrixT, nbytes: int = INT_NBYTES) -> bytes:
    b = np.array(b, dtype=MATRIX_ENTRY_TYPE)
    if b.shape != (2, 2):
        raise ValueError("Group matrix must be 2x2")

    out = bytearray()
    for i in range(2):
        for j in range(2):
            out += encode_int_fixed(b[i, j], nbytes)

    return bytes(out)


def decode_matrix(data: bytes, nbytes: int = INT_NBYTES) -> MatrixT:
    if len(data) != 4 * nbytes:
        raise ValueError("Invalid group encoding length")

    entries = []
    for i in range(0, len(data), nbytes):
        entries.append(decode_int_fixed(data[i : i + nbytes]))

    return np.array(entries, dtype=MATRIX_ENTRY_TYPE).reshape(2, 2)


def encode_set(s: SetElementT) -> bytes:
    return encode_matrix(s)


def decode_set(bs: bytes) -> SetElementT:
    return decode_matrix(bs)


def encode_group(g: GroupElementT) -> bytes:
    return encode_matrix(g)


def decode_group(bs: bytes) -> GroupElementT:
    return decode_matrix(bs)

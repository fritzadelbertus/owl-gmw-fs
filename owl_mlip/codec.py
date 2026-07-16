from owl_mlip.params import LOGN
from owl_mlip.params import NTT_MODULUS as Q
from owl_mlip.typing import GroupElementT, SetElementT


def encode_poly_matrix(mat: list[list[list[int]]]) -> bytes:
    coeffs = []

    for row in mat:
        for poly in row:
            for c in poly:
                coeffs.append(c % Q)

    encoded = b""
    for c in coeffs:
        encoded += c.to_bytes(4, byteorder="big")

    return encoded


def decode_poly_matrix(encoded: bytes) -> list[list[list[int]]]:
    coeffs = []
    n = 2**LOGN

    for i in range(0, len(encoded), 4):
        c = int.from_bytes(encoded[i : i + 4], byteorder="big")
        coeffs.append(c % Q)

    # coeffs_array = np.array(coeffs, dtype=int)
    zero_entry = [0 for _ in range(2**LOGN)]
    mat = [[zero_entry, zero_entry], [zero_entry, zero_entry]]

    idx = 0
    for i in range(2):
        for j in range(2):
            mat[i][j] = coeffs[idx : idx + n]
            idx += n

    return mat


def encode_set(s: SetElementT) -> bytes:
    return encode_poly_matrix(s)


def decode_set(bs: bytes) -> SetElementT:
    return decode_poly_matrix(bs)


def encode_group(g: GroupElementT) -> bytes:
    return encode_poly_matrix(g)


def decode_group(bs: bytes) -> GroupElementT:
    return decode_poly_matrix(bs)

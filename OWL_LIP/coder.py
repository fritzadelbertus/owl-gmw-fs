import numpy as np

SIG_DIGITS = 12
INT_NBYTES = 16


def encode_int_fixed(x, nbytes=INT_NBYTES):
    return int(x).to_bytes(nbytes, byteorder="big", signed=True)


def decode_int_fixed(b):
    return int.from_bytes(b, byteorder="big", signed=True)


def encode_group_matrix(B, nbytes=INT_NBYTES):
    B = np.array(B, dtype=object)
    if B.shape != (2, 2):
        raise ValueError("Group matrix must be 2x2")

    out = bytearray()
    for i in range(2):
        for j in range(2):
            out += encode_int_fixed(B[i, j], nbytes)

    return bytes(out)


def decode_group_matrix(data, nbytes=INT_NBYTES):
    if len(data) != 4 * nbytes:
        raise ValueError("Invalid group encoding length")

    entries = []
    for i in range(0, len(data), nbytes):
        entries.append(decode_int_fixed(data[i:i+nbytes]))

    return np.array(entries, dtype=object).reshape(2, 2)


def encode_alg(x, nbytes=INT_NBYTES):
    out = bytearray()
    for coeff in x:
        out += int(coeff).to_bytes(nbytes, "big", signed=True)
    return bytes(out)

def decode_alg(data, nbytes=INT_NBYTES):
    if len(data) != 4 * nbytes:
        raise ValueError("Invalid algebraic number encoding")
    return tuple(
        int.from_bytes(data[i:i+nbytes], "big", signed=True)
        for i in range(0, 4*nbytes, nbytes)
    )

def encode_set_matrix(Q, nbytes=INT_NBYTES):
    out = bytearray()
    for i in range(2):
        for j in range(2):
            out += encode_alg(Q[i][j], nbytes)
    return bytes(out)

def decode_set_matrix(data, nbytes=INT_NBYTES):
    if len(data) != 2 * 2 * 4 * nbytes:
        raise ValueError("Invalid set matrix encoding")

    pos = 0
    Q = [[None, None], [None, None]]

    for i in range(2):
        for j in range(2):
            Q[i][j] = decode_alg(data[pos:pos+4*nbytes], nbytes)
            pos += 4*nbytes

    return Q
import numpy as np

def encode_poly_matrix(mat, q=8380417):
    """
    Encode a 2x2 polynomial matrix where each polynomial has length n.
    """
    coeffs = []
    
    for row in mat:
        for poly in row:
            for c in poly:
                coeffs.append(c % q)

    encoded = bytearray()
    for c in coeffs:
        encoded += c.to_bytes(4, byteorder='big')

    return encoded


def decode_poly_matrix(encoded, logn, q=8380417):
    """
    Decode a 2x2 polynomial matrix with polynomial length n.
    """
    coeffs = []
    n = 2**logn

    for i in range(0, len(encoded), 4):
        c = int.from_bytes(encoded[i:i+4], byteorder='big')
        coeffs.append(c % q)

    coeffs = np.array(coeffs, dtype=int)

    mat = np.zeros((2, 2, n), dtype=int)

    idx = 0
    for i in range(2):
        for j in range(2):
            mat[i, j] = coeffs[idx:idx+n]
            idx += n

    return mat
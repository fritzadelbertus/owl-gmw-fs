import numpy as np

from owl_lip.params import MATRIX_ENTRY_TYPE, P
from owl_lip.typing import MatrixT


def matmul_fp(a: MatrixT, b: MatrixT) -> MatrixT:
    result = np.empty((2, 2), dtype=MATRIX_ENTRY_TYPE)

    for i in range(2):
        for j in range(2):
            value = sum(int(a[i, k]) * int(b[k, j]) for k in range(2))
            result[i, j] = value % P

    return result

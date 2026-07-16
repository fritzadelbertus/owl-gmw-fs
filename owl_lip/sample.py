import numpy as np

from owl_lip.field import matmul_fp
from owl_lip.params import MATRIX_ENTRY_TYPE, ORBIT, L, P
from owl_lip.typing import GroupElementT, MatrixT, SetElementT
from rngcontext import RngContext


def canonical_pm(b: MatrixT) -> MatrixT:
    normalized = np.empty((2, 2), dtype=MATRIX_ENTRY_TYPE)

    for i in range(2):
        for j in range(2):
            normalized[i, j] = int(b[i, j]) % P

    for x_raw in normalized.flat:
        x = int(x_raw)

        if x != 0:
            if x <= P // 2:
                return normalized

            return np.array(
                [[(-int(normalized[i, j])) % P for j in range(2)] for i in range(2)],
                dtype=MATRIX_ENTRY_TYPE,
            )

    return normalized


def sample_k(rng: RngContext) -> int:
    return int(rng.random(1)[0] % 5) - 2


def sample_group(rng: RngContext) -> GroupElementT:
    b = np.eye(2, dtype=MATRIX_ENTRY_TYPE)
    moves = rng.random(L) & 1

    for move in moves[:-1]:
        k = sample_k(rng) % P
        e = (
            np.array([[1, k], [0, 1]], dtype=MATRIX_ENTRY_TYPE)
            if move == 0
            else np.array([[1, 0], [k, 1]], dtype=MATRIX_ENTRY_TYPE)
        )
        b = matmul_fp(b, e)

    if moves[-1] == 1:
        s = np.array([[0, 1], [1, 0]], dtype=MATRIX_ENTRY_TYPE)
        b = matmul_fp(b, s)

    return canonical_pm(b)


def sample_set(rng: RngContext) -> SetElementT:
    b = sample_group(rng)
    return matmul_fp(matmul_fp(b.T, ORBIT), b)

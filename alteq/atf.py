from ALTEQ.field import (
    multiplication_mod_p,
    reduction_mod_p,
    reduction_strict,
    set_inversion_mod_p,
)
from ALTEQ.params import LEN, PRIME, C, N


def compress_atf(atf_in, vec_size: int, nb_atf: int):
    catf_out = [0] * (LEN * vec_size)
    index = 0

    for i in range(N - 2):
        for j in range(i + 1, N - 1):
            for k in range(j + 1, N):
                for r in range(nb_atf):
                    catf_out[index] = reduction_strict(
                        atf_in[(i * N * N + j * N + k) * vec_size + r]
                    )
                    index += 1
    return catf_out


def decompress_atf(atf_in, vec_size: int):
    datf_out = [0] * (N * N * N * vec_size)

    index = 0

    for i in range(N - 2):
        for j in range(i + 1, N - 1):
            for k in range(j + 1, N):
                for r in range(vec_size):
                    datf_out[(i * N * N + j * N + k) * vec_size + r] = atf_in[index + r]
                index += vec_size

    return datf_out


def _mul_(bid, x, y, vec_size, buf, atf, col):
    for r in range(vec_size):
        buf[(bid) * vec_size + r] = multiplication_mod_p(
            atf[(x) * vec_size + r], col[(y) * vec_size + r]
        )
    return buf


def _madd_(bid, x, y, vec_size, buf, atf, col):
    for r in range(vec_size):
        buf[(bid) * vec_size + r] += multiplication_mod_p(
            atf[(x) * vec_size + r], col[(y) * vec_size + r]
        )
    return buf


def _msub_(bid, x, y, vec_size, buf, atf, col):
    for r in range(vec_size):
        buf[(bid) * vec_size + r] += multiplication_mod_p(
            atf[(x) * vec_size + r], PRIME - col[(y) * vec_size + r]
        )
    return buf


def _mod_red_(x, bid, vec_size, buf, atf):
    for r in range(vec_size):
        atf[(x) * vec_size + r] = reduction_mod_p(buf[(bid) * vec_size + r])
    return atf


def acting_on_atf(atf, col, j, vec_size):
    buf = [0] * (N * N * vec_size)

    for k in range(j):
        for m in range(j + 1, N):
            buf = _mul_(k * N + m, k * (N * N) + j * N + m, j, vec_size, buf, atf, col)

    for k in range(j - 1):
        for m in range(j + 1, N):
            for i in range(k + 1, j):
                buf = _madd_(
                    k * N + m, k * (N * N) + i * N + m, i, vec_size, buf, atf, col
                )

    for k in range(1, j):
        for m in range(j + 1, N):
            for i in range(k):
                buf = _msub_(
                    k * N + m, i * (N * N) + k * N + m, i, vec_size, buf, atf, col
                )

    for k in range(j):
        for m in range(j + 1, N - 1):
            for i in range(m + 1, N):
                buf = _msub_(
                    k * N + m, k * (N * N) + m * N + i, i, vec_size, buf, atf, col
                )

    for k in range(j):
        for m in range(j + 2, N):
            for i in range(j + 1, m):
                buf = _madd_(
                    k * N + m, k * (N * N) + i * N + m, i, vec_size, buf, atf, col
                )

    for k in range(j):
        for m in range(j + 1, N):
            atf = _mod_red_(k * (N * N) + j * N + m, k * N + m, vec_size, buf, atf)

    for k in range(j + 1, N - 1):
        for m in range(k + 1, N):
            buf = _mul_(k * N + m, j * (N * N) + k * N + m, j, vec_size, buf, atf, col)

    for k in range(j + 1, N - 1):
        for m in range(k + 1, N):
            for i in range(j):
                buf = _madd_(
                    k * N + m, i * (N * N) + k * N + m, i, vec_size, buf, atf, col
                )

    for k in range(j + 1, N - 2):
        for m in range(k + 1, N - 1):
            for i in range(m + 1, N):
                buf = _madd_(
                    k * N + m, k * (N * N) + m * N + i, i, vec_size, buf, atf, col
                )

    for k in range(j + 1, N - 2):
        for m in range(k + 2, N):
            for i in range(k + 1, m):
                buf = _msub_(
                    k * N + m, k * (N * N) + i * N + m, i, vec_size, buf, atf, col
                )

    for k in range(j + 2, N - 1):
        for m in range(k + 1, N):
            for i in range(j + 1, k):
                buf = _madd_(
                    k * N + m, i * (N * N) + k * N + m, i, vec_size, buf, atf, col
                )

    for k in range(j + 1, N - 1):
        for m in range(k + 1, N):
            atf = _mod_red_(j * (N * N) + k * N + m, k * N + m, vec_size, buf, atf)

    for k in range(j - 1):
        for m in range(k + 1, j):
            buf = _mul_(k * N + m, k * (N * N) + m * N + j, j, vec_size, buf, atf, col)

    for k in range(j - 2):
        for m in range(k + 1, j - 1):
            for i in range(m + 1, j):
                buf = _madd_(
                    k * N + m, k * (N * N) + m * N + i, i, vec_size, buf, atf, col
                )

    for k in range(j - 2):
        for m in range(k + 2, j):
            for i in range(k + 1, m):
                buf = _msub_(
                    k * N + m, k * (N * N) + i * N + m, i, vec_size, buf, atf, col
                )

    for k in range(1, j - 1):
        for m in range(k + 1, j):
            for i in range(k):
                buf = _madd_(
                    k * N + m, i * (N * N) + k * N + m, i, vec_size, buf, atf, col
                )

    for k in range(j - 1):
        for m in range(k + 1, j):
            for i in range(j + 1, N):
                buf = _madd_(
                    k * N + m, k * (N * N) + m * N + i, i, vec_size, buf, atf, col
                )

    for k in range(j - 1):
        for m in range(k + 1, j):
            atf = _mod_red_(k * (N * N) + m * N + j, k * N + m, vec_size, buf, atf)

    return atf


def inverting_on_atf(atf_in, columns):
    atf = decompress_atf(atf_in, C)

    diagonal = [0] * (N * C)
    column = [0] * (N * C)
    for i in range(N):
        for r in range(C):
            diagonal[i * C + r] = columns[i * (N + 1) * C + r]

    diagonal = set_inversion_mod_p(diagonal)

    for j in range(N - 1, -1, -1):
        for i in range(N):
            if i != j:
                for r in range(C):
                    column[i * C + r] = reduction_mod_p(
                        multiplication_mod_p(
                            (PRIME - diagonal[j * C + r]), columns[(j * N + i) * C + r]
                        )
                    )
        for r in range(C):
            column[j * C + r] = diagonal[j * C + r]
        atf = acting_on_atf(atf, column, j, C)

    return compress_atf(atf, C, C)


def acting_on_atfs(atf_in, columns, vec_size):
    atf = decompress_atf(atf_in, vec_size)

    for j in range(N):
        column = columns[j * N * vec_size : (j + 1) * N * vec_size]
        atf = acting_on_atf(atf, column, j, vec_size)

    return compress_atf(atf, vec_size, vec_size)

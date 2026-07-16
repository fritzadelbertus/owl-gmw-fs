from owl_mlip.ntrugen.fft import inv_fft
from owl_mlip.ntrugen.ntrugen_hawk import ntru_solve
from owl_mlip.params import LOGN, params
from owl_mlip.poly import adjoint, infnorm, isinvertible, l2norm, poly_add, poly_mul_ntt
from owl_mlip.typing import GroupElementT, SetElementT
from rngcontext import RngContext, shake256x4


def regeneratefg(kgseed: bytes, n: int) -> tuple[list[int], list[int]]:
    b = int(n / 64)
    assert b == 4 or b == 8 or b == 16

    y = shake256x4(kgseed, 2 * n * b // 64)

    ybits = [0] * b * 2 * n
    for j, z in enumerate(y):
        for bi in range(64):
            ybits[j * 64 + bi] = (z >> bi) & 1

    f = [0] * n
    for i in range(n):
        sum = 0
        for j in range(b):
            sum += ybits[i * b + j]
        f[i] = sum - b // 2

    g = [0] * n
    for i in range(n):
        sum = 0
        for j in range(b):
            sum += ybits[(i + n) * b + j]
        g[i] = sum - b // 2

    return (f, g)


def validfg(f: list[int], g: list[int], n: int) -> bool:
    invertibility = isinvertible(f, 2) and isinvertible(g, 2)
    l2norm_check = (l2norm(f) + l2norm(g)) > (2 * n * (params(LOGN, "sigmakrsec") ** 2))
    return invertibility and l2norm_check


def sample_basis_matrix(
    rng: RngContext,
) -> tuple[list[int], list[int], list[int], list[int]]:
    n = 1 << LOGN
    kgseed = rng.random(params(LOGN, "lenkgseed")).tobytes()
    f, g = regeneratefg(kgseed, n)

    while not validfg(f, g, n):
        kgseed = rng.random(params(LOGN, "lenkgseed")).tobytes()
        f, g = regeneratefg(kgseed, n)

    try:
        upper_f, upper_g = ntru_solve(f, g)
    except ValueError:
        return sample_basis_matrix(rng)

    if infnorm(upper_f) > 127 or infnorm(upper_g) > 127:
        return sample_basis_matrix(rng)

    return f, g, upper_f, upper_g


def sample_group(rng: RngContext) -> GroupElementT:
    f, g, upper_f, upper_g = sample_basis_matrix(rng)

    p = 8380417
    f = [i % p for i in f]
    upper_f = [i % p for i in upper_f]
    g = [i % p for i in g]
    upper_g = [i % p for i in upper_g]
    return [[f, upper_f], [g, upper_g]]


def sample_set(rng: RngContext) -> SetElementT:
    f, g, upper_f, upper_g = sample_basis_matrix(rng)

    f_adjoint = adjoint(f)
    g_adjoint = adjoint(g)

    p0 = (1 << 16) + 1
    q00 = poly_add(poly_mul_ntt(f, f_adjoint, p0), poly_mul_ntt(g, g_adjoint, p0))
    p1, p2 = (2147473409, 2147389441)
    if not isinvertible(q00, p1) or not isinvertible(q00, p2):
        return sample_set(rng)
    invq00 = inv_fft(q00)
    if invq00[0] >= params(LOGN, "beta0"):
        return sample_set(rng)

    upper_f_adjoint = adjoint(upper_f)
    upper_g_adjoint = adjoint(upper_g)

    q01 = poly_add(
        poly_mul_ntt(upper_f, f_adjoint, p0), poly_mul_ntt(upper_g, g_adjoint, p0)
    )
    q10 = adjoint(q01)
    p = 8380417
    q11 = poly_add(
        poly_mul_ntt(upper_f, upper_f_adjoint, p),
        poly_mul_ntt(upper_g, upper_g_adjoint, p),
    )

    if any(abs(q11i) >= 2 ** (params(LOGN, "high11")) for q11i in q11[1:]):
        return sample_set(rng)

    q00 = [i % p for i in q00]
    q01 = [i % p for i in q01]
    q10 = [i % p for i in q10]
    q11 = [i % p for i in q11]
    return [[q00, q01], [q10, q11]]

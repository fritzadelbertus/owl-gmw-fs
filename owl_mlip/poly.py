# This file contains several helper functions of polynomial arithmetic

import math

from sympy.ntheory import primitive_root


###############################
# NTT / invNTT
###############################
def brv(x: int, b: int = 32) -> int:
    r = 0
    for i in range(b):
        r |= ((x >> i) & 1) << (b - 1 - i)
    return r


def ntt(f: list[int], p: int, zetas: list[int]) -> list[int]:
    upper_f = [int(f) for f in f]
    deg = len(f)
    length = deg // 2
    k = 1
    while length > 0:
        for s in range(0, deg, 2 * length):
            zeta = zetas[k]
            k += 1
            for j in range(s, s + length):
                t = (upper_f[j + length] * zeta) % p
                upper_f[j + length] = (upper_f[j] - t) % p
                upper_f[j] = (upper_f[j] + t) % p
        length = length // 2
    return upper_f


def intt(f: list[int], p: int, zetas: list[int]) -> list[int]:
    # if zetas is None:
    #    _, zetas = get_roots(p, len(f))
    upper_f = [int(f) for f in f]
    deg = len(f)
    length = 1
    k = deg - 1
    while length < deg:
        for s in reversed(range(0, deg, 2 * length)):
            zeta = zetas[k]
            k -= 1
            for j in range(s, s + length):
                t = upper_f[j]
                upper_f[j] = (t + upper_f[j + length]) % p
                upper_f[j + length] = (t - upper_f[j + length]) % p
                upper_f[j + length] = (upper_f[j + length] * zeta) % p
        length = length * 2

    ideg = pow(deg, p - 2, p)
    return [(f * ideg) % p for f in upper_f]


def nttadj(u: list[int], p: int) -> list[int]:
    zetas, izetas = get_roots(p, len(u))
    ui = intt(u, p, izetas)
    return ntt(adjoint(ui), p, zetas)


def get_roots(p: int, n: int) -> tuple[list[int], list[int]]:
    logn = int(math.log2(n))
    g0 = primitive_root(p)
    b = (p - 1) // (2 * n)
    g0 = (g0**b) % p

    zetas = [pow(g0, brv(u, logn), p) for u in range(1 << logn)]
    izetas = [pow(z, p - 2, p) for z in zetas]
    return zetas, izetas


########################################
# Polynomial operations
########################################
def poly_mul_ntt(p1: list[int], p2: list[int], p: int) -> list[int]:
    n = len(p1)
    zetas, izetas = get_roots(p, n)
    p1ntt = ntt(p1, p, zetas)
    p2ntt = ntt(p2, p, zetas)
    rntt = [(a * b) % p for a, b in zip(p1ntt, p2ntt, strict=True)]
    m = intt(rntt, p, izetas)
    for i in range(n):
        if m[i] > (p - 1) // 2:
            m[i] -= p
    return m


def poly_mul_schoolbook(p1, p2, p):
    deg = len(p1)
    r = [0 for _ in range(deg)]
    for i in range(deg):
        for j in range(deg):
            ij = i + j
            if ij >= deg:
                r[ij % deg] -= (p1[i] * p2[j]) % p
            else:
                r[ij % deg] += (p1[i] * p2[j]) % p
    return [x % p for x in r]


def adjoint(u: list[int]) -> list[int]:
    ustar = u.copy()
    n = len(u)
    for i in range(1, n):
        ustar[i] = -u[n - i]
    return ustar


def poly_sub(p0: list[int], p1: list[int]) -> list[int]:
    return [p0 - p1 for (p0, p1) in zip(p0, p1, strict=True)]


def poly_add(p0: list[int], p1: list[int]) -> list[int]:
    return [p0 + p1 for (p0, p1) in zip(p0, p1, strict=True)]


########################################
# Polynomial properties
########################################
def infnorm(poly: list[int]) -> int:
    return max([abs(p) for p in poly])


def bytes_to_poly(h, n):
    h0 = [None] * n
    for i in range(n):
        h0[i] = (h[i // 8] >> (i % 8)) & 1
    return h0


def l2norm(x: list[int]) -> int:
    return sum([a * a for a in x])


def isinvertible(poly: list[int], p: int) -> bool:
    if p == 2:
        return (sum(poly) % 2) == 1
    _, zetas = get_roots(p, len(poly))
    polyntt = ntt(poly, p, zetas)
    return all(c != 0 for c in polyntt)

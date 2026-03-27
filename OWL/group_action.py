from OWL.hawk.poly import adjoint, isinvertible, poly_mul_ntt, poly_add, poly_sub, infnorm, l2norm
from OWL.hawk.rngcontext import RngContext, SHAKE256x4
from OWL.hawk.params import PARAMS
import OWL.hawk.ntrugen
from OWL.hawk.ntrugen.ntrugen_hawk import ntru_solve

import numpy as np

from OWL.params import LOGN as logn
# ===================================================
# NTT Operations

def poly_add_p(p0,p1,p):
    return [(p0 + p1)%p for (p0, p1) in zip(p0, p1)]

def matrix_mult(A,B):
    p = 8380417
    c00 = poly_add_p(poly_mul_ntt(A[0][0], B[0][0],p), poly_mul_ntt(A[0][1], B[1][0],p),p)
    c01 = poly_add_p(poly_mul_ntt(A[0][0], B[0][1],p), poly_mul_ntt(A[0][1], B[1][1],p),p)
    c10 = poly_add_p(poly_mul_ntt(A[1][0], B[0][0],p), poly_mul_ntt(A[1][1], B[1][0],p),p)    
    c11 = poly_add_p(poly_mul_ntt(A[1][0], B[0][1],p), poly_mul_ntt(A[1][1], B[1][1],p),p)
    return [[c00,c01],[c10,c11]]

# ==============================================================================================================
# The Group Action

def action(B, Q):
    B_star = [[adjoint(B[0][0]),adjoint(B[1][0])],[adjoint(B[0][1]),adjoint(B[1][1])]]
    return matrix_mult(matrix_mult(B_star, Q), B)


# ==============================================================================================================
# Group and Set Sampling

# Regenerate f,g from basis matrix (from HAWK)
def regeneratefg(kgseed, n):
    b = int(n / 64)
    assert b == 4 or b == 8 or b == 16

    y = SHAKE256x4(kgseed, 2 * n * b // 64)  # array of 64-bit values

    # map y to a sequence of bits
    ybits = [None] * b * 2 * n
    for j, y in enumerate(y):
        for bi in range(64):
            ybits[j * 64 + bi] = (y >> bi) & 1

    f = [int(0)] * n
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

# Sample a random basis matrix
def group_sampler(rng=None): # Samples a random matrix
    n = 1 << logn

    if rng is None:
        rng = RngContext(np.random.randint(0, 256, 40, dtype=np.uint8))

    # Line 1: kgseen <- Rnd(len_bits(kgseed))
    kgseed = rng.random(PARAMS(logn, "lenkgseed"))

    # Line 2: (f, g) <- Regeneratefg(kgseed)
    f, g = regeneratefg(kgseed.tobytes(), n)

    # Line 3: if isInvertible(f, 2) false or isInvertible(f, 2) is false then
    if not isinvertible(f, 2) or not isinvertible(g, 2):
        # Line 4: restart
        return group_sampler(rng)

    # Line 5: if ||(f,g)||2 <= 2*n*sigkrsec**2:
    if (l2norm(f) + l2norm(g)) <= (2 * n * (PARAMS(logn, "sigmakrsec") ** 2)):
        return group_sampler(rng)

    try:
        # Line 13: r <- NTRUSolve(f, g, 1)
        # Line 16: (F, G) <- r
        F, G = ntru_solve(f, g)
    except ValueError:
        # Line 14&15: if r = \bot then restart
        return group_sampler(rng)

    # Line 17: if infnorm( (F,G) ) > 127 then
    if infnorm(F) > 127 or infnorm(G) > 127:
        # Line 18: restart
        return group_sampler(rng)

    # Line 28: return (priv, pub)
    p = 8380417
    f = [i%p for i in f]
    F = [i%p for i in F]
    g = [i%p for i in g]
    G = [i%p for i in G]
    return [[f,F],[g,G]]

# Sample a random gram matrix
def set_sampler(rng=None):
    n = 1 << logn

    if rng is None:
        rng = RngContext(np.random.randint(0, 256, 40, dtype=np.uint8))

    # Line 1: kgseen <- Rnd(len_bits(kgseed))
    kgseed = rng.random(PARAMS(logn, "lenkgseed"))

    # Line 2: (f, g) <- Regeneratefg(kgseed)
    f, g = regeneratefg(kgseed.tobytes(), n)

    # Line 3: if isInvertible(f, 2) false or isInvertible(f, 2) is false then
    if not isinvertible(f, 2) or not isinvertible(g, 2):
        # Line 4: restart
        return set_sampler(rng)

    fadj = adjoint(f)
    gadj = adjoint(g)

    # Line 5: if ||(f,g)||2 <= 2*n*sigkrsec**2:
    if (l2norm(f) + l2norm(g)) <= (2 * n * (PARAMS(logn, "sigmakrsec") ** 2)):
        return set_sampler(rng)

    p = (1 << 16) + 1

    # Line 7: q00 <- f f* + g g*
    q00 = poly_add(poly_mul_ntt(f, fadj, p), poly_mul_ntt(g, gadj, p))

    # Line 8: (p1, p2) <- (2147473409, 2147389441)
    p1, p2 = (2147473409, 2147389441)

    # Line 9: if isInvertible(q00, p1) false or isInvertible(q00, p2) is false then
    if not isinvertible(q00, p1) or not isinvertible(q00, p2):
        # Line 10: restart
        return set_sampler(rng)

    # Line 11: if (1/q00)[0] >= beta0 then
    invq00 = OWL.hawk.ntrugen.fft.inv_fft(q00)
    if invq00[0] >= PARAMS(logn, "beta0"):
        # Line 12: restart
        return set_sampler(rng)

    try:
        # Line 13: r <- NTRUSolve(f, g, 1)
        # Line 16: (F, G) <- r
        F, G = ntru_solve(f, g)
    except ValueError:
        # Line 14&15: if r = \bot then restart
        return set_sampler(rng)

    # Line 17: if infnorm( (F,G) ) > 127 then
    if infnorm(F) > 127 or infnorm(G) > 127:
        # Line 18: restart
        return set_sampler(rng)

    Fadj = adjoint(F)
    Gadj = adjoint(G)

    # Line 19: q01 <- F f* + G g*
    q01 = poly_add(poly_mul_ntt(F, fadj, p), poly_mul_ntt(G, gadj, p))
    q10 = adjoint(q01)
    p = 8380417
    # Line 20: q11 <- F F* + G G*
    q11 = poly_add(poly_mul_ntt(F, Fadj, p), poly_mul_ntt(G, Gadj, p))

    # Line 21: if |q11[i]| >= 2^high11 for any i > 0 then
    if any(abs(q11i) >= 2 ** (PARAMS(logn, "high11")) for q11i in q11[1:]):
        # Line 22: restart
        return set_sampler(rng)
    q00 = [i%p for i in q00]
    q01 = [i%p for i in q01]
    q10 = [i%p for i in q10]
    q11 = [i%p for i in q11]
    return [[q00,q01],[q10,q11]]

# ==============================================================================================================
# Group Operations

# Calculate group inverse
def group_inverse(A):
    p = 8380417
    a = A[0][0]
    minus_b = [(i*(-1))%p for i in A[0][1]]
    minus_c = [(i*(-1))%p for i in A[1][0]]
    d = A[1][1]
    return [[d, minus_b],[minus_c,a]]

# Calculate group operation (matrix multiplication)
def group_operator(A,B):
    return matrix_mult(B,A)


poly_zeros = [0 for i in range(2**logn)]
poly_unit = [0 for i in range(2**logn)]
poly_unit[0] = 1
group_identity = [[poly_unit, poly_zeros],[poly_zeros, poly_unit]]
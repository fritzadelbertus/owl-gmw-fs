from HAWK.poly import adjoint, isinvertible, poly_mul_ntt, poly_add, infnorm, l2norm, poly_sub, bytes_to_poly, brv, get_roots, ntt, nttadj
from HAWK.ntrugen.ntrugen_hawk import ntru_solve
import HAWK.ntrugen
from HAWK.params import PARAMS
from HAWK.codec import encode_public, encode_private, decode_private, encode_sign, decode_public, decode_sign
import hashlib
import numpy as np
from HAWK.rngcontext import RngContext, SHAKE256x4

from sympy.ntheory import npartitions

def hawkkeygen(logn, rng=None):
    """
    hawkkeygen (see Alg 13)

    Regenerate a priv/pub key pair for Hawk

    Inputs:
        - logn : log2 of polynomial degree
        - rng : rng context used

    Outputs:
        - priv : encoded private key (uint8)
        - pub : encoded public key (uint8)
    """
    if rng is None:
        rng = RngContext(np.random.randint(0, 256, 40, dtype=np.uint8))

    _, _, _, _, _, _, _, priv, pub = hawkkeygen_unpacked(logn, rng)
    return priv, pub


def hawkkeygen_unpacked(logn, rng=None):
    """
    hawkkeygen_unpacked (see Alg 13)

    This is used as a subroutine of hawkkeygen used to exposed directly polynomials

    Inputs:
        - logn : log2 of polynomial degree
        - rng : rng context used

    Outputs:
        - f : polynomial (int8)
        - g : polynomial (int8)
        - F : polynomial (int8)
        - G : polynomial (int8)
        - q00: polynomial (int16)
        - q01: polynomial (int16)
        - kgseed: seed to regenerate (f,g) (uint8)
        - priv: encoded private key
        - pub: encoded public key
    """
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
        return hawkkeygen_unpacked(logn, rng)

    fadj = adjoint(f)
    gadj = adjoint(g)

    # Line 5: if ||(f,g)||2 <= 2*n*sigkrsec**2:
    if (l2norm(f) + l2norm(g)) <= (2 * n * (PARAMS(logn, "sigmakrsec") ** 2)):
        return hawkkeygen_unpacked(logn, rng)

    p = (1 << 16) + 1

    # Line 7: q00 <- f f* + g g*
    q00 = poly_add(poly_mul_ntt(f, fadj, p), poly_mul_ntt(g, gadj, p))

    # Line 8: (p1, p2) <- (2147473409, 2147389441)
    p1, p2 = (2147473409, 2147389441)

    # Line 9: if isInvertible(q00, p1) false or isInvertible(q00, p2) is false then
    if not isinvertible(q00, p1) or not isinvertible(q00, p2):
        # Line 10: restart
        return hawkkeygen_unpacked(logn, rng)

    # Line 11: if (1/q00)[0] >= beta0 then
    invq00 = HAWK.ntrugen.fft.inv_fft(q00)
    if invq00[0] >= PARAMS(logn, "beta0"):
        # Line 12: restart
        return hawkkeygen_unpacked(logn, rng)

    try:
        # Line 13: r <- NTRUSolve(f, g, 1)
        # Line 16: (F, G) <- r
        F, G = ntru_solve(f, g)
    except ValueError:
        # Line 14&15: if r = \bot then restart
        return hawkkeygen_unpacked(logn, rng)

    # Line 17: if infnorm( (F,G) ) > 127 then
    if infnorm(F) > 127 or infnorm(G) > 127:
        # Line 18: restart
        return hawkkeygen_unpacked(logn, rng)

    Fadj = adjoint(F)
    Gadj = adjoint(G)

    # Line 19: q01 <- F f* + G g*
    q01 = poly_add(poly_mul_ntt(F, fadj, p), poly_mul_ntt(G, gadj, p))
    p = 8380417
    # Line 20: q11 <- F F* + G G*
    q11 = poly_add(poly_mul_ntt(F, Fadj, p), poly_mul_ntt(G, Gadj, p))

    # Line 21: if |q11[i]| >= 2^high11 for any i > 0 then
    if any(abs(q11i) >= 2 ** (PARAMS(logn, "high11")) for q11i in q11[1:]):
        # Line 22: restart
        return hawkkeygen_unpacked(logn, rng)

    # Line 23: pub <- EncodePublic(q00, q01)
    pub = encode_public(logn, q00, q01)

    # Line 24: if pub = \bot then
    if pub is None:
        # Line 25: restart
        return hawkkeygen_unpacked(logn, rng)

    # Line 16: hpub <- SHAKE256(pub)
    shake256 = hashlib.shake_256()
    shake256.update(pub.tobytes())
    hpub = np.frombuffer(shake256.digest(PARAMS(logn, "lenhpub")), dtype=np.uint8)

    # Line 27: priv <- EncodePrivate(kgseen, F mod 2, G mod 2, hpub)
    priv = encode_private(kgseed, F, G, hpub)

    # Line 28: return (priv, pub)
    return f, g, F, G, q00, q01, kgseed, priv, pub


def regeneratefg(kgseed, n):
    """
    regeneratefg (see Alg 12)

    Regenerates (f,g) deterministically from a seed

    Inputs:
        - kgseed : seed
        - n : polynomial degree of f and g

    Outputs:
        - f : polynomial of int
        - g : polynomial of int
    """
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




def samplersign(s, t, n, T0, T1):
    """
    SampleSign (see Alg 14)

    Inputs:
        - s randomness seed (bytes)
        - t center vector
        - n degree of polynomials
        - T0 / T1 CDT tables

    Outputs:
        - d: discret Gaussian polynomial
    """
    y = SHAKE256x4(s, int(5 * n / 2))
    d = [None] * 2 * n

    for j in range(4):
        for i in range(int(n / 8)):
            for k in range(3 + 1):
                r = 16 * i + 4 * j + k
                a = y[j + 4 * (5 * i + k)]
                b = (y[j + 4 * (5 * i + 4)] // (2 ** (16 * k))) % 2**15
                c = (a % 2**63) + 2**63 * b
                (v0, v1) = (0, 0)
                z = 0
                while T0[z] != 0 and T1[z] != 0:
                    if c < T0[z]:
                        v0 = v0 + 1
                    if c < T1[z]:
                        v1 = v1 + 1
                    z = z + 1
                if t[r] == 0:
                    v = 2 * v0
                else:
                    v = 2 * v1 + 1
                if a >= 2**63:
                    v = -v
                d[r] = v
    return d


def symbreak(w):
    """
    symbreak (see Section 3.5.2)

    Inputs:
        - w input polynomial

    Outputs:
        - True iff w[i] > 0 and w[j] = 0 for all j < i
    """
    for x in w:
        if x != 0:
            if x > 0:
                return True
            else:
                return False
    return False


def hawksign(logn, priv, msg, rng=None):
    """
    hawksign (see Alg. 15)

    Inputs:
        - logn log2 of polynomial degree
        - priv private key
        - msg  mesasge to sign
        - rng RngContext used during signing

    outputs:
        - valid signature of msg under the private key
    """
    # Line 1: (kgseed, F mod 2, G mod 2, hpub) <- DecodePrivate(priv)
    kgseed, Fmod2, Gmod2, hpub = decode_private(logn, priv)

    while True:
        salt, s1 = hawksign_unpacked(logn, Fmod2, Gmod2, kgseed, hpub, msg, rng)

        # Line 18: EncodeSignature(salt, s1)
        sig = encode_sign(logn, salt, s1)

        # Line 19: if sig != bottom then
        if sig is not None:
            # Line 20: return sig
            return sig


def hawksign_unpacked(logn, Fmod2, Gmod2, kgseed, hpub, msg, rng=None):
    """
    hawksign_unpacked (see Alg. 15)

    This is the subroutine of hawksign with decoded / unpacked private key.

    Inputs:
        - logn log2 of polynomial degree
        - Fmod2
        - Gmod2
        - kgseed seed to reconstruct (f,g)
        - hpub hash of public key
        - msg message to sign
        - rng RngContext used during signing

    Outputs:
        - salt to generate (h0,h1)
        - s1 signature polynomial
    """
    if rng is None:
        rng = RngContext(np.random.randint(0, 256, 32))

    n = 1 << logn

    # Line 2: (f,g) <- Regenerate(kgseed)
    f, g = regeneratefg(kgseed.tobytes(), n)

    # Line 3: M <- SHAKE256(m || hpub)[0:512]
    shake256 = hashlib.shake_256()
    shake256.update(msg.tobytes() + hpub.tobytes())
    M = shake256.digest(64)

    # Line 4: a <- 0
    a = 0

    # prime used for NTTs
    p = (1 << 16) + 1

    # Line 5: loop
    while True:
        # Line 6: salt <- SHAKE256(M||kgseed||EncodeInt(a,32)||Rnd(saltlen))[0:saltlenbits]
        shake256 = hashlib.shake_256()
        shake256.update(
            M
            + kgseed.tobytes()
            + a.to_bytes(4, "little")
            + rng.random(PARAMS(logn, "lensalt")).tobytes()
        )
        salt = np.frombuffer(shake256.digest(PARAMS(logn, "lensalt")), dtype=np.uint8)

        # Line 7: (h0,h1) <- SHAKE256(M||salt)[0:2n]
        shake256 = hashlib.shake_256()
        shake256.update(M + salt.tobytes())

        h = shake256.digest(n // 4)
        h0 = bytes_to_poly(h[: n // 8], n)
        h1 = bytes_to_poly(h[n // 8 :], n)

        # Line 8: (t0,t1) <- ((h0*f + f1*F) mod 2, (h0*g + f1*G) mod 2)
        t0 = [
            x % 2 for x in poly_add(poly_mul_ntt(h0, f, p), poly_mul_ntt(h1, Fmod2, p))
        ]
        t1 = [
            x % 2 for x in poly_add(poly_mul_ntt(h0, g, p), poly_mul_ntt(h1, Gmod2, p))
        ]

        # Line 9: s <- M || kgseed || EncodeInt(a+1,32) || Rnd(320)
        seed = rng.random(320 // 8)
        s = M + kgseed.tobytes() + (a + 1).to_bytes(4, "little") + seed.tobytes()

        # Line 10: (d0, d1) <- SampleSign(s,(t0,t1))
        d = samplersign(s, t0 + t1, n, PARAMS(logn, "T0"), PARAMS(logn, "T1"))
        d0 = d[:n]
        d1 = d[n:]

        # Line 11: a <- a + 2
        a = a + 2

        # Line 12: if ||(d0,d1)||2 > 8 * n * sigm_verify^2
        if l2norm(d0) + l2norm(d1) > 8 * n * (PARAMS(logn, "sigmaverify") ** 2):
            # Line 13: continue loop
            continue

        # Line 14: f * d1 - g * d0
        w1 = poly_sub(poly_mul_ntt(f, d1, p), poly_mul_ntt(g, d0, p))

        # Line 15: if sym-break(w1) ==
        if not symbreak(w1):
            # Line 16: w1 <- -w1
            w1 = [-x for x in w1]

        # Line 17: s <- (h1 - w1) / 2
        sig = [x // 2 for x in poly_sub(h1, w1)]

        return salt, sig


def hawkverify(logn, pub, msg, sig):
    """
    hawkverify (see Alg 20)

    Performs signature verification for Hawk

    Inputs:
        - logn : log2 of polynomial degree
        - pub : public key
        - msg : message
        - sig : candidate signature

    Outputs:
        - True if sig is a valid signature of msg, False otherwise
    """
    n = 1 << logn
	
	# Line 1: r <- DecodeSignature(sig)
    r = decode_sign(logn, sig)
	# Line 2: if r == ⊥ then
    if r is None:
		#Line 3: return False
        return False
	# Line 4: (salt,s1) <- r 
    salt, s1 = r

	# Line 5: DecodePublic(pub)
    r = decode_public(logn, pub)
	# Line 6: if r == ⊥ then
    if r is None:
		#Line 7: return False
        return False
	# Line 8: (q00,q01) <- r
    q00, q01 = r
	
	# Line 9: hpub <- SHAKE256(pub)
    shake256 = hashlib.shake_256()
    shake256.update(pub.tobytes())
    hpub = shake256.digest(PARAMS(logn, "lenhpub"))

	# Line 10: M <- SHAKE256( m || hpub)
    shake256 = hashlib.shake_256()
    shake256.update(msg.tobytes() + hpub)
    M = shake256.digest(64)

	# Line 11: (h0,h1) <- SHAKE256(M||salt)[0:2n]
    shake256 = hashlib.shake_256()
    shake256.update(M + salt.tobytes())
    h = shake256.digest(n // 4)
    h0 = bytes_to_poly(h[: n // 8], n)
    h1 = bytes_to_poly(h[n // 8 :], n)

    f = hawkverify_unpacked(logn, s1, q00, q01, h0, h1)
    return f


def hawkverify_unpacked(logn, s1, q00, q01, h0, h1):
    """
    hawkverify_unpacked (see Alg 20)

    This is used a subroutine for hawkverify, operating on polynomials and not encoded sig / pub
    """
    n = 1 << logn

	# Line 12: w1 <- h1 - 2 * s1
    w1 = poly_sub(h1, 2 * np.array(s1))

	# Line 13: if sym-break(w1) == false then 
    if symbreak(w1) == False:
		# Line 14: Return False 
        return False
		
	# Line 15: w0 <- RebuildS0(q00,q01,w1,h0)
    w0 = rebuilds0(logn, q00, q01, w1, h0)
	# Line 16: if w0 == ⊥ then
    if w0 is None:
		# Line 17: return False 
        return False

    p1, p2 = (2147473409, 2147389441)
	
	# Line 18: r1 <- PolyQnorm(q00,q01,w0,w1,p1)
    r1 = polyQnorm(q00, q01, w0, w1, p1)
	# Line 19: r2 <- PolyQnorm(q00,q01,w0,w1,p2)
    r2 = polyQnorm(q00, q01, w0, w1, p2)

	# Line 20: if r1 != r2 or r1 != 0 mod n then
    if r1 != r2 or r1 % n != 0:
		# Line 21: return False 
        return False
	
	# Line 22 : r1 <- r1 / n
    r1 = r1 // n

	# Line 23: if r1 > 8 * n * sigmaverfy**2 then 
    if r1 > 8 * n * PARAMS(logn, "sigmaverify") ** 2:
		# Line 24: return False 
        return False

	# Line 25: return True
    return True


def rebuilds0(logn, q00, q01, w1, h0):
    """
    rebuilds0 (see Alg 18)

    Rebuild s0 from public key and signature polynomials

    Inputs:
        - logn : log2 of polynomial degree
        - q00 : public polynomial (pub key)
        - q01 : public polynomial (pub key)
        - w1 : public polynomial (sig)
        - h0 : public polynomial (h0)

    Outputs:
        - w0: recovered s0
    """

    n = 1 << logn
    cw1 = 2 ** (29 - (1 + PARAMS(logn, "highs1")))
    cq00 = 2 ** (29 - PARAMS(logn, "high00"))
    cq01 = 2 ** (29 - PARAMS(logn, "high01"))
    cs0 = (2 * cw1 * cq01) / (n * cq00)
    w1scalled = (np.array(w1, dtype=np.int32)) * cw1
    w1fft = fft(w1scalled)
    z00 = q00.copy()
    if z00[0] < 0:
        return False
    z00[0] = 0
    q00fft = fft(np.array(z00, dtype=np.int32) * cq00)
    q01fft = fft(np.array(q01, dtype=np.int32) * cq01)
    alpha = (2 * cq00 * np.int64(q00[0])) // n
    for u in range(0, n // 2):
        x_re = np.int64(q01fft[u]) * np.int64(w1fft[u])
        x_re -= np.int64(q01fft[u + n // 2]) * np.int64(w1fft[u + n // 2])

        x_im = np.int64(q01fft[u]) * np.int64(w1fft[u + n // 2])
        x_im += np.int64(q01fft[u + n // 2]) * np.int64(w1fft[u])

        (x_re, z_re) = (np.abs(x_re), sgn(x_re))
        (x_im, z_im) = (np.abs(x_im), sgn(x_im))

        v = alpha + q00fft[u]
        if v <= 0 or v >= 2**32 or x_re >= v * 2**32 or x_im >= v * 2**32:
            return False

        y_re = np.int32(x_re // v)
        y_im = np.int32(x_im // v)

        q01fft[u] = y_re - 2 * z_re * y_re
        q01fft[u + n // 2] = y_im - 2 * z_im * y_im

    t = invfft(q01fft)
    w0 = np.zeros(n, dtype=np.int32)
    for u in range(n):
        v = cs0 * np.int32(h0[u]) + np.int32(t[u])
        z = (v + cs0) // (2 * cs0)
        if z < -(2 ** PARAMS(logn, "highs0")) or z >= 2 ** (PARAMS(logn, "highs0")):
            return False
        w0[u] = np.int32(h0[u]) - 2 * z

    return w0


def polyQnorm(q00, q01, w0, w1, p):
    """
    polyQnorm (see Alg 19)

    Computes the Q norm of w

    Inputs:
        - q00 : polynomial
        - q01 : polynomial
        - w0 : polynomial
        - w1 : polynomial
        - p : prime

    Outputs:
        - n ||W||_Q**2 mod p
    """
    zetas, _ = get_roots(p, len(q00))
    q00ntt = ntt(q00, p, zetas)
    q01ntt = ntt(q01, p, zetas)
    w0ntt = ntt(w0, p, zetas)
    w1ntt = ntt(w1, p, zetas)

    d = [(w1 * pow(int(q00), p - 2, p)) % p for (w1, q00) in zip(w1ntt, q00ntt)]
    c = [(d * w1adj) % p for (d, w1adj) in zip(d, ntt(adjoint(w1), p))]
    acc = sum(c)
    e = [(int(w0) + d * int(q01)) % p for (w0, d, q01) in zip(w0ntt, d, q01ntt)]
    c = [
        (int(q00) * int(e) * int(eadj)) % p
        for (q00, e, eadj) in zip(q00ntt, e, nttadj(e, p))
    ]
    acc += sum(c)
    return acc % p


###########################
## FFT and invFFT fixed point
###########################
def delta(k):
    """
    delta (see Section 3.6)

    Computes roots

    Inputs:
        - k : index

    Outputs:
        - re : real part
        - im : imaginary part
    """
    d = np.exp(2 * 1j * np.pi * brv(k, b=10) / 2048)
    d = d * (2**31)
    re = np.int32(np.round(d.real))
    im = np.int32(np.round(d.imag))
    return re, im


def sgn(x):
    """
    sgn

    Returns the sign bit of x

    Inputs:
        - x : integer

    Outputs:
        - 1 if x < 0, 0 otherwise
    """
    if x < 0:
        return 1
    return 0


def fft(a):
    """
    fft (see Alg 16)

    Computes the fixed point fft

    Inputs:
        - a : polynomial

    Outputs:
        - afft : fft representation of a

    """
    n = len(a)
    afft = a.copy()
    t = n // 2
    m = 2
    while m < n:
        v0 = 0
        for u in range(0, m // 2):
            e_re, e_im = delta(u + m)
            e_re = np.int64(e_re)
            e_im = np.int64(e_im)

            for v in range(v0, v0 + t // 2):
                x1_re = np.int64(afft[v])
                x1_im = np.int64(afft[v + n // 2])

                x2_re = np.int64(afft[v + t // 2])
                x2_im = np.int64(afft[v + t // 2 + n // 2])

                t_re = x2_re * e_re - x2_im * e_im
                t_im = x2_re * e_im + x2_im * e_re

                afft[v] = np.int32((2**31 * x1_re + t_re) // 2**32)
                afft[v + n // 2] = np.int32((2**31 * x1_im + t_im) // 2**32)
                afft[v + t // 2] = np.int32((2**31 * x1_re - t_re) // 2**32)
                afft[v + t // 2 + n // 2] = np.int32(
                    (2**31 * x1_im - t_im) // 2**32
                )
            v0 = v0 + t
        t = t // 2
        m = 2 * m

    return afft


def invfft(afft):
    """
    invfft (see Alg 17)

    Computes the fixed point invfft

    Inputs:
        - afft : fft representaiton of a

    Outputs:
        - a : polynomial a
    """
    n = len(afft)
    a = afft.copy()
    t = 2
    m = n // 2
    while m > 1:
        v0 = 0
        for u in range(0, m // 2):
            eta_re, eta_im = delta(u + m)

            eta_re = np.int64(eta_re)
            eta_im = np.int64(-eta_im)

            for v in range(v0, v0 + t // 2):
                x1_re = a[v]
                x1_im = a[v + n // 2]

                x2_re = a[v + t // 2]
                x2_im = a[v + t // 2 + n // 2]

                t1_re = x1_re + x2_re
                t1_im = x1_im + x2_im
                t2_re = x1_re - x2_re
                t2_im = x1_im - x2_im

                a[v] = t1_re // 2
                a[v + n // 2] = t1_im // 2
                a[v + t // 2] = np.int32(
                    (np.int64(t2_re) * eta_re - np.int64(t2_im) * eta_im) // 2**32
                )
                a[v + t // 2 + n // 2] = np.int32(
                    (np.int64(t2_re) * eta_im + np.int64(t2_im) * eta_re) // 2**32
                )
            v0 = v0 + t
        t = 2 * t
        m = m // 2

    return a

from OWL.hawk.poly import adjoint, isinvertible, poly_mul_ntt, poly_add, poly_sub, infnorm, l2norm
from OWL.hawk.ntrugen.ntrugen_hawk import ntru_solve
import OWL.hawk.ntrugen
from OWL.hawk.params import PARAMS
import numpy as np
from OWL.hawk.rngcontext import RngContext
import hashlib

from OWL.coder import encode_poly_matrix, decode_poly_matrix
from OWL.group_action import regeneratefg, action, group_inverse, group_operator

from OWL.params import LOGN as logn
from OWL.params import ONE

# ===================================================================================================

CHLG_SIZE = ONE["LAMBDA"]//4

def pair_sampler_sign(rng=None):
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
        return pair_sampler_sign(rng)

    fadj = adjoint(f)
    gadj = adjoint(g)

    # Line 5: if ||(f,g)||2 <= 2*n*sigkrsec**2:
    if (l2norm(f) + l2norm(g)) <= (2 * n * (PARAMS(logn, "sigmakrsec") ** 2)):
        return pair_sampler_sign(rng)

    p = (1 << 16) + 1

    # Line 7: q00 <- f f* + g g*
    q00 = poly_add(poly_mul_ntt(f, fadj, p), poly_mul_ntt(g, gadj, p))

    # Line 8: (p1, p2) <- (2147473409, 2147389441)
    p1, p2 = (2147473409, 2147389441)

    # Line 9: if isInvertible(q00, p1) false or isInvertible(q00, p2) is false then
    if not isinvertible(q00, p1) or not isinvertible(q00, p2):
        # Line 10: restart
        return pair_sampler_sign(rng)

    # Line 11: if (1/q00)[0] >= beta0 then
    invq00 = OWL.hawk.ntrugen.fft.inv_fft(q00)
    if invq00[0] >= PARAMS(logn, "beta0"):
        # Line 12: restart
        return pair_sampler_sign(rng)

    try:
        # Line 13: r <- NTRUSolve(f, g, 1)
        # Line 16: (F, G) <- r
        F, G = ntru_solve(f, g)
    except ValueError:
        # Line 14&15: if r = \bot then restart
        return pair_sampler_sign(rng)

    # Line 17: if infnorm( (F,G) ) > 127 then
    if infnorm(F) > 127 or infnorm(G) > 127:
        # Line 18: restart
        return pair_sampler_sign(rng)

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
        return pair_sampler_sign(rng)


    # Line 28: return (priv, pub)
    return q00, q01, q11, f,g,F,G

def owl_Gen(rng=None):            
    #print("Generating Key...")
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
        return owl_Gen(rng)

    fadj = adjoint(f)
    gadj = adjoint(g)

    # Line 5: if ||(f,g)||2 <= 2*n*sigkrsec**2:
    if (l2norm(f) + l2norm(g)) <= (2 * n * (PARAMS(logn, "sigmakrsec") ** 2)):
        return owl_Gen(rng)

    p = (1 << 16) + 1

    # Line 7: q00 <- f f* + g g*
    q00 = poly_add(poly_mul_ntt(f, fadj, p), poly_mul_ntt(g, gadj, p))

    # Line 8: (p1, p2) <- (2147473409, 2147389441)
    p1, p2 = (2147473409, 2147389441)

    # Line 9: if isInvertible(q00, p1) false or isInvertible(q00, p2) is false then
    if not isinvertible(q00, p1) or not isinvertible(q00, p2):
        # Line 10: restart
        return owl_Gen(rng)

    # Line 11: if (1/q00)[0] >= beta0 then
    invq00 = OWL.hawk.ntrugen.fft.inv_fft(q00)
    if invq00[0] >= PARAMS(logn, "beta0"):
        # Line 12: restart
        return owl_Gen(rng)

    try:
        # Line 13: r <- NTRUSolve(f, g, 1)
        # Line 16: (F, G) <- r
        F, G = ntru_solve(f, g)
    except ValueError:
        # Line 14&15: if r = \bot then restart
        return owl_Gen(rng)

    # Line 17: if infnorm( (F,G) ) > 127 then
    if infnorm(F) > 127 or infnorm(G) > 127:
        # Line 18: restart
        return owl_Gen(rng)

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
        return owl_Gen(rng)

    #print("Key Generated!")

    # Line 28: return (priv, pub)
    return encode_poly_matrix([[[i%p for i in q00],[i%p for i in q01]],[[i%p for i in q10],[i%p for i in q11]]]), kgseed


def owl_Sign(privateKeyInBits, publicKeyInBits, messageInByte):
    #print("Signing Message...")
    n = 1 << logn
    f, g = regeneratefg(privateKeyInBits.tobytes(), n)
    F, G = ntru_solve(f, g)
    p = 8380417
    private_key = [[[i%p for i in f],[i%p for i in F]],[[i%p for i in g],[i%p for i in G]]]

    q00, q01, q11, hf,hg,hF,hG = pair_sampler_sign()
    q10 = adjoint(q01)
    t_i = [[[i%p for i in q00],[i%p for i in q01]],[[i%p for i in q10],[i%p for i in q11]]]

    h_i = [[[i%p for i in hf],[i%p for i in hF]],[[i%p for i in hg],[i%p for i in hG]]]

    hash_input = messageInByte + encode_poly_matrix(t_i)
    
    cha = hashlib.shake_256(hash_input).digest(CHLG_SIZE)

    f_i = group_operator(h_i, group_inverse(private_key))
    sign = cha + encode_poly_matrix(f_i)

    #print("Message Signed!")
    return sign

def owl_Vrfy(publicKey, messageInByte, sign):
    #print("Verifying Message...")

    public_key = decode_poly_matrix(publicKey, logn)
    # Get the b_i as integer
    cha = sign[:CHLG_SIZE]
    right_part = sign[CHLG_SIZE:]
    f_i = decode_poly_matrix(right_part, logn)

    t_i = action(f_i, public_key)
    hash_input = messageInByte + encode_poly_matrix(t_i)
    cha_2 = hashlib.shake_256(hash_input).digest(CHLG_SIZE)
    if cha == cha_2:
        return 0
    else:
        return 1

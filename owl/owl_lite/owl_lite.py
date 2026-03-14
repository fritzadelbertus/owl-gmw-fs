
from hawk_keygen import keygen, regeneratefg
from hawk.codec import decode_private
from hawk.poly import adjoint, isinvertible, poly_mul_ntt, poly_add, poly_sub, infnorm, l2norm

from expand import expand_challenge
import hashlib
# HAWK PARAMETERS
logn = 8
HAWK_SK_SIZE = 96
HAWK_PK_SIZE = 450

# GMWFS PARAMETERS
LAMBDA = 128
C = 7 
K = 22
ROUND = 84

MSG_HASH_SIZE = LAMBDA//4
CHLG_SIZE = LAMBDA//4

def poly_add_p(p0,p1,p):
    return [(p0 + p1)%p for (p0, p1) in zip(p0, p1)]

def owl_lite_keygen():
    print("Generating Key...")
    priv = [0] * C
    pub = [0] * C
    for i in range(C):
        priv[i*HAWK_SK_SIZE:(i+1)*HAWK_SK_SIZE], pub[i*HAWK_PK_SIZE:(i+1)*HAWK_PK_SIZE] =  keygen(logn)
    print("Key Generated!")
    return priv,pub

def owl_lite_sign(message, secret_key):
    # for keys in secret_key:
    #     _ , Fmod2, Gmod2, _ = decode_private(logn, key)
    # secret_keys = []
    print("Signing Message...")
    groups = []
    sets = []
    for i in range(ROUND):
        groups[i*HAWK_SK_SIZE:(i+1)*HAWK_SK_SIZE], sets[i*HAWK_PK_SIZE:(i+1)*HAWK_PK_SIZE] =  keygen(logn)
    
    set_bytes = [i.tobytes(4, byteorder='big') for i in sets]

    hash_message = hashlib.shake_256(message).digest(MSG_HASH_SIZE)
    hash_input = hash_message+b''.join(set_bytes)
    signed_message = hashlib.shake_256(hash_input).digest(CHLG_SIZE)
        
    chg_c, chg_nc, chg_val = expand_challenge(signed_message, CHLG_SIZE)
    
    unpacked_secret_keys = [0] * C
    nc_bases = [0] * K
    
    n = 1 << logn
    p = 8380417
    for r in range(K):
        secret_seed, b01, b11 = decode_private(groups[chg_nc[r]*HAWK_SK_SIZE:(chg_nc[r]+1)*HAWK_SK_SIZE]) 
        b00, b10 = regeneratefg(secret_seed.tobytes(), n)  
        if unpacked_secret_keys[chg_val[r]] == 0:
            secret_seed, F, G = decode_private(secret_key[chg_val[r]*HAWK_SK_SIZE:(chg_val[r]+1)*HAWK_SK_SIZE]) 
            f, g = regeneratefg(secret_seed.tobytes(), n)  
            a00 = G
            a01 = [(i*(-1))%p for i in F]
            a10 = [(i*(-1))%p for i in g]
            a11 = f
            unpacked_secret_keys[chg_val[r]] = [a00,a01,a10,a11]
            
        else:
            a00,a01,a10,a11 = unpacked_secret_keys[chg_val[r]]
        c00 = poly_add_p(poly_mul_ntt(a00, b00,p), poly_mul_ntt(a01, b10,p),p)
        c01 = poly_add_p(poly_mul_ntt(a00, b01,p), poly_mul_ntt(a01, b11,p),p)
        c10 = poly_add_p(poly_mul_ntt(a10, b00,p), poly_mul_ntt(a11, b10,p),p)    
        c11 = poly_add_p(poly_mul_ntt(a10, b01,p), poly_mul_ntt(a11, b11,p),p)
        nc_bases[r] = [[c00,c01,c10,c11]]
    
    c_bases = [0] * (HAWK_SK_SIZE*(ROUND-K))
    for r in range(ROUND-K):
        c_bases[r*HAWK_SK_SIZE:(r+1)*HAWK_SK_SIZE] = groups[chg_c[r]*HAWK_SK_SIZE:(chg_c[r]+1)*HAWK_SK_SIZE]

    return (signed_message, c_bases, nc_bases)

def owl_lite_verify():
    return

priv,pub = owl_lite_keygen()
message = "hello"
byte_list = message.encode('utf-8') 
signature = owl_lite_sign(byte_list, priv)
# print(alteq_verify(byte_list, pk, signature))

print(signature)

print(len(priv))
# print(len(pub))
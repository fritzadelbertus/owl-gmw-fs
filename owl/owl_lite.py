from gmwfs import Gen, Sign, Vrfy

from hawk.poly import adjoint, isinvertible, poly_mul_ntt, poly_add, poly_sub, infnorm, l2norm
from hawk.ntrugen.ntrugen_hawk import ntru_solve
import hawk.ntrugen
from hawk.params import PARAMS
import numpy as np
from hawk.rngcontext import RngContext, SHAKE256x4
import hashlib

# GMW-FS Requirements
# Parameters
# 1. c and C = 2^c              : integer                                           C is the number of group/set elements in the private/public key
# 2. r                          : integer                                           The number of rounds
# 3. hash_function              : (bitstring) => bitstring                          The hash function
# 4. messageInBits              : bitstring                                         The message to be signed

# Algebraic Requirements
# 1. action                     : (group element, set element) => set element       The group action function that returns an element of the set
# 2. group_sampler              : () => group element                               A function to sample and element of the group
# 3. set_sampler                : () => set element                                 A function to sample and element of the set
# 4. group_identity             : group element                                     The identity of the group
# 5. group_inverse              : (group_element) => group element                  A function that returns the inverse of the given group element
# 6. group_operator             : (group element, group element) => group element   The binary operator of the group

# Converters
# 1. setElementLengthInBits     : integer                                           The length of a set element in bitstring
# 2. groupElementLengthInBits   : integer                                           The length of a group element in bitstring
# 3. setToBits                  : (set element) => bitstring                        A function to convert a set element to bitstring
# 4. groupToBits                : (group element) => bitstring                      A function to convert a group element to bitstring
# 5. bitsToSet                  : (bitstring) => set                                A function to convert bitstring to the set element
# 6. bitsToGroup                : (bitstring) => group                              A function to convert bitstring to the group element

# Input
# messageInBits

# 3. hash_function              
# (bitstring) => bitstring                          
# The hash function
def bitsToBytes(bit_string):
    # Pad with zeros to make length multiple of 8
    length = len(bit_string)
    return int(bit_string, 2).to_bytes(length // 8, 'big')

def owl_hash(stringOfBits): 
    hash_obj = hashlib.sha3_256(bitsToBytes(stringOfBits))
    hash_bytes = hash_obj.digest()
    return ''.join(f'{byte:08b}' for byte in hash_bytes)

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

# # Algebraic Requirements
# # 1. action                     
# # (group element, set element) => set element       
# # The group action function that returns an element of the set
def action(B, Q):
    B_star = [[adjoint(B[0][0]),adjoint(B[1][0])],[adjoint(B[0][1]),adjoint(B[1][1])]]
    return matrix_mult(matrix_mult(B_star, Q), B)

# Regeneratefg for group_sample and set_sampler

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

# # 2. group_sampler              
# # () => group element                               
# # A function to sample and element of the group
# rng = np.random.default_rng(42)

def group_sampler(logn=8, rng=None): # Samples a random matrix
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
        return group_sampler(logn, rng)

    # Line 5: if ||(f,g)||2 <= 2*n*sigkrsec**2:
    if (l2norm(f) + l2norm(g)) <= (2 * n * (PARAMS(logn, "sigmakrsec") ** 2)):
        return group_sampler(logn, rng)

    try:
        # Line 13: r <- NTRUSolve(f, g, 1)
        # Line 16: (F, G) <- r
        F, G = ntru_solve(f, g)
    except ValueError:
        # Line 14&15: if r = \bot then restart
        return group_sampler(logn, rng)

    # Line 17: if infnorm( (F,G) ) > 127 then
    if infnorm(F) > 127 or infnorm(G) > 127:
        # Line 18: restart
        return group_sampler(logn, rng)

    # Line 28: return (priv, pub)
    p = 8380417
    f = [i%p for i in f]
    F = [i%p for i in F]
    g = [i%p for i in g]
    G = [i%p for i in G]
    return [[f,F],[g,G]]

# # 3. set_sampler
# # () => set element                                 
# # A function to sample and element of the set
def set_sampler(logn=8, rng=None):
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
        return set_sampler(logn, rng)

    fadj = adjoint(f)
    gadj = adjoint(g)

    # Line 5: if ||(f,g)||2 <= 2*n*sigkrsec**2:
    if (l2norm(f) + l2norm(g)) <= (2 * n * (PARAMS(logn, "sigmakrsec") ** 2)):
        return set_sampler(logn, rng)

    p = (1 << 16) + 1

    # Line 7: q00 <- f f* + g g*
    q00 = poly_add(poly_mul_ntt(f, fadj, p), poly_mul_ntt(g, gadj, p))

    # Line 8: (p1, p2) <- (2147473409, 2147389441)
    p1, p2 = (2147473409, 2147389441)

    # Line 9: if isInvertible(q00, p1) false or isInvertible(q00, p2) is false then
    if not isinvertible(q00, p1) or not isinvertible(q00, p2):
        # Line 10: restart
        return set_sampler(logn, rng)

    # Line 11: if (1/q00)[0] >= beta0 then
    invq00 = hawk.ntrugen.fft.inv_fft(q00)
    if invq00[0] >= PARAMS(logn, "beta0"):
        # Line 12: restart
        return set_sampler(logn, rng)

    try:
        # Line 13: r <- NTRUSolve(f, g, 1)
        # Line 16: (F, G) <- r
        F, G = ntru_solve(f, g)
    except ValueError:
        # Line 14&15: if r = \bot then restart
        return set_sampler(logn, rng)

    # Line 17: if infnorm( (F,G) ) > 127 then
    if infnorm(F) > 127 or infnorm(G) > 127:
        # Line 18: restart
        return set_sampler(logn, rng)

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
        return set_sampler(logn, rng)
    q00 = [i%p for i in q00]
    q01 = [i%p for i in q01]
    q10 = [i%p for i in q10]
    q11 = [i%p for i in q11]
    return [[q00,q01],[q10,q11]]

# # 4. group_identity             
# # group element                                     
# # The identity of the group
zeros256 = [0 for i in range(256)]
poly_unit = [0 for i in range(256)]
poly_unit[0] = 1
group_identity = [[poly_unit, zeros256],[zeros256, poly_unit]]

# # 5. group_inverse              
# # (group_element) => group element                  
# # A function that returns the inverse of the given group element
def group_inverse(A):
    p = 8380417
    a = A[0][0]
    minus_b = [(i*(-1))%p for i in A[0][1]]
    minus_c = [(i*(-1))%p for i in A[1][0]]
    d = A[1][1]
    return [[d, minus_b],[minus_c,a]]

# # 6. group_operator             
# # (group element, group element) => group element   
# # The binary operator of the group
def group_operator(A,B):
    return matrix_mult(B,A)

# # ===================================================================================================

# # Converters

def bytes_to_bits(byte_data):
    # Convert each byte to 8-bit binary, then concatenate
    return ''.join(f'{b:08b}' for b in byte_data)

def bits_to_bytes(bit_str):
    # Make sure length is multiple of 8
    assert len(bit_str) % 8 == 0, "Bit string length must be multiple of 8"
    return bytes(int(bit_str[i:i+8], 2) for i in range(0, len(bit_str), 8))

def encode_poly_matrix(mat):
    q=8380417
    """Encode a 2x2 polynomial matrix into a bytearray using fixed-width 23-bit per coefficient"""
    coeffs = []
    for row in mat:
        for poly in row:
            for c in poly:
                # For negative integers, convert to positive mod q
                c_mod = c % q
                coeffs.append(c_mod)
    
    # Pack 23-bit integers into bytes
    # We'll use 32-bit integers for simplicity (padding unused bits)
    encoded = bytearray()
    for c in coeffs:
        encoded += c.to_bytes(4, byteorder='big')  # 4 bytes = 32 bits
    return encoded

def decode_poly_matrix(encoded):
    q=8380417
    coeffs = []
    for i in range(0, len(encoded), 4):
        c = int.from_bytes(encoded[i:i+4], byteorder='big')
        coeffs.append(c % q)  # reduce modulo q just in case
    coeffs = np.array(coeffs, dtype=int)
    
    mat = np.zeros((2,2,256), dtype=int)
    idx = 0
    for i in range(2):
        for j in range(2):
            mat[i,j] = coeffs[idx:idx+256]
            idx += 256
    return mat

# # 1. setElementLengthInBits 
# # integer                                           
# # The length of a set element in bitstring
setElementLengthInBits = 2*2 * 256 * 32

# # 2. groupElementLengthInBits 
# # integer                                           
# # The length of a group element in bitstring
groupElementLengthInBits = 2*2 * 256 * 32

# # 3. setToBits                  
# # (set element) => bitstring                        
# # A function to convert a set element to bitstring
def setToBits(Q):
    return bytes_to_bits(encode_poly_matrix(Q))

# # 4. groupToBits                
# # (group element) => bitstring                      
# # A function to convert a group element to bitstring
def groupToBits(B):
    return bytes_to_bits(encode_poly_matrix(B))

# # 5. bitsToSet                  
# # (bitstring) => set                                
# # A function to convert bitstring to the set element
def bitsToSet(bitstring):
    return decode_poly_matrix(bits_to_bytes(bitstring))

# # 6. bitsToGroup                
# # (bitstring) => group                              
# # A function to convert bitstring to the group element
def bitsToGroup(bitstring):
    return decode_poly_matrix(bits_to_bytes(bitstring))


# OWL - GMWFS Construction
LAMBDA = 128
CHLG_SIZE = LAMBDA//4
C = 1
K = 22
ROUND = 1

def owl_Gen():
    # C                 : integer                                       Number of Set Elements in the private and public key
    # group_identity    : group element                                 The identity of the group
    # group_sampler     : () => group element                           A function to sample and element of the group
    # set_sampler       : () => set element                             A function to sample and element of the set
    # action            : (group element, set element) => set element   The group action function that returns an element of the set
    # setToBits         : (set element) => bitstring                    A function to convert a set element to bitstring
    # groupToBits       : (group element) => bitstring                  A function to convert a group element to bitstring

    #print("Generating Key...")
    B_i = group_sampler()
    public_key = []
    Q_C = set_sampler()
    
    public_key.append(action(B_i, Q_C))
    public_key.append(Q_C)
    
    privateKeyInBits = groupToBits(group_inverse(B_i))
    publicKeyInBits = ''.join([setToBits(s) for s in public_key])
    #print("Key Generated!")
    return publicKeyInBits, privateKeyInBits

def owl_Sign(privateKeyInBits, publicKeyInBits, messageInBits):
    # privateKeyInBits              : bitstring                                         private key in bits
    # publicKeyInBits               : bitstring                                         public key in bits
    # messageInBits                 : bitstring                                         message in bits
    # r                             : integer                                           number of rounds
    # c                             : integer                                           log2 of C (the num of private/public keys)
    # hash_function                 : (bitstring) => bitstring                          the hash function
    # group_sampler                 : () => group element                               A function to sample and element of the group
    # action                        : (group element, set element) => set element       The group action function that returns an element of the set
    # group_operator                : (group element, group element) => group element   The binary operator of the group
    # group_inverse                 : (group element) => group element                  A function that returns the inverse of the given group element
    # bitsToGroup                   : (bitstring) => group element                      A function to convert bitstring to the group element
    # bitsToSet                     : (bitstring) => set element                        A function to convert bitstring to the set element
    # setToBits                     : (set element) => bitstring                        A function to convert a set element to bitstring
    # groupToBits                   : (group element) => bitstring                      A function to convert a group element to bitstring
    # groupElemenentLengthInBits    : integer                                           The length of a group element in bitstring
    # setElementLengthInBits        : integer                                           The length of a set element in bitstring

    #print("Signing Message...")
    
    # =========================== Bit Operations ==================================

    # Split the private key and public key to set elements and group elements (still in bits)
    publicKeySetElementsInBits = [publicKeyInBits[i:i+setElementLengthInBits] for i in range(0, len(publicKeyInBits), setElementLengthInBits)]

    # Convert to actual group and set element
    public_key = [bitsToSet(element) for element in publicKeySetElementsInBits]
    private_key = bitsToGroup(privateKeyInBits)

    # =========================== Non Bit Operations ========================================
    h_i = group_sampler()
    t_i = action(h_i, public_key[1])

    # =========================== Bit Operations ============================================
    hash_input = messageInBits + setToBits(t_i)
    
    cha = hashlib.shake_256(bitsToBytes(hash_input)).digest(CHLG_SIZE)
    signed_message = ''.join(f'{byte:08b}' for byte in cha)

    f_i = group_operator(h_i, private_key)

    # ========================== Bit operations ===========================================
    sign = signed_message + groupToBits(f_i)

    #print("Message Signed!")
    return sign

def owl_Vrfy(publicKeyInBits, messageInBits, sign):
    # publicKeyInBits           : bitstring                                     public key in bits
    # messageInBits             : bitstring                                     message in bits
    # sign                      : bitstring                                     The signature in bits
    # r                         : integer                                       number of rounds
    # c                         : integer                                       log2 of C (the num of private/public keys)
    # setElementLengthInBits    : integer                                       The length of a set element in bitstring
    # groupElementLengthInBits  : integer                                       The length of a group element in bitstring
    # hash_function             : (bitstring) => bitstring                      The hash function
    # bitsToGroup               : (bitstring) => group element                  A function to convert bitstring to the group element
    # bitsToSet                 : (bitstring) => set element                    A function to convert bitstring to the set element
    # setToBits                 : (set element) => bistring                     A function to convert a set element to bitstring
    # action                    : (group element, set element) => set element   The group action function that returns an element of the set

    #print("Verifying Message...")

    # Convert public key in bits to public key in set elements (S)
    publicKeySetElementsInBits = [publicKeyInBits[i:i+setElementLengthInBits] for i in range(0, len(publicKeyInBits), setElementLengthInBits)]
    public_key = [bitsToSet(element) for element in publicKeySetElementsInBits]

    # Get the b_i as integer
    signed_message = sign[:CHLG_SIZE*8]
    cha = bytes(int(signed_message[i:i+8], 2) for i in range(0, len(signed_message), 8))

    # ============================== Non bit Operation
    right_part = sign[CHLG_SIZE*8:]
    f_i = bitsToGroup(right_part)

    t_i = action(f_i, public_key[0])
    # ================================ Bit Operations
    hash_input = messageInBits + setToBits(t_i)

    cha_2 = hashlib.shake_256(bitsToBytes(hash_input)).digest(CHLG_SIZE)
    if cha == cha_2:
        return 0
    else:
        return 1

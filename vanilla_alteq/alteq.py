from gmwfs import Gen, Sign, Vrfy
import hashlib
import numpy as np
import random
import math
import galois as gl

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

# Parameters
# 1. c and C = 2^c 
# integer                                           
# C is the number of group/set elements in the private/public key
c = 8
C = 2**c

# 2. r                          
# integer                                           
# The number of rounds
r = 32

# 3. hash_function              
# (bitstring) => bitstring                          
# The hash function
def bitsToBytes(bit_string):
    # Pad with zeros to make length multiple of 8
    length = len(bit_string)
    return int(bit_string, 2).to_bytes(length // 8, 'big')

def alteq_hash(stringOfBits): 
    hash_obj = hashlib.sha3_256(bitsToBytes(stringOfBits))
    hash_bytes = hash_obj.digest()
    return ''.join(f'{byte:08b}' for byte in hash_bytes)

# ==============================================================================================================

# Algebraic Requirements
q = 2**128
n = 4
gf = gl.GF(q, irreducible_poly = gl.irreducible_poly(2,128))

# 1. action                     
# (group element, set element) => set element       
# The group action function that returns an element of the set
def alteq_action(matrix, atf):
    det = gf(int(np.linalg.det(matrix)) % q)
    result = gf.Zeros(((n,n,n)))
    for k in range(n-1, -1, -1):
        for j in range(k-1, -1, -1):
            for i in range(j-1, -1, -1):
                result[k][i][j] = det*atf[k][i][j]
    return result

# 2. group_sampler              
# () => group element                               
# A function to sample and element of the group
rng = np.random.default_rng(42)

def group_sampler(): # Samples a random unimodular matrix
    matrix = gf.Identity(n)
    randRNG = np.random.default_rng(random.randint(0, 2**32 - 1))
    flag = True
    while flag:
        matrix = gf.Random((n, n), seed=randRNG)
        if np.linalg.det(matrix) != 0:
            flag = False
    return matrix

# 3. set_sampler
# () => set element                                 
# A function to sample and element of the set
def set_sampler(): # Samples a random hermitian matrix
    result = gf.Zeros(((n,n,n)))
    randRNG = np.random.default_rng(random.randint(0, 2**32 - 1))
    for k in range(n-1, -1, -1):
        for j in range(k-1, -1, -1):
            for i in range(j-1, -1, -1):
                result[k][i][j] = gf.Random(seed = randRNG)
    return result

# 4. group_identity             
# group element                                     
# The identity of the group
group_identity = gf.Identity(n)

# 5. group_inverse              
# (group_element) => group element                  
# A function that returns the inverse of the given group element
def group_inverse(matrix):
    return np.linalg.inv(matrix)

# 6. group_operator             
# (group element, group element) => group element   
# The binary operator of the group
def group_operator(matrix1, matrix2):
    return matrix1 @ matrix2

# ===================================================================================================

# Converters
# 1. setElementLengthInBits 
# integer                                           
# The length of a set element in bitstring
setElementLengthInBits = math.comb(n, 3) * q.bit_length()

# 2. groupElementLengthInBits 
# integer                                           
# The length of a group element in bitstring
groupElementLengthInBits = n*n * q.bit_length()

# 3. setToBits                  
# (set element) => bitstring                        
# A function to convert a set element to bitstring
def setToBits(atf):
    compressed = []
    for k in range(n-1, -1, -1):
        for j in range(k-1, -1, -1):
            for i in range(j-1, -1, -1):
                compressed.append(atf[k][i][j])
    res = ""
    num_bits = q.bit_length()
    for i in compressed:
        res += format(i, f'0{num_bits}b')
    return res

# 4. groupToBits                
# (group element) => bitstring                      
# A function to convert a group element to bitstring
def groupToBits(matrix):
    num_bits = q.bit_length()
    res = ""
    for i in range(n):
        for j in range(n):
            res += format(matrix[i][j], f'0{num_bits}b')
    return res

# 5. bitsToSet                  
# (bitstring) => set                                
# A function to convert bitstring to the set element
def bitsToSet(bitstring):
    res = []
    num_bits = q.bit_length()
    for j in [bitstring[i:i+num_bits] for i in range(0, len(bitstring), num_bits)]:
        res.append(int(j, 2))
    atf = gf.Zeros(((n,n,n)))
    c = 0
    for k in range(n-1, -1, -1):
        for j in range(k-1, -1, -1):
            for i in range(j-1, -1, -1):
                atf[k][i][j] = gf(res[c])    
                c += 1
    return atf

# 6. bitsToGroup                
# (bitstring) => group                              \
# A function to convert bitstring to the group element
def bitsToGroup(bitstring):
    res = []
    num_bits = q.bit_length()
    for j in [bitstring[i:i+num_bits] for i in range(0, len(bitstring), num_bits)]:
        res.append(int(j, 2))
    matrix = gf.Identity(n)
    c = 0
    for i in range(n):
        for j in range(n):
            matrix[i][j] = gf(res[c])
            c += 1
    return matrix


# ALTEQ - GMWFS Construction

def alteq_Gen():
    return Gen(C = C, group_identity=group_identity, 
               group_sampler=group_sampler, set_sampler=set_sampler,
               action=alteq_action, groupToBits=groupToBits,
               setToBits=setToBits)

def alteq_Sign(privateKeyInBits, publicKeyInBits, messageInBits):
    return Sign(privateKeyInBits=privateKeyInBits, publicKeyInBits=publicKeyInBits,
                messageInBits=messageInBits, setElementLengthInBits=setElementLengthInBits,
                groupElementLengthInBits=groupElementLengthInBits, r=r, c=c,
                hash_function=alteq_hash, group_sampler=group_sampler,
                action=alteq_action, group_operator=group_operator,
                group_inverse=group_inverse, bitsToGroup=bitsToGroup,
                bitsToSet=bitsToSet, setToBits=setToBits, groupToBits=groupToBits)

def alteq_Vrfy(publicKeyInBits, messageInBits, sign):
    return Vrfy(publicKeyInBits=publicKeyInBits, messageInBits=messageInBits, sign=sign, r=r, c=c,
                setElementLengthInBits=setElementLengthInBits, groupElementLengthInBits=groupElementLengthInBits,
                hash_function=alteq_hash, bitsToGroup=bitsToGroup, bitsToSet=bitsToSet,
                setToBits=setToBits, action=alteq_action)


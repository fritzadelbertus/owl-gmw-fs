from gmwfs import Gen, Sign, Vrfy
import hashlib
import numpy as np
import random
import galois

# WARNING! NOT USED ANYMORE

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

def owl_hash(stringOfBits): 
    hash_obj = hashlib.sha3_256(bitsToBytes(stringOfBits))
    hash_bytes = hash_obj.digest()
    return ''.join(f'{byte:08b}' for byte in hash_bytes)

# ==============================================================================================================

# Algebraic Requirements
q = 2**128
n = 2
gf = galois.GF(q, irreducible_poly = galois.irreducible_poly(2,128))

# 1. action                     
# (group element, set element) => set element       
# The group action function that returns an element of the set
def owl_action(unimod, hermit):
    return unimod.T @ hermit @ unimod

# 2. group_sampler              
# () => group element                               
# A function to sample and element of the group
rng = np.random.default_rng(42)

def group_sampler(steps=10): # Samples a random matrix
    A = gf(np.eye(n, dtype=int))
    randRNG = np.random.default_rng(random.randint(0, 2**32 - 1))
    for _ in range(steps):
        i, j = random.sample(range(n), 2)

        # Row addition: Ri ← Ri + k * Rj
        k = gf.Random(seed=randRNG)  # random nonzero field element
        if k != 0:
            A[i] += k * A[j]

        # Random row swap
        if random.random() < 0.25:
            A[[i, j]] = A[[j, i]]

    return A

# 3. set_sampler
# () => set element                                 
# A function to sample and element of the set
def set_sampler(steps=10): # Samples a random ATF
    A = gf(np.eye(n, dtype=int))
    randRNG = np.random.default_rng(random.randint(0, 2**32 - 1))
    for _ in range(steps):
        i, j = random.sample(range(n), 2)

        # Row addition: Ri ← Ri + k * Rj
        k = gf.Random(seed=randRNG)  # random nonzero field element
        if k != 0:
            A[i] += k * A[j]

        # Random row swap
        if random.random() < 0.25:
            A[[i, j]] = A[[j, i]]

    return A.T @ A

# 4. group_identity             
# group element                                     
# The identity of the group
group_identity = gf.Identity(n)

# 5. group_inverse              
# (group_element) => group element                  
# A function that returns the inverse of the given group element
def group_inverse(unimod):
    return np.linalg.inv(unimod)

# 6. group_operator             
# (group element, group element) => group element   
# The binary operator of the group
def group_operator(unimod1, unimod2):
    return unimod2 @ unimod1

# ===================================================================================================

# Converters
# 1. setElementLengthInBits 
# integer                                           
# The length of a set element in bitstring
setElementLengthInBits = n*n * q.bit_length()

# 2. groupElementLengthInBits 
# integer                                           
# The length of a group element in bitstring
groupElementLengthInBits = n*n * q.bit_length()

# 3. setToBits                  
# (set element) => bitstring                        
# A function to convert a set element to bitstring
def setToBits(hermit):
    num_bits = q.bit_length()
    res = ""
    for i in range(n):
        for j in range(n):
            res += format(hermit[i][j], f'0{num_bits}b')
    return res

# 4. groupToBits                
# (group element) => bitstring                      
# A function to convert a group element to bitstring
def groupToBits(unimod):
    num_bits = q.bit_length()
    res = ""
    for i in range(n):
        for j in range(n):
            res += format(unimod[i][j], f'0{num_bits}b')
    return res

# 5. bitsToSet                  
# (bitstring) => set                                
# A function to convert bitstring to the set element
def bitsToSet(bitstring):
    res = []
    num_bits = q.bit_length()
    for j in [bitstring[i:i+num_bits] for i in range(0, len(bitstring), num_bits)]:
        res.append(int(j, 2))
    hermit = gf.Identity(n)
    c = 0
    for i in range(n):
        for j in range(n):
            hermit[i][j] = gf(res[c])
            c += 1
    return hermit

# 6. bitsToGroup                
# (bitstring) => group                              
# A function to convert bitstring to the group element
def bitsToGroup(bitstring):
    res = []
    num_bits = q.bit_length()
    for j in [bitstring[i:i+num_bits] for i in range(0, len(bitstring), num_bits)]:
        res.append(int(j, 2))
    unimod = gf.Identity(n)
    c = 0
    for i in range(n):
        for j in range(n):
            unimod[i][j] = gf(res[c])
            c += 1
    return unimod


# OWL - GMWFS Construction

def owl_Gen():
    return Gen(C = C, group_identity=group_identity, 
               group_sampler=group_sampler, set_sampler=set_sampler,
               action=owl_action, groupToBits=groupToBits,
               setToBits=setToBits)

def owl_Sign(privateKeyInBits, publicKeyInBits, messageInBits):
    return Sign(privateKeyInBits=privateKeyInBits, publicKeyInBits=publicKeyInBits,
                messageInBits=messageInBits, setElementLengthInBits=setElementLengthInBits,
                groupElementLengthInBits=groupElementLengthInBits, r=r, c=c,
                hash_function=owl_hash, group_sampler=group_sampler,
                action=owl_action, group_operator=group_operator,
                group_inverse=group_inverse, bitsToGroup=bitsToGroup,
                bitsToSet=bitsToSet, setToBits=setToBits, groupToBits=groupToBits)

def owl_Vrfy(publicKeyInBits, messageInBits, sign):
    return Vrfy(publicKeyInBits=publicKeyInBits, messageInBits=messageInBits, sign=sign, r=r, c=c,
                setElementLengthInBits=setElementLengthInBits, groupElementLengthInBits=groupElementLengthInBits,
                hash_function=owl_hash, bitsToGroup=bitsToGroup, bitsToSet=bitsToSet,
                setToBits=setToBits, action=owl_action)

# Input
# def stringToBits(s):
#     byte_data = s.encode('utf-8')
#     bit_string = ''.join(f'{byte:08b}' for byte in byte_data)

#     return bit_string

# message = "Hello World!"
# messageInBits = stringToBits(message)

# res = 0
# for i in range(100):
#     test = np.linalg.det(group_sampler())
#     if test > gf(100):
#         res += 1
# print(res)
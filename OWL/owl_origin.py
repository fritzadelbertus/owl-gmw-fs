import hashlib

from OWL.expand import expand_challenge
from OWL.coder import encode_poly_matrix, decode_poly_matrix
from OWL.group_action import action, group_sampler, set_sampler, group_inverse, group_operator

from OWL.params import LOGN as logn
from OWL.params import ORIGIN

# ===================================================================================================
CHLG_SIZE = ORIGIN["LAMBDA"]//4
C = ORIGIN["C"]
K = ORIGIN["K"]
ROUND = ORIGIN["ROUND"]

setElementLengthInByte = 2*2 * (2**logn) * 4
groupElementLengthInByte = 2*2 * (2**logn) * 4


def owl_Gen():              
    #print("Generating Key...")
    private_key = []
    public_key = []
    Q_C = set_sampler()
    for i in range(C):
        B_i = group_sampler()
        public_key.append(action(B_i, Q_C))
        private_key.append(group_inverse(B_i))
    public_key.append(Q_C)
    
    privateKeyInBits = [encode_poly_matrix(g) for g in private_key]
    publicKeyInBits = [encode_poly_matrix(s) for s in public_key]
    #print("Key Generated!")
    return publicKeyInBits, privateKeyInBits

def owl_Sign(privateKey, publicKey, messageInByte):
    #print("Signing Message...")
 
    # Convert to actual group and set element
    public_key = [decode_poly_matrix(element, logn) for element in publicKey]
    private_key = [decode_poly_matrix(element, logn) for element in privateKey]

    
    h_i = [group_sampler() for i in range(ROUND)]
    t_i = [action(h_i[i], public_key[C]) for i in range(ROUND)]

    hash_input = messageInByte
    for ts in t_i:
        hash_input += encode_poly_matrix(ts)
    
    cha = hashlib.shake_256(hash_input).digest(CHLG_SIZE)
    chg_c, chg_nc, chg_val = expand_challenge(cha, CHLG_SIZE)

    f_i = [0]*ROUND

    for r in range(K):
        f_i[chg_nc[r]] = group_operator(h_i[chg_nc[r]],private_key[chg_val[r]])

    for r in range(ROUND-K):
        f_i[chg_c[r]] = h_i[chg_c[r]]

    sign = cha
    for f in f_i:
        sign += encode_poly_matrix(f)

    #print("Message Signed!")
    return sign

def owl_Vrfy(publicKey, messageInByte, sign):

    #print("Verifying Message...")

    public_key = [decode_poly_matrix(element, logn) for element in publicKey]

    # Get the b_i as integer
    cha = sign[:CHLG_SIZE]
    chg_c, chg_nc, chg_val = expand_challenge(cha, CHLG_SIZE)

    right_part = sign[CHLG_SIZE:]
    f_i_byte = [right_part[i:i+groupElementLengthInByte] for i in range(0, len(right_part), groupElementLengthInByte)]
    f_i = [decode_poly_matrix(f,logn) for f in f_i_byte]

    t_i = [0] * ROUND
    for r in range(K):
        t_i[chg_nc[r]] = action(f_i[chg_nc[r]], public_key[chg_val[r]])

    for r in range(ROUND-K):
        t_i[chg_c[r]] = action(f_i[chg_c[r]], public_key[C])

    hash_input = messageInByte
    for ts in t_i:
        hash_input += encode_poly_matrix(ts)

    cha_2 = hashlib.shake_256(hash_input).digest(CHLG_SIZE)
    if cha == cha_2:
        return 0
    else:
        return 1

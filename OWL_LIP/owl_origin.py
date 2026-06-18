import hashlib

from OWL_LIP.expand import expand_challenge
from OWL_LIP.group_action import action, group_sampler, set_sampler, group_identity, group_inverse, group_operator
from OWL_LIP.coder import encode_group_matrix, decode_group_matrix, encode_set_matrix, decode_set_matrix, INT_NBYTES

from OWL_LIP.params import ORIGIN

# ===================================================================================================
CHLG_SIZE = ORIGIN["LAMBDA"]//4
C = ORIGIN["C"]
K = ORIGIN["K"]
ROUND = ORIGIN["ROUND"]


groupElementLengthInByte = 4 * INT_NBYTES


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
    
    privateKey = [encode_group_matrix(g) for g in private_key]
    publicKey = [encode_set_matrix(s) for s in public_key]
    #print("Key Generated!")
    return publicKey, privateKey

def owl_Sign(privateKey, publicKey, messageInByte):
    #print("Signing Message...")
    # Convert to actual group and set element
    public_key = [decode_set_matrix(element) for element in publicKey]
    private_key = [decode_group_matrix(element) for element in privateKey]

    
    h_i = [group_sampler() for i in range(ROUND)]
    t_i = [action(h_i[i], public_key[C]) for i in range(ROUND)]
    hash_input = messageInByte
    for ts in t_i:
        hash_input += encode_set_matrix(ts)
    
    cha = hashlib.shake_256(hash_input).digest(CHLG_SIZE)
    chg_c, chg_nc, chg_val = expand_challenge(cha, CHLG_SIZE)

    f_i = [0]*ROUND

    for r in range(K):
        f_i[chg_nc[r]] = group_operator(h_i[chg_nc[r]],private_key[chg_val[r]])

    for r in range(ROUND-K):
        f_i[chg_c[r]] = h_i[chg_c[r]]

    sign = cha
    for f in f_i:
        sign += encode_group_matrix(f)

    #print("Message Signed!")
    return sign

def owl_Vrfy(publicKey, messageInByte, sign):

    #print("Verifying Message...")

    public_key = [decode_set_matrix(element) for element in publicKey]

    # Get the b_i as integer
    cha = sign[:CHLG_SIZE]
    chg_c, chg_nc, chg_val = expand_challenge(cha, CHLG_SIZE)

    right_part = sign[CHLG_SIZE:]
    f_i_byte = [right_part[i:i+groupElementLengthInByte] for i in range(0, len(right_part), groupElementLengthInByte)]
    f_i = [decode_group_matrix(f) for f in f_i_byte]

    t_i = [0] * ROUND
    for r in range(K):
        t_i[chg_nc[r]] = action(f_i[chg_nc[r]], public_key[chg_val[r]])

    for r in range(ROUND-K):
        t_i[chg_c[r]] = action(f_i[chg_c[r]], public_key[C])

    hash_input = messageInByte
    for ts in t_i:
        hash_input += encode_set_matrix(ts)

    cha_2 = hashlib.shake_256(hash_input).digest(CHLG_SIZE)
    return cha == cha_2

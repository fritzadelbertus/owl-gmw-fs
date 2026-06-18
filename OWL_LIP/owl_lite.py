import hashlib

from OWL.coder import encode_poly_matrix, decode_poly_matrix
from OWL.group_action import action, group_sampler, set_sampler, group_inverse, group_operator

from OWL.params import LOGN as logn
from OWL.params import LITE

# ==============================================================================================================

setElementLengthInBits = 2*2 * (2**logn) * 4
groupElementLengthInBits = 2*2 * (2**logn) * 4

CHLG_SIZE = LITE["LAMBDA"]//4

def owl_Gen():          

    #print("Generating Key...")
    B_i = group_sampler()
    public_key = []
    Q_C = set_sampler()
    
    public_key.append(action(B_i, Q_C))
    public_key.append(Q_C)
    
    privateKeyInBits = encode_poly_matrix(group_inverse(B_i))
    publicKeyInBits = [encode_poly_matrix(s) for s in public_key]
    #print("Key Generated!")
    return publicKeyInBits, privateKeyInBits

def owl_Sign(privateKey, publicKey, messageInByte):
    #print("Signing Message...")
    
    # Convert to actual group and set element
    public_key = [decode_poly_matrix(element, logn) for element in publicKey]
    private_key = decode_poly_matrix(privateKey, logn)

    h_i = group_sampler()
    t_i = action(h_i, public_key[1])

    hash_input = messageInByte + encode_poly_matrix(t_i)
    
    cha = hashlib.shake_256(hash_input).digest(CHLG_SIZE)

    f_i = group_operator(h_i, private_key)

    sign = cha + encode_poly_matrix(f_i)

    #print("Message Signed!")
    return sign

def owl_Vrfy(publicKey, messageInByte, sign):

    #print("Verifying Message...")

    public_key = [decode_poly_matrix(element, logn) for element in publicKey]

    # Get the b_i as integer
    cha = sign[:CHLG_SIZE]

    right_part = sign[CHLG_SIZE:]
    f_i = decode_poly_matrix(right_part,logn)

    t_i = action(f_i, public_key[0])
    hash_input = messageInByte + encode_poly_matrix(t_i)

    cha_2 = hashlib.shake_256(hash_input).digest(CHLG_SIZE)
    if cha == cha_2:
        return 0
    else:
        return 1

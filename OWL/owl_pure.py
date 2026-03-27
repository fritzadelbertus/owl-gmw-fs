import hashlib

from OWL.group_action import action, group_sampler, set_sampler, group_identity, group_inverse, group_operator
from OWL.coder import encode_poly_matrix, decode_poly_matrix

from OWL.params import LOGN as logn
from OWL.params import PURE
 
C = PURE["C"]+1
r = PURE["ROUND"]
c = C.bit_length() - 1

setElementLengthInBits = 2*2 * (2**logn) * 32
groupElementLengthInBits = 2*2 * (2**logn) * 32
# ==============================================================================================================

def bytes_to_bits(byte_data):
    # Convert each byte to 8-bit binary, then concatenate
    return ''.join(f'{b:08b}' for b in byte_data)

def bits_to_bytes(bit_str):
    # Make sure length is multiple of 8
    assert len(bit_str) % 8 == 0, "Bit string length must be multiple of 8"
    return bytes(int(bit_str[i:i+8], 2) for i in range(0, len(bit_str), 8))

def groupToBits(B):
    return bytes_to_bits(encode_poly_matrix(B))

def bitsToGroup(bitstring):
    return decode_poly_matrix(bits_to_bytes(bitstring), logn)

# GMW-FS

def owl_Gen():      
    #print("Generating Key...")
    private_key = [group_identity]
    for i in range(1,C):
        private_key.append(group_sampler())
    public_key = [set_sampler()]
    for i in range(1,C):
        public_key.append(action(private_key[i], public_key[0]))
    privateKeyInBits = [encode_poly_matrix(g) for g in private_key]
    publicKeyInBits = [encode_poly_matrix(s) for s in public_key]
    #print("Key Generated!")
    return publicKeyInBits, privateKeyInBits
    

def owl_Sign(privateKey, publicKey, messageInByte):
    #print("Signing Message...")
    
    public_key = [decode_poly_matrix(element, logn) for element in publicKey]
    private_key = [decode_poly_matrix(element, logn) for element in privateKey]

    h_i = [group_sampler() for i in range(r)]
    t_i = [action(h_i[i], public_key[0]) for i in range(r)]

    hash_input = messageInByte
    for ts in t_i:
        hash_input += encode_poly_matrix(ts)
    
    hash_byte = hashlib.sha3_256(hash_input).digest()
    hash_bits = bytes_to_bits(hash_byte)

    b_i = [int(hash_bits[i:i+c],2) for i in range(0, len(hash_bits), c)]
    f_i = [group_operator(h_i[i], group_inverse(private_key[b_i[i]])) for i in range(r)]

    sign = hash_bits
    for f in f_i:
        sign += groupToBits(f)

    #print("Message Signed!")
    return sign


def owl_Vrfy(publicKey, messageInByte, sign):
    
    #print("Verifying Message...")

    public_key = [decode_poly_matrix(element, logn) for element in publicKey]
    # Get the b_i as integer
    b_i_bits = sign[:r*c]
    b_i = [int(b_i_bits[i:i+c],2) for i in range(0, len(b_i_bits), c)]

    right_part = sign[r*c:]
    f_i_bits = [right_part[i:i+groupElementLengthInBits] for i in range(0, len(right_part), groupElementLengthInBits)]
    f_i = [bitsToGroup(f) for f in f_i_bits]
    t_i = [action(f_i[i], public_key[b_i[i]]) for i in range(r)]

    hash_input = messageInByte
    for ts in t_i:
        hash_input += encode_poly_matrix(ts)

    hash_byte = hashlib.sha3_256(hash_input).digest()
    hash_bits = bytes_to_bits(hash_byte)

    if hash_bits == sign[:r*c]:
        return 0
    else:
        return 1



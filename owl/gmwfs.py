# GMW-FS

def Gen(C, group_identity, group_sampler, set_sampler, action, groupToBits, setToBits):
    # C                 : integer                                       Number of Set Elements in the private and public key
    # group_identity    : group element                                 The identity of the group
    # group_sampler     : () => group element                           A function to sample and element of the group
    # set_sampler       : () => set element                             A function to sample and element of the set
    # action            : (group element, set element) => set element   The group action function that returns an element of the set
    # setToBits         : (set element) => bitstring                    A function to convert a set element to bitstring
    # groupToBits       : (group element) => bitstring                  A function to convert a group element to bitstring

    #print("Generating Key...")
    private_key = [group_identity]
    for i in range(1,C):
        private_key.append(group_sampler())
    public_key = [set_sampler()]
    for i in range(1,C):
        public_key.append(action(private_key[i], public_key[0]))
    
    privateKeyInBits = ''.join([groupToBits(g) for g in private_key])
    publicKeyInBits = ''.join([setToBits(s) for s in public_key])
    #print("Key Generated!")
    return publicKeyInBits, privateKeyInBits
    

def Sign(privateKeyInBits, publicKeyInBits, messageInBits, setElementLengthInBits, groupElementLengthInBits,
         r, c, hash_function, group_sampler, action, group_operator, group_inverse,
         bitsToGroup, bitsToSet, setToBits, groupToBits):
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
    privateKeySetElementsInBits = [privateKeyInBits[i:i+groupElementLengthInBits] for i in range(0, len(privateKeyInBits), groupElementLengthInBits)]

    # Convert to actual group and set element
    public_key = [bitsToSet(element) for element in publicKeySetElementsInBits]
    private_key = [bitsToGroup(element) for element in privateKeySetElementsInBits]

    # =========================== Non Bit Operations ========================================
    h_i = [group_sampler() for i in range(r)]
    t_i = [action(h_i[i], public_key[0]) for i in range(r)]

    # =========================== Bit Operations ============================================
    hash_input = messageInBits
    for ts in t_i:
        hash_input += setToBits(ts)
    
    hash_bits = hash_function(hash_input)

    # ========================== Non Bit Operations ========================================
    b_i = [int(hash_bits[i:i+c],2) for i in range(0, len(hash_bits), c)]
    f_i = [group_operator(h_i[i], group_inverse(private_key[b_i[i]])) for i in range(r)]

    # ========================== Bit operations ===========================================
    sign = hash_bits
    for f in f_i:
        sign += groupToBits(f)

    #print("Message Signed!")
    return sign


def Vrfy(publicKeyInBits, messageInBits, sign, r, c, 
         setElementLengthInBits, groupElementLengthInBits, hash_function, 
         bitsToGroup, bitsToSet, setToBits, action):
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
    b_i_bits = sign[:r*c]
    b_i = [int(b_i_bits[i:i+c],2) for i in range(0, len(b_i_bits), c)]

    # ============================== Non bit Operation
    right_part = sign[r*c:]
    f_i_bits = [right_part[i:i+groupElementLengthInBits] for i in range(0, len(right_part), groupElementLengthInBits)]
    f_i = [bitsToGroup(f) for f in f_i_bits]
    t_i = [action(f_i[i], public_key[b_i[i]]) for i in range(r)]

    # ================================ Bit Operations
    hash_input = messageInBits
    for ts in t_i:
        hash_input += setToBits(ts)

    hash_bits = hash_function(hash_input)

    if hash_bits == sign[:r*c]:
        return 0
    else:
        return 1


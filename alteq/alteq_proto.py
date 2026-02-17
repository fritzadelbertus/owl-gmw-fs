import numpy as np
import galois as gl
import matplotlib.pyplot as plt
import random
import math

atf_counter = 1
gl_counter = 1

rng = np.random.default_rng(42)
class ATF:
    def __init__(self, q, n, name=""):
        self.q = q
        self.n = n
        self.name = name
        self.poly = gl.irreducible_poly(2,128)
        self.gf = gl.GF(q, irreducible_poly = self.poly)
        self.tensor = self.gf.Zeros(((n,n,n)))
        self.comp_length = int(n*(n-1)*(n-2)/6)
        self.compressed = self.gf.Zeros(self.comp_length)

    def getCompressed(self):
        return self.compressed

    def randomFill(self, seed=rng):
        c = 0
        for k in range(self.n-1, -1, -1):
            for j in range(k-1, -1, -1):
                for i in range(j-1, -1, -1):
                    rand = self.gf.Random(seed = seed)
                    self.tensor[k][i][j] = rand
                    self.compressed[c] = rand
                    c += 1

    def fill(self, arrbits):
        c = 0
        for k in range(self.n-1, -1, -1):
            for j in range(k-1, -1, -1):
                for i in range(j-1, -1, -1):
                    self.tensor[k][i][j] = self.gf(arrbits[c])
                    self.compressed[c] = self.gf(arrbits[c])
                    c += 1

    def getTensor(self):
        return self.tensor
    
    def evaluate(self, t):
        res = self.gf(0)
        for k in range(self.n-1, -1, -1):
            for j in range(k-1, -1, -1):
                for i in range(j-1, -1, -1):
                    det = self.gf(int((np.linalg.det(np.array([
                        [t[0][i], t[1][i], t[2][i]],
                        [t[0][j], t[1][j], t[2][j]],
                        [t[0][k], t[1][k], t[2][k]]
                    ])))%self.q))
                    res = res + det*self.tensor[k][i][j]
        return int(res)
    
    def showTensor(self):
        tensor = np.array(self.tensor)
        num_matrices = self.n
        fig, axes = plt.subplots(1, num_matrices, figsize=(4 * num_matrices, 4))
        fig.canvas.manager.set_window_title(self.name)
        if num_matrices == 1:
            axes = [axes]

        for i, (ax, mat) in enumerate(zip(axes, tensor)):
            ax.set_axis_off()
            table = ax.table(cellText=mat, loc='center', cellLoc='center')
            table.set_fontsize(14)
            table.scale(2, 2)
            ax.set_title(f"{i+1}-th frontal slice")

        plt.tight_layout()
        plt.show()

        
    def showCompressed(self):
        print(self.compressed)
    
    def getBits(self):
        res = ""
        a = self.getCompressed()
        num_bits = self.q.bit_length()
        for i in a:
            res += format(i, f'0{num_bits}b')
        return res

    def acts(self, m):
        det = self.gf(int(np.linalg.det(m.getMatrix())) % self.q)
        newATF = ATF(self.q, self.n, f"({m.name} acts on {self.name})")
        c = 0
        for k in range(self.n-1, -1, -1):
            for j in range(k-1, -1, -1):
                for i in range(j-1, -1, -1):
                    newATF.tensor[k][i][j] = det*self.tensor[k][i][j]
                    newATF.compressed[c] = det*self.tensor[k][i][j]
                    c += 1
        return newATF

class GL:
    def __init__(self, q, n, name=""):
        self.q = q
        self.n = n
        self.name = name
        self.poly = gl.irreducible_poly(2,128)
        self.gf = gl.GF(q, irreducible_poly = self.poly)
        self.entries = self.gf.Identity(n)
    
    def randomFill(self, seed=rng):
        flag = True
        while flag:
            self.entries = self.gf.Random((self.n, self.n), seed=seed)
            if np.linalg.det(self.entries) != 0:
                flag = False
    
    def fill(self, arrbits):
        c = 0
        for i in range(self.n):
            for j in range(self.n):
                self.entries[i][j] = arrbits[c]
                c += 1
    
    def getMatrix(self):
        return self.entries
    
    def show(self):
        entries = np.array(self.entries)

        fig, ax = plt.subplots()
        fig.canvas.manager.set_window_title(self.name)
        ax.set_axis_off() 

        table = ax.table(cellText=entries, loc='center', cellLoc='center')

        table.scale(2, 2)  # Cell size: (width, height)
        table.set_fontsize(14)

        plt.tight_layout()
        plt.show()

    def getBits(self):
        m = self.getMatrix()
        num_bits = self.q.bit_length()
        res = ""
        for i in range(self.n):
            for j in range(self.n):
                res += format(m[i][j], f'0{num_bits}b')
        return res




def createRandomATF(q,n):
    global atf_counter
    atf = ATF(q,n, f"random_atf_{atf_counter}")
    randRNG = np.random.default_rng(random.randint(0, 2**32 - 1))
    atf.randomFill(randRNG)
    atf_counter += 1
    return atf

def createRandomGL(q,n):
    global gl_counter
    m = GL(q,n, f"random_gl_{gl_counter}")
    randRNG = np.random.default_rng(random.randint(0, 2**32 - 1))
    m.randomFill(randRNG)
    gl_counter += 1
    return m

# atf1 = createRandomATF(7,4)
# m = createRandomGL(7,4)
# tens = [
#     [1,2,3,4],
#     [3,2,1,2],
#     [4,3,2,1]
# ]


# atf2 = atf1.acts(m)

# atf1.showTensor()
# m.show()
# atf2.showTensor()

# evaltens = atf1.evaluate(tens)
# print(evaltens)

# m.show()


def bitsToMatrix(bits, q, n):
    res = []
    num_bits = q.bit_length()
    for j in [bits[i:i+num_bits] for i in range(0, len(bits), num_bits)]:
        res.append(int(j, 2))
    newMatrix = GL(q,n)
    newMatrix.fill(res)
    return newMatrix


# bitsToMatrix(m.getBits(), 7, 4).show()

def atfToBits(atf):
    res = ""
    a = atf.getCompressed()
    num_bits = atf.q.bit_length()
    for i in a:
        res += format(i, f'0{num_bits}b')
    return res

def bitsToAtf(bits, q, n):
    res = []
    num_bits = q.bit_length()
    for j in [bits[i:i+num_bits] for i in range(0, len(bits), num_bits)]:
        res.append(int(j, 2))
    newATF = ATF(q,n)
    newATF.fill(res)
    return newATF

# atfToBits(atf1)
# res = bitsToAtf(atfToBits(atf1), atf1.q, atf1.n)
# print(res.getTensor())

c = 8
r = 32
C = 2**c
q = 2**128
n = 4

def generateIdentity(q,n):
    GF = gl.GF(q)
    return GF(np.identity(n, dtype=int))


# def generatePrivateKey(C, q, n):
#     pub = [GL(q,n)]
#     for i in range(1,C):
#         pub.append(createRandomGL(q,n))
#     return pub

# private_key = generatePrivateKey(C,7,4)

################################################ KEY GENERATION

def generaterKey(C,q,n):
    print("Generating Key...")
    pri = [GL(q,n)]
    for i in range(1,C):
        pri.append(createRandomGL(q,n))
    pub = [createRandomATF(q,n)]
    for i in range(1,C):
        pub.append(pub[0].acts(pri[i]))
    print("Key Generated!")
    return pub, pri

pub,pri = generaterKey(C,q,n)


# for matrix in pri:
#     print(matrix.getMatrix())

# for atfs in pub:
#     print(atfs.getCompressed())

public_key = ""
for atfs in pub:
    public_key += atfs.getBits()

private_key = ""
for matrix in pri:
    private_key += matrix.getBits()

# print(public_key)
# print(private_key)

######################################## SIGNING
import hashlib

def stringToBits(s):
    byte_data = s.encode('utf-8')  # b'hello'

    # Step 2: Convert each byte to 8-bit binary and join
    bit_string = ''.join(f'{byte:08b}' for byte in byte_data)

    return bit_string

def bitsToBytes(bit_string):
    # Pad with zeros to make length multiple of 8
    length = len(bit_string)
    return int(bit_string, 2).to_bytes(length // 8, 'big')

def signMessage(private_key, public_key, message, r, c, q,n):
    print("Signing Message...")
    num_bits_priv = n*n * q.bit_length()
    private_keys = []

    num_bits_pub = math.comb(n, 3) * q.bit_length()
    public_keys = []
    for keys in [public_key[i:i+num_bits_pub] for i in range(0, len(public_key), num_bits_pub)]:
        public_keys.append(bitsToAtf(keys, q, n))

    for keys in [private_key[i:i+num_bits_priv] for i in range(0, len(private_key), num_bits_priv)]:
        private_keys.append(bitsToMatrix(keys, q, n))
    
    # print(public_keys)
    # for key in private_keys:
    #     print(key.getMatrix())

    h_i = []
    for i in range(r):
        h_i.append(createRandomGL(q,n))
    t_i = []
    for i in range(r):
        t_i.append(public_keys[0].acts(h_i[i]))
    hash_input = stringToBits(message)
    for ts in t_i:
        hash_input += ts.getBits()
    
    hash_obj = hashlib.sha3_256(bitsToBytes(hash_input))
    hash_bytes = hash_obj.digest()
    hash_bits = ''.join(f'{byte:08b}' for byte in hash_bytes)
    b_i = [int(hash_bits[i:i+c],2) for i in range(0, len(hash_bits), c)]
    f_i = []
    for i in range(r):
        # print(private_keys[b_i[i]].getMatrix().shape)
        t = GL(q,n)
        t.entries = h_i[i].getMatrix() @ np.linalg.inv(private_keys[b_i[i]].getMatrix())
        f_i.append(t)
    sign = hash_bits

    for f in f_i:
        sign += f.getBits()
    print("Message Signed!")
    return sign



# Input data (must be bytes)
message = "hello world"
sign = signMessage(private_key, public_key, message, r, c, q, n)

########################################## VERIFY

def verifyMessage(public_key, message, sign, r, c, q, n):
    print("Verifying Message...")
    num_bits_pub = math.comb(n, 3) * q.bit_length()
    public_keys = []
    for keys in [public_key[i:i+num_bits_pub] for i in range(0, len(public_key), num_bits_pub)]:
        public_keys.append(bitsToAtf(keys, q, n))

    b_i_bits = sign[:r*c]
    b_i = [int(b_i_bits[i:i+c],2) for i in range(0, len(b_i_bits), c)]

    f_i_length = n*n*q.bit_length()
    right_part = sign[r*c:]
    f_i_bits = [right_part[i:i+f_i_length] for i in range(0, len(right_part), f_i_length)]
    f_i = []
    for f in f_i_bits:
        # print(len(f))
        f_i.append(bitsToMatrix(f, q, n))
    t_i = []
    for i in range(r):
        t_i.append(public_keys[b_i[i]].acts(f_i[i]))

    hash_input = stringToBits(message)
    for ts in t_i:
        hash_input += ts.getBits()

    hash_obj = hashlib.sha3_256(bitsToBytes(hash_input))
    hash_bytes = hash_obj.digest()
    hash_bits = ''.join(f'{byte:08b}' for byte in hash_bytes)

    if hash_bits == sign[:r*c]:
        print("Verification Completed! Signature is Real!")
    else:
        print("Verification Completed! Signature is Fake!")

# if sign[-1] == "0":
#     fakesign = sign[:-1]+"1"
# else:
#     fakesign = sign[:-1]+"0"
# verifyMessage(public_key, message, sign, r, c, q, n)
# verifyMessage(public_key, message, fakesign, r, c, q, n)

# print(sign)
# print(fakesign)

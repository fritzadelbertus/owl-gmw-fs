# Test

from alteq import alteq_Gen, alteq_Sign, alteq_Vrfy

def stringToBits(s):
    byte_data = s.encode('utf-8')
    bit_string = ''.join(f'{byte:08b}' for byte in byte_data)

    return bit_string

message = "Hello World!"
messageInBits = stringToBits(message)

pub,pri = alteq_Gen()


sign = alteq_Sign(pri, pub, messageInBits)
# with open("public_key.txt", "w") as f:
#     f.write(pub)
# with open("private_key.txt", "w") as f:
#     f.write(pri)
# with open("signature1.txt", "w") as f:
#     f.write(sign)

if sign[-1] == "0":
    fakesign = sign[:-1]+"1"
else:
    fakesign = sign[:-1]+"0"
alteq_Vrfy(pub, messageInBits, sign)
alteq_Vrfy(pub, messageInBits, fakesign)


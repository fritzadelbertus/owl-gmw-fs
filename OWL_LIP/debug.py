# Test

from owl_pure import owl_Gen, owl_Sign, owl_Vrfy

def stringToByte(s):
    byte_data = s.encode('utf-8')

    return byte_data
message = "Hello World!"
messageInByte = stringToByte(message)

pub,pri = owl_Gen()

sign= owl_Sign(pri, pub, messageInByte)
# sign,t= owl_Sign(pri, pub, messageInBits)
# with open("public_key.txt", "w") as f:
#     f.write(pub)
# with open("private_key.txt", "w") as f:
#     f.write(pri)
# with open("signature1.txt", "w") as f:
#     f.write(sign)
print(len(sign))
# fakesign = bytearray(sign)   # mutable copy

# fakesign[-1] ^= 1            # flip the least significant bit

# fakesign = bytes(fakesign)   # back to bytes
# print(owl_Vrfy(pub, messageInBits, sign,t))
# print(owl_Vrfy(pub, messageInBits, fakesign,t))

print(owl_Vrfy(pub, messageInByte, sign))
# print(owl_Vrfy(pub, messageInByte, fakesign))


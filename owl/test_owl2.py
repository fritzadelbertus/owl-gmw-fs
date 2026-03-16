# Test

from owl_one import owl_Gen, owl_Sign, owl_Vrfy

def stringToBits(s):
    byte_data = s.encode('utf-8')
    bit_string = ''.join(f'{byte:08b}' for byte in byte_data)

    return bit_string

message = "Hello World!"
messageInBits = stringToBits(message)

pub,pri = owl_Gen()

pk_bits = len(pub)
sk_bits = len(pri)

pk_bytes = pk_bits // 8
sk_bytes = sk_bits // 8

print("Public key size:")
print(f"{pk_bits} bits ({pk_bytes} bytes)")

print("Private key size:")
print(f"{sk_bits} bits ({sk_bytes} bytes)")
# sign= owl_Sign(pri, pub, messageInBits)

sign,t= owl_Sign(pri, pub, messageInBits)
sig_bits = len(sign)
sig_bytes = sig_bits // 8

print("Signature size:")
print(f"{sig_bits} bits ({sig_bytes} bytes)")
# with open("public_key.txt", "w") as f:
#     f.write(pub)
# with open("private_key.txt", "w") as f:
#     f.write(pri)
# with open("signature1.txt", "w") as f:
#     f.write(sign)
print(len(sign))
if sign[-1] == "0":
    fakesign = sign[:-1]+"1"
else:
    fakesign = sign[:-1]+"0"
print(owl_Vrfy(pub, messageInBits, sign,t))
print(owl_Vrfy(pub, messageInBits, fakesign,t))

# print(owl_Vrfy(pub, messageInBits, sign))
# print(owl_Vrfy(pub, messageInBits, fakesign))


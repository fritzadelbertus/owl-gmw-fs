from alteq import alteq_keygen, alteq_sign, alteq_verify

pk,sk = alteq_keygen()
message = "hello"
byte_list = message.encode('utf-8') 
signature = alteq_sign(byte_list, sk)
print(alteq_verify(byte_list, pk, signature))

# print(pk)
# print(sk)
# print(type(pk[0][1]))
# print(len(pk[0])*4+32)
# print(len(sk))
# print([1,2,3,4,5,6][2:4])
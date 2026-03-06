from vanilla_alteq import alteq_keygen

pk,sk = alteq_keygen()
print(pk)
print(sk)
print(len(pk[0])*4+32)
print(len(sk))
# print([1,2,3,4,5,6][2:4])
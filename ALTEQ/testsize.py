import os
import pickle

from alteq import alteq_keygen, alteq_sign


def size_test():

    pk, sk = alteq_keygen()

    msg = os.urandom(32)

    sig = alteq_sign(msg, sk)

    pk_size = len(pickle.dumps(pk))
    sk_size = len(sk)
    sig_size = len(pickle.dumps(sig))

    print("Public key size :", pk_size, "bytes")
    print("Secret key size :", sk_size, "bytes")
    print("Signature size  :", sig_size, "bytes")


size_test()
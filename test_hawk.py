import time
import numpy as np
from HAWK.hawk import hawkkeygen, hawksign, hawkverify

def measure_sizes(logn):
    msg = np.random.randint(0,256,32,dtype=np.uint8)

    priv, pub = hawkkeygen(logn)

    sig = hawksign(logn, priv, msg)

    print("Private key size:", len(priv), "bytes")
    print("Public key size:", len(pub), "bytes")
    print("Signature size:", len(sig), "bytes")

def benchmark(logn, trials=20):

    msg = np.random.randint(0,256,32,dtype=np.uint8)

    # keygen
    t0 = time.time()
    for _ in range(trials):
        priv, pub = hawkkeygen(logn)
    keygen_time = (time.time() - t0)/trials

    # sign
    priv, pub = hawkkeygen(logn)
    t0 = time.time()
    for _ in range(trials):
        sig = hawksign(logn, priv, msg)
    sign_time = (time.time() - t0)/trials

    # verify
    sig = hawksign(logn, priv, msg)
    t0 = time.time()
    for _ in range(trials):
        hawkverify(logn, pub, msg, sig)
    verify_time = (time.time() - t0)/trials

    print("KeyGen:", keygen_time*1000, "ms")
    print("Sign:", sign_time*1000, "ms")
    print("Verify:", verify_time*1000, "ms")



def correctness_test(logn, trials=20):

    for _ in range(trials):

        msg = np.random.randint(0,256,32,dtype=np.uint8)

        priv, pub = hawkkeygen(logn)

        sig = hawksign(logn, priv, msg)

        assert hawkverify(logn, pub, msg, sig)

    print("All tests passed!")

for logn in [9]:
    print("logn =", logn)
    correctness_test(logn)
    measure_sizes(logn)
    benchmark(logn)
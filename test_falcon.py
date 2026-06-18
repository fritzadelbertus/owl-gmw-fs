import time
import numpy as np
from FALCON.falcon import Falcon

falcon = Falcon(512)
sk,vk = falcon.keygen()

def measure_sizes():
    msg = np.random.randint(0,256,32,dtype=np.uint8)
    m = msg.tobytes()

    priv, pub = falcon.keygen()

    sig = falcon.sign(priv, m)

    print("Private key size:", len(priv), "bytes")
    print("Public key size:", len(pub), "bytes")
    print("Signature size:", len(sig), "bytes")

def benchmark(trials=20):

    msg = np.random.randint(0,256,32,dtype=np.uint8)
    m = msg.tobytes()

    # keygen
    t0 = time.time()
    for _ in range(trials):
        priv, pub = falcon.keygen()
    keygen_time = (time.time() - t0)/trials

    # sign
    priv, pub = falcon.keygen()
    t0 = time.time()
    for _ in range(trials):
        sig = falcon.sign(priv, m)
    sign_time = (time.time() - t0)/trials

    # verify
    sig = falcon.sign(priv, m)
    t0 = time.time()
    for _ in range(trials):
        falcon.verify(pub, m, sig)
    verify_time = (time.time() - t0)/trials

    print("KeyGen:", keygen_time*1000, "ms")
    print("Sign:", sign_time*1000, "ms")
    print("Verify:", verify_time*1000, "ms")



def correctness_test(trials=20):

    for _ in range(trials):

        msg = np.random.randint(0,256,32,dtype=np.uint8)
        m = msg.tobytes()

        priv, pub = falcon.keygen()

        sig = falcon.sign(priv, m)

        assert falcon.verify(pub, m, sig)

    print("All tests passed!")

    
correctness_test()
measure_sizes()
benchmark()
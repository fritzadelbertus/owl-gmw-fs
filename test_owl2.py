import time
import numpy as np
from OWL.owl_one import owl_Gen, owl_Sign, owl_Vrfy

def measure_sizes():
    msg = np.random.randint(0,256,32,dtype=np.uint8)
    m = msg.tobytes()

    pub, priv = owl_Gen()

    sig = owl_Sign(priv, pub, m)

    print("Private key size:", len(priv), "bytes")
    print("Public key size:", len(pub), "bytes")
    print("Signature size:", len(sig), "bytes")

def benchmark(trials=20):

    msg = np.random.randint(0,256,32,dtype=np.uint8)
    m = msg.tobytes()
    # keygen
    t0 = time.time()
    for _ in range(trials):
        priv, pub = owl_Gen()
    keygen_time = (time.time() - t0)/trials

    # sign
    pub, priv = owl_Gen()
    t0 = time.time()
    for _ in range(trials):
        sig = owl_Sign(priv, pub, m)
    sign_time = (time.time() - t0)/trials

    # verify
    sig = owl_Sign(priv, pub, m)
    t0 = time.time()
    for _ in range(trials):
        owl_Vrfy(pub, m, sig)
    verify_time = (time.time() - t0)/trials

    print("KeyGen:", keygen_time*1000, "ms")
    print("Sign:", sign_time*1000, "ms")
    print("Verify:", verify_time*1000, "ms")



def correctness_test(trials=20):

    for _ in range(trials):

        msg = np.random.randint(0,256,32,dtype=np.uint8)
        m = msg.tobytes()
        pub, priv = owl_Gen()

        sig = owl_Sign(priv, pub, m)

        assert owl_Vrfy(pub, m, sig)

    print("All tests passed!")



correctness_test()
measure_sizes()
benchmark()
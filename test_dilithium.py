import time

import numpy as np

from DIL.dilithium import Dilithium2


def measure_sizes():
    msg = np.random.randint(0, 256, 32, dtype=np.uint8)
    m = msg.tobytes()

    pub, priv = Dilithium2.keygen()

    sig = Dilithium2.sign(priv, m)

    print("Private key size:", len(priv), "bytes")
    print("Public key size:", len(pub), "bytes")
    print("Signature size:", len(sig), "bytes")


def benchmark(trials=100):

    msg = np.random.randint(0, 256, 32, dtype=np.uint8)
    m = msg.tobytes()

    # keygen
    t0 = time.time()
    for _ in range(trials):
        pub, priv = Dilithium2.keygen()
    keygen_time = (time.time() - t0) / trials

    # sign
    pub, priv = Dilithium2.keygen()
    t0 = time.time()
    for _ in range(trials):
        sig = Dilithium2.sign(priv, m)
    sign_time = (time.time() - t0) / trials

    # verify
    sig = Dilithium2.sign(priv, m)
    t0 = time.time()
    for _ in range(trials):
        Dilithium2.verify(pub, m, sig)
    verify_time = (time.time() - t0) / trials

    print("KeyGen:", keygen_time * 1000, "ms")
    print("Sign:", sign_time * 1000, "ms")
    print("Verify:", verify_time * 1000, "ms")


def correctness_test(trials=100):

    for _ in range(trials):
        msg = np.random.randint(0, 256, 32, dtype=np.uint8)
        m = msg.tobytes()

        pub, priv = Dilithium2.keygen()

        sig = Dilithium2.sign(priv, m)

        assert Dilithium2.verify(pub, m, sig)

    print("All tests passed!")


correctness_test()
measure_sizes()
benchmark()

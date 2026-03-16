import time
import tracemalloc
import os

from owl_one import owl_Gen, owl_Sign, owl_Vrfy


TEST_ROUNDS = 5


# ===============================
# Utilities
# ===============================

def bytes_to_bits(data):
    return ''.join(f'{b:08b}' for b in data)


def memory_snapshot():
    current, peak = tracemalloc.get_traced_memory()
    return current / (1024*1024), peak / (1024*1024)


# ===============================
# Correctness Test
# ===============================

def test_correctness():

    print("=== Correctness Test ===")

    pk, sk = owl_Gen()

    msg = os.urandom(32)
    msg_bits = bytes_to_bits(msg)

    sig = owl_Sign(sk, pk, msg_bits)

    result = owl_Vrfy(pk, msg_bits, sig)

    if result == 0:
        print("PASS: signature verified")
    else:
        print("FAIL: signature rejected")


# ===============================
# Forgery Tests
# ===============================

def test_forgery():

    print("\n=== Forgery Tests ===")

    pk, sk = owl_Gen()

    msg = os.urandom(32)
    msg_bits = bytes_to_bits(msg)

    sig = owl_Sign(sk, pk, msg_bits)


    # 1️⃣ Message modification
    forged_msg = bytearray(msg)
    forged_msg[0] ^= 1
    forged_msg_bits = bytes_to_bits(bytes(forged_msg))

    res = owl_Vrfy(pk, forged_msg_bits, sig)

    print("Message forgery:", "PASS" if res == 1 else "FAIL")


    # 2️⃣ Signature bit modification
    forged_sig = list(sig)
    forged_sig[0] = '1' if forged_sig[0] == '0' else '0'
    forged_sig = ''.join(forged_sig)

    res = owl_Vrfy(pk, msg_bits, forged_sig)

    print("Signature forgery:", "PASS" if res == 1 else "FAIL")


    # 3️⃣ Public key modification
    forged_pk = list(pk)
    forged_pk[0] = '1' if forged_pk[0] == '0' else '0'
    forged_pk = ''.join(forged_pk)

    res = owl_Vrfy(forged_pk, msg_bits, sig)

    print("Public key forgery:", "PASS" if res == 1 else "FAIL")


# ===============================
# Speed Benchmark
# ===============================

def benchmark():

    print("\n=== Speed Benchmark ===")

    keygen_time = 0
    sign_time = 0
    verify_time = 0

    for _ in range(TEST_ROUNDS):

        msg = os.urandom(32)
        msg_bits = bytes_to_bits(msg)

        t0 = time.perf_counter()
        pk, sk = owl_Gen()
        keygen_time += time.perf_counter() - t0


        t0 = time.perf_counter()
        sig = owl_Sign(sk, pk, msg_bits)
        sign_time += time.perf_counter() - t0


        t0 = time.perf_counter()
        owl_Vrfy(pk, msg_bits, sig)
        verify_time += time.perf_counter() - t0


    print(f"Keygen avg: {keygen_time/TEST_ROUNDS:.4f}s")
    print(f"Sign   avg: {sign_time/TEST_ROUNDS:.4f}s")
    print(f"Verify avg: {verify_time/TEST_ROUNDS:.4f}s")


# ===============================
# Memory Usage Test
# ===============================

def memory_test():

    print("\n=== Memory Usage ===")

    tracemalloc.start()

    pk, sk = owl_Gen()

    msg = os.urandom(32)
    msg_bits = bytes_to_bits(msg)

    sig = owl_Sign(sk, pk, msg_bits)

    owl_Vrfy(pk, msg_bits, sig)

    current, peak = memory_snapshot()

    print(f"Current memory: {current:.2f} MB")
    print(f"Peak memory: {peak:.2f} MB")

    tracemalloc.stop()


# ===============================
# Main
# ===============================

if __name__ == "__main__":

    test_correctness()

    test_forgery()

    benchmark()

    memory_test()
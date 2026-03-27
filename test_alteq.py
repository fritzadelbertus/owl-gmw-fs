import time
import tracemalloc
import copy
import os

from ALTEQ.alteq import alteq_keygen, alteq_sign, alteq_verify


TEST_ROUNDS = 1


def memory_snapshot():
    current, peak = tracemalloc.get_traced_memory()
    return current / (1024*1024), peak / (1024*1024)


def test_correctness():
    print("=== Correctness Test ===")

    pk, sk = alteq_keygen()

    msg = os.urandom(32)

    sig = alteq_sign(msg, sk)

    result = alteq_verify(msg, pk, sig)

    if result == 0:
        print("PASS: signature verified")
    else:
        print("FAIL: signature rejected")


def test_forgery():
    print("\n=== Forgery Tests ===")

    pk, sk = alteq_keygen()
    msg = os.urandom(32)

    sig = alteq_sign(msg, sk)

    # Modify message
    forged_msg = bytearray(msg)
    forged_msg[0] ^= 1

    res = alteq_verify(bytes(forged_msg), pk, sig)

    print("Message forgery:", "PASS" if res == 1 else "FAIL")

    # Modify signature hash
    forged_sig = list(sig)
    forged_sig[0] = bytearray(forged_sig[0])
    forged_sig[0][0] ^= 1
    forged_sig[0] = bytes(forged_sig[0])

    res = alteq_verify(msg, pk, tuple(forged_sig))

    print("Hash forgery:", "PASS" if res == 1 else "FAIL")

    # Modify matrices
    forged_sig = list(copy.deepcopy(sig))
    forged_sig[2][0] ^= 1

    res = alteq_verify(msg, pk, tuple(forged_sig))

    print("Matrix forgery:", "PASS" if res == 1 else "FAIL")


def benchmark():
    print("\n=== Speed Benchmark ===")

    keygen_time = 0
    sign_time = 0
    verify_time = 0

    # --- Keygen benchmark ---
    for _ in range(TEST_ROUNDS):
        t0 = time.perf_counter()
        alteq_keygen()
        keygen_time += time.perf_counter() - t0

    # --- Sign/Verify benchmark ---
    pk, sk = alteq_keygen()

    messages = [os.urandom(32) for _ in range(TEST_ROUNDS)]

    for msg in messages:
        t0 = time.perf_counter()
        sig = alteq_sign(msg, sk)
        sign_time += time.perf_counter() - t0

        t0 = time.perf_counter()
        alteq_verify(msg, pk, sig)
        verify_time += time.perf_counter() - t0

    print(f"Keygen avg: {keygen_time/TEST_ROUNDS:.6f}s")
    print(f"Sign   avg: {sign_time/TEST_ROUNDS:.6f}s")
    print(f"Verify avg: {verify_time/TEST_ROUNDS:.6f}s")

    # print(f"Keygen ops/sec: {TEST_ROUNDS/keygen_time:.2f}")
    # print(f"Sign   ops/sec: {TEST_ROUNDS/sign_time:.2f}")
    # print(f"Verify ops/sec: {TEST_ROUNDS/verify_time:.2f}")

if __name__ == "__main__":

    test_correctness()

    # test_forgery()

    benchmark()

import time
import tracemalloc
import os

from OWL.owl_pure import owl_Gen, owl_Sign, owl_Vrfy


TEST_ROUNDS = 1
# See params.py to modify OWL parameters

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

    sig = owl_Sign(sk, pk, msg)

    result = owl_Vrfy(pk, msg, sig)

    if result == 0:
        print("PASS: signature verified")
    else:
        print("FAIL: signature rejected")


# ===============================
# Forgery Tests
# ===============================

# def test_forgery():

#     print("\n=== Forgery Tests ===")

#     pk, sk = owl_Gen()

#     msg = os.urandom(32)

#     sig = owl_Sign(sk, pk, msg)


#     # Message modification
#     forged_msg = bytearray(msg)
#     forged_msg[0] ^= 1

#     res = owl_Vrfy(pk, bytes(forged_msg), sig)

#     print("Message forgery:", "PASS" if res == 1 else "FAIL")


#     # Signature bit modification
#     forged_sig = bytearray(sig)   # mutable copy
#     forged_sig[-1] ^= 1     

#     res = owl_Vrfy(pk, msg, bytes(forged_sig))

#     print("Signature forgery:", "PASS" if res == 1 else "FAIL")


    #  Public key modification
    # forged_pk = list(pk)
    # forged_pk[0] = '1' if forged_pk[0] == '0' else '0'
    # forged_pk = ''.join(forged_pk)

    # res = owl_Vrfy(forged_pk, msg_bits, sig)

    # print("Public key forgery:", "PASS" if res == 1 else "FAIL")


# ===============================
# Speed Benchmark
# ===============================
def benchmark():

    print("\n=== Speed Benchmark ===")

    keygen_time = 0
    sign_time = 0
    verify_time = 0

    # =========================
    # Keygen benchmark
    # =========================
    for _ in range(TEST_ROUNDS):
        t0 = time.perf_counter()
        owl_Gen()
        keygen_time += time.perf_counter() - t0

    # =========================
    # Sign / Verify benchmark
    # =========================
    pk, sk = owl_Gen()

    messages = [os.urandom(32) for _ in range(TEST_ROUNDS)]

    # --- Sign ---
    signatures = []
    for msg in messages:
        t0 = time.perf_counter()
        sig = owl_Sign(sk, pk, msg)   
        sign_time += time.perf_counter() - t0
        signatures.append(sig)

    # --- Verify ---
    for msg, sig in zip(messages, signatures):
        t0 = time.perf_counter()
        owl_Vrfy(pk, msg, sig)       
        verify_time += time.perf_counter() - t0

    print(f"Keygen avg: {keygen_time/TEST_ROUNDS:.6f}s")
    print(f"Sign   avg: {sign_time/TEST_ROUNDS:.6f}s")
    print(f"Verify avg: {verify_time/TEST_ROUNDS:.6f}s")

    # Throughput
    # print(f"Keygen ops/sec: {TEST_ROUNDS/keygen_time:.2f}")
    # print(f"Sign   ops/sec: {TEST_ROUNDS/sign_time:.2f}")
    # print(f"Verify ops/sec: {TEST_ROUNDS/verify_time:.2f}")

# ===============================
# Main
# ===============================

if __name__ == "__main__":

    test_correctness()

    # test_forgery()

    benchmark()

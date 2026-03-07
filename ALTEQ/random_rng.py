import hashlib

def shake_rng(seed: bytes):
    """Return a generator of bytes from SHAKE-256."""
    shake = hashlib.shake_256(seed)
    i = 0
    while True:
        # generate chunks of 64 bytes at a time
        chunk = shake.digest(64)
        for b in chunk:
            yield b

def get_random_value(rng_gen, max_value: int) -> int: 
    """ Generate a random integer in [0, max_value) using SHAKE-256 output. 
    Combines 8 bytes to produce a 64-bit integer, then reduce modulo max_value. """ 
    val = 0 
    for _ in range(8): 
        val = (val << 8) | next(rng_gen) 
    return val % max_value

def deterministic_sample(rng_gen, n: int, k: int) -> list[int]:
    """
    Deterministically sample k unique indices from range(n) using SHAKE-256.
    Equivalent to random.sample(range(n), k) but deterministic.
    """
    indices = list(range(n))

    # Fisher-Yates shuffle for first k elements
    for i in range(k):
        # get a random index in [i, n-1]
        rand_byte = 0
        for _ in range(8):
            rand_byte = (rand_byte << 8) | next(rng_gen)
        j = i + (rand_byte % (n - i))
        # swap indices[i] and indices[j]
        indices[i], indices[j] = indices[j], indices[i]

    # return first k shuffled indices
    return indices[:k]
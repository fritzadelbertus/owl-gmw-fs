import numpy as np

from rngcontext import RngContext


def get_random_value(rng_gen: RngContext, max_value: int) -> int:
    # Generate a random integer in [0, max_value) using SHAKE-256 output.
    random_bytes = rng_gen.random(8)
    value = int.from_bytes(random_bytes.tobytes(), byteorder="big")
    return value % max_value


def deterministic_sample(rng_gen: RngContext, n: int, k: int) -> list[int]:
    # Deterministically sample k unique indices from range(n).
    indices = list(range(n))

    # Fisher-Yates shuffle for first k elements
    for i in range(k):
        j = i + get_random_value(rng_gen, n - i)
        indices[i], indices[j] = indices[j], indices[i]

    return indices[:k]


def expand_challenge(
    seed: bytes, c: int, k: int, r: int, seed_size: int
) -> tuple[list[int], list[int], list[int]]:
    # Generate two list of index and a list of key index to sign h_i
    if seed_size < 32:
        raise ValueError("expand_challenge: seed size to low")
    seed_array = np.frombuffer(seed, dtype=np.uint8)
    rng_gen = RngContext(seed_array)

    if r - k < k:
        chg = [0] * r
        # pick ROUND-K coefficients to be C
        for i in deterministic_sample(rng_gen, r, r - k):
            chg[i] = c
        # fill remaining coefficients with values < C
        for i in range(r):
            if chg[i] == 0:
                chg[i] = get_random_value(rng_gen, c)
    else:
        # initialize all coefficients to C
        chg = [c] * r
        # pick K coefficients to be < C
        for i in deterministic_sample(rng_gen, r, k):
            chg[i] = get_random_value(rng_gen, c)

    # separate outputs
    chg_c = [i for i, v in enumerate(chg) if v == c]
    chg_nc = [i for i, v in enumerate(chg) if v < c]
    chg_val = [v for v in chg if v < c]

    return chg_c, chg_nc, chg_val

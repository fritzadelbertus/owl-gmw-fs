import hashlib
import os
from params import N, PRIME, LEN, ROUND, K, C
from random_rng import shake_rng, get_random_value, deterministic_sample

# SUBJECT TO CHANGES
# STORAGE SIZE IN expand_columns AND expand_atf IS NOT DETERMINED YET


def random_seed(seed_size: int) -> bytes:
    r = os.urandom(seed_size)
    return r

def expand_seeds(seed: bytes, n_seeds: int, out_seed_size: int) -> list[bytes]:
    """
    Expand a seed into several seeds with specific sizes using SHAKE256.

    Parameters:
    - seed: bytes, master seed
    - n_seeds: int, number of seeds to expand
    - out_seed_size: int, size of each seed

    Returns:
    - seeds: list of bytes, expanded seeds
    """
    total_size = n_seeds * out_seed_size
    
    shake = hashlib.shake_256(seed)
    expanded = shake.digest(total_size)

    seeds = [
        expanded[i*out_seed_size:(i+1)*out_seed_size]
        for i in range(n_seeds)
    ]

    return seeds

def expand_columns(seed: bytes, storage_size: int = 256) -> list[int]:
    """
    Expand a seed into an N x N column matrix modulo PRIME using SHAKE256.

    Parameters:
    - seed: bytes, master seed
    - storage_size: int, number of 32-bit integers generated per SHAKE call

    Returns:
    - col: N x N list of integers flatten
    """
    # initialize matrix
    col = [0] * (N * N)

    # storage buffer (like 'storage' in C code)
    storage = []
    k = 0  # index into storage
    i = 0  # row index
    j = 0  # column index

    # initial hash to fill storage
    shake = hashlib.shake_256(seed)
    # generate storage_size 32-bit integers
    raw_bytes = shake.digest(storage_size * 4)
    storage = [int.from_bytes(raw_bytes[4*r:4*(r+1)], 'little') for r in range(storage_size)]

    while i < N:
        if k >= storage_size:
            # rehash the storage buffer to refresh pseudo-random numbers
            shake = hashlib.shake_256(bytes([x & 0xFF for x in storage]))
            raw_bytes = shake.digest(storage_size * 4)
            storage = [int.from_bytes(raw_bytes[4*r:4*(r+1)], 'little') for r in range(storage_size)]
            k = 0

        col[i*N+j] = storage[k]
        k += 1

        # apply constraints: val < PRIME, diagonal non-zero
        if col[i*N+j] < PRIME and (i != j or col[i*N+j] > 0):
            j += 1
            if j == N:
                j = 0
                i += 1

    return col

def expand_atfs(seed: bytes, nb_copies: int, vec_size: int, storage_size: int = 256) -> list[int]:
    """
    Expand a seed into a vectorized ATFs (<PRIME)
    
    Parameters:
    - seed: bytes, seed for SHAKE256
    - LEN: int, number of positions to fill
    - PRIME: int, modulo prime
    - nb_copies: int, number of copies per position
    - vec_size: int, stride between vectors (number of ATF copies)
    - storage_size: int, how many integers to generate per hash
    """
    atf = [0] * (LEN * vec_size)  # output buffer
    storage = []
    k = 0  # index in storage
    i = 0  # position in atf
    
    # initial pseudo-random buffer
    shake = hashlib.shake_256(seed)
    raw_bytes = shake.digest(storage_size * 4)
    storage = [int.from_bytes(raw_bytes[4*r:4*(r+1)], 'little') for r in range(storage_size)]
    
    j = 0
    while i < LEN:
        if j >= storage_size:
            # refresh storage by rehashing
            shake = hashlib.shake_256(bytes([x & 0xFF for x in storage]))
            raw_bytes = shake.digest(storage_size * 4)
            storage = [int.from_bytes(raw_bytes[4*r:4*(r+1)], 'little') for r in range(storage_size)]
            j = 0

        val = storage[j]
        j += 1

        if val < PRIME:
            # copy value nb_copies times with stride vec_size
            for r in range(nb_copies):
                atf[i * vec_size + r] = val
            i += 1

    return atf

def expand_challenge(seed, seed_size):
    if seed_size < 32:
        raise Exception("expand_challenge Error: seed size to low")
    
    rng_gen = shake_rng(seed)
    
    chg = [0] * ROUND
    if ROUND - K < K:
        # pick ROUND-K coefficients to be C
        for r in deterministic_sample(rng_gen, ROUND, ROUND-K):
            chg[r] = C
        # fill remaining coefficients with values < C
        for i in range(ROUND):
            if chg[i] == 0:
                chg[i] = get_random_value(rng_gen, C)
    else:
        # initialize all coefficients to C
        chg = [C] * ROUND
        # pick K coefficients to be < C
        for r in deterministic_sample(rng_gen, ROUND, K):
            chg[r] = get_random_value(rng_gen, C)

    # separate outputs
    chg_c = [i for i, v in enumerate(chg) if v == C]
    chg_nc = [i for i, v in enumerate(chg) if v < C]
    chg_val = [v for v in chg if v < C]

    return chg_c, chg_nc, chg_val
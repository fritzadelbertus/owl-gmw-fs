from params import (SK_SEED_SIZE, MAT_SK_SEED_SIZE,
    C, K, N, ROUND, PK_SEED_SIZE, MSG_HASH_SIZE, USE_SALT, 
    SALT_SIZE, SIG_SEED_SIZE, EXPCOL_SIG_SEED_SIZE,
    CHLG_SIZE)

from expand import (random_seed, expand_seeds, 
    expand_columns, expand_atfs, expand_challenge)

from atf import inverting_on_atf, acting_on_atfs

from matrix import columns_matrix, columns_decomposition

import hashlib


def alteq_keygen():

    secret_key = random_seed(SK_SEED_SIZE)
    # dummy key
    # secret_key = b'\xa1\xf3\x8c\x9d\x02\xaf\x11\x5e\xde\x7a\x88\xbb\xef\x01\x9d\x30\x4c\x7e\x12\xab\x56\x9f\xcd\xef\x02\x11\x33\x45\x67\x89\xab\xcd\xef'

    seeds = expand_seeds(secret_key, 1, (MAT_SK_SEED_SIZE*C) + PK_SEED_SIZE)[0]

    columns = [x for r in range(C) for x in expand_columns(seeds[r*MAT_SK_SEED_SIZE:(r+1)*MAT_SK_SEED_SIZE])]
    
    atfs = expand_atfs(seeds[C*MAT_SK_SEED_SIZE:], C, C)
    atfs = inverting_on_atf(atfs, columns)
    
    public_key = (atfs,seeds[C*MAT_SK_SEED_SIZE:])
    return (public_key, secret_key)

def alteq_sign(message:list[bytes], secret_key: bytes):
    if len(message) <= 0:
        raise Exception("can't sign empty message")
        return
    
    seeds_sk = expand_seeds(secret_key, 1, (MAT_SK_SEED_SIZE*C) + PK_SEED_SIZE)[0]

    atfs = expand_atfs(seeds_sk[C*MAT_SK_SEED_SIZE:], ROUND, ROUND)

    hash_message = hashlib.shake_256(message).digest(MSG_HASH_SIZE)

    success = 0

    while (success == 0):
        if USE_SALT:
            seed = random_seed(SIG_SEED_SIZE)
            temp_seeds = expand_seeds(seed, ROUND, SIG_SEED_SIZE)
            salt = random_seed(SALT_SIZE)
            seeds = [0]*ROUND
            for r in range(ROUND):
                seeds[r] = (temp_seeds[r]+salt+bytes([r]))
        else:
            seed = random_seed(SIG_SEED_SIZE)
            seeds = expand_seeds(seed, ROUND, SIG_SEED_SIZE)

        round_columns = [x for r in range(ROUND) for x in expand_columns(seeds[r])]

        atfs = acting_on_atfs(atfs, round_columns)

        atf_bytes = [i.to_bytes(4, byteorder='big') for i in atfs]

        hash_input = hash_message+b''.join(atf_bytes)

        signed_message = hashlib.shake_256(hash_input).digest(CHLG_SIZE)
        
        chg_c, chg_nc, chg_val = expand_challenge(signed_message, CHLG_SIZE)

        temp_columns = [0] * (N*N*K)
        for r in range(K):
            temp_columns[r*(N*N):(r+1)*N*N] = round_columns[chg_nc[r]*N*N:(chg_nc[r]+1)*N*N]
        for i in range(N*N):
            for r in range(K):
                round_columns[i*K+r] = temp_columns[r*(N*N)+i]
        print(chg_val)
        if C >= K:
            temp_columns = [x for r in chg_val for x in expand_columns(seeds_sk[r*MAT_SK_SEED_SIZE:(r+1)*MAT_SK_SEED_SIZE])]
        else:
            temp_columns = [0] * (N*N*K)
            expanded = [-1] * C
            for r in range(K):
                if expanded[chg_val[r]] != -1:
                    temp_columns[r*(N*N):(r+1)*N*N] = temp_columns[expanded[chg_val[r]]*(N*N):(expanded[chg_val[r]]+1)*N*N]
                else:
                    temp_columns[r*(N*N):(r+1)*N*N] = expand_columns(seeds_sk[chg_val[r]*MAT_SK_SEED_SIZE:(chg_val[r]+1)*MAT_SK_SEED_SIZE])
                    expanded[chg_val[r]] = r
        challenge_columns = [0] * (N*N*K)
        for i in range(N*N):
            for r in range(K):
                challenge_columns[i*K+r] = temp_columns[r*(N*N)+i]
        
        final_matrices = columns_matrix(round_columns, challenge_columns)
        

        success = columns_decomposition(final_matrices, K)

        c_seeds = [0] * (ROUND-K)
        i = 0
        for indices in chg_c:
            c_seeds[i] = seeds[indices]
            i += 1

        if USE_SALT:
            return (signed_message, salt, c_seeds, final_matrices)
        else :
            return (signed_message, c_seeds, final_matrices)

    


        

    



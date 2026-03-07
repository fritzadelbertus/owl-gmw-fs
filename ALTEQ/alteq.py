from params import (SK_SEED_SIZE, MAT_SK_SEED_SIZE,
    C, K, N, PRIME, ROUND, LEN, PK_SEED_SIZE, MSG_HASH_SIZE, USE_SALT, 
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
    vectorized_columns = [0] * (N*N*C)
    for i in range(N*N):
        for r in range(C):
            vectorized_columns[i*C + r] = columns[r*N*N+i]

    atfs = expand_atfs(seeds[C*MAT_SK_SEED_SIZE:], C, C)
    atfs = inverting_on_atf(atfs, vectorized_columns)
    
    public_key = (atfs,seeds[C*MAT_SK_SEED_SIZE:])
    return (public_key, secret_key)

def alteq_sign(message:list[bytes], secret_key: bytes):
    if len(message) <= 0:
        raise Exception("can't sign empty message")
        return
    
    seeds_sk = expand_seeds(secret_key, 1, (MAT_SK_SEED_SIZE*C) + PK_SEED_SIZE)[0]

    

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
        vectorized_round_columns = [0] * (N*N*ROUND)
        for i in range(N*N):
            for r in range(ROUND):
                vectorized_round_columns[i*ROUND + r] = round_columns[r*N*N+i]

        atfs = expand_atfs(seeds_sk[C*MAT_SK_SEED_SIZE:], ROUND, ROUND)

        atfs = acting_on_atfs(atfs, vectorized_round_columns)

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
        round_columns = round_columns[:N*N*K]

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
        # print(challenge_columns)
        # print(round_columns)
        final_matrices = columns_matrix(round_columns, challenge_columns)

        success, final_matrices = columns_decomposition(final_matrices, K)

    c_seeds = [0] * (ROUND-K)
    for r in range(ROUND-K):
        c_seeds[r] = seeds[chg_c[r]]

    return (signed_message, c_seeds, final_matrices, atfs)

def alteq_verify(message, public_key, signature):
    signed_message, c_seeds, final_matrices, real_atfs = signature

    chg_c, chg_nc, chg_val = expand_challenge(signed_message, CHLG_SIZE)
    print(chg_nc)
    print(chg_val)
    columns = [0] * (N*N*ROUND)
    c_columns = [x for r in c_seeds for x in expand_columns(r)]

    # for i in range(N*N):
    #     for r in range(ROUND-K):
    #         columns[i*ROUND + chg_c[r]] = c_columns[r*N*N+i]
    
    # for i in range(N*N):
    #     for r in range(K):
    #         columns[i*ROUND + chg_nc[r]] = final_matrices[i*K+r]
    
    # # Forgeries Checking
    # for i in range(N*N):
    #     for r in range(ROUND-K, ROUND):
    #         if columns[i*ROUND+r] >= PRIME:
    #             return 1
    # for i in range(N):
    #     for r in range(ROUND-K, ROUND):
    #         if columns[i*(N+1)*ROUND+r] == 0:
    #             return 1
    # print("hello")
    # public_atfs, pk_seed = public_key
    
    # atfs = expand_atfs(pk_seed, ROUND, ROUND)
    # # print(public_atfs)
    # for i in range(LEN):
    #     for r in range(K):
    #         atfs[i*ROUND+chg_nc[r]] = public_atfs[i*C+chg_val[r]]
    
    # atfs = acting_on_atfs(atfs, columns)
    # count = 0
    # rcount = 0
    # list_reals = []
    # list_weirds = []
    # for i in range(len(atfs)):
    #     if real_atfs[i] != atfs[i]:
    #         count += 1
    #         list_reals.append(real_atfs[i])
    #         list_weirds.append(atfs[i])
    #         if i%ROUND in chg_nc:
    #             rcount += 1
    # print(count)

    # hash_message = hashlib.shake_256(message).digest(MSG_HASH_SIZE)

    
    # atf_bytes = [i.to_bytes(4, byteorder='big') for i in atfs]

    # hash_input = hash_message+b''.join(atf_bytes)

#==================================================

    for i in range(N*N):
        for r in range(ROUND-K):
            columns[i*ROUND + r] = c_columns[r*N*N+i]
    
    for i in range(N*N):
        for r in range(K):
            columns[i*ROUND + ROUND-K + r] = final_matrices[i*K+r]
    
    # Forgeries Checking
    for i in range(N*N):
        for r in range(ROUND-K, ROUND):
            if columns[i*ROUND+r] >= PRIME:
                return 1
    for i in range(N):
        for r in range(ROUND-K, ROUND):
            if columns[i*(N+1)*ROUND+r] == 0:
                return 1
    print("hello")
    public_atfs, pk_seed = public_key

    atfs = expand_atfs(pk_seed, ROUND-K, ROUND)
    for i in range(LEN):
        for r in range(K):
            atfs[i*ROUND+ROUND-K+r] = public_atfs[LEN*chg_val[r]+i]
    
    atfs = acting_on_atfs(atfs, columns)

    final_atfs = [0] * (LEN*ROUND)

    for r in range(ROUND-K):
        final_atfs[chg_c[r]*LEN:(chg_c[r]+1)*LEN] = atfs[r*LEN:(r+1)*LEN]
    for r in range(ROUND-K, ROUND):
        final_atfs[(chg_nc[r - ROUND + K])*LEN:(chg_nc[r - ROUND + K]+1)*LEN] = atfs[r*LEN:(r+1)*LEN]

    count = 0
    for i in range(len(atfs)):
        if real_atfs[i] != atfs[i]:
            count += 1
    print(count)
    hash_message = hashlib.shake_256(message).digest(MSG_HASH_SIZE)

    atf_bytes = [i.to_bytes(4, byteorder='big') for i in final_atfs]

    hash_input = hash_message+b''.join(atf_bytes)
    chk = hashlib.shake_256(hash_input).digest(CHLG_SIZE)
    print(chk)
    print(signed_message)
    return 0 if chk == signed_message else 1
        

        

    



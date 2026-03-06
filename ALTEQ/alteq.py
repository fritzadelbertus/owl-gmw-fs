from params import (SK_SEED_SIZE, MAT_SK_SEED_SIZE,
    C, PK_SEED_SIZE)

from expand import (random_seed, expand_seeds, 
    expand_columns, expand_atfs)

from atf import inverting_on_atf

def alteq_keygen():

    secret_key = random_seed(SK_SEED_SIZE)
    # dummy key
    # secret_key = b'\xa1\xf3\x8c\x9d\x02\xaf\x11\x5e\xde\x7a\x88\xbb\xef\x01\x9d\x30\x4c\x7e\x12\xab\x56\x9f\xcd\xef\x02\x11\x33\x45\x67\x89\xab\xcd\xef'

    seeds = expand_seeds(secret_key, 1, (MAT_SK_SEED_SIZE*C) + PK_SEED_SIZE)[0]

    columns = [x for r in range(C) for x in expand_columns(seeds[r*MAT_SK_SEED_SIZE:(r+1)*MAT_SK_SEED_SIZE])]
    
    atfs = expand_atfs(seeds[C*MAT_SK_SEED_SIZE:], C, C)
    atfs = inverting_on_atf(atfs, columns)
    
    public_key = (atfs,seeds[C*MAT_SK_SEED_SIZE:])
    return(public_key, secret_key)


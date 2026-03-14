from random_rng import shake_rng, deterministic_sample, get_random_value

def expand_challenge(seed, seed_size, ROUND, C):
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
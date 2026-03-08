from params import PRIME, LOG_Q, N, C

def multiplication_mod_p(a: int, b: int) -> int:
    r = a * b
    r = r - PRIME * (r >> LOG_Q)
    return r

def reduction_mod_p(a: int) -> int:
    r = a
    r -= PRIME * (r >> 32)
    r -= PRIME * (r >> 32)
    return r & 0xFFFFFFFF

def reduction_strict(x: int) -> int:
    if x >= PRIME:
        return x - PRIME
    return x & 0xFFFFFFFF

def inversion_modulo_p(a: int) -> int:
    return pow(a, PRIME - 2, PRIME)

def set_inversion_mod_p(set_list):
    mul = [0] * (N * C)
    inv0 = [0] * C

    for r in range(C):
        mul[r] = set_list[r]

    for i in range(1, N):
        for r in range(C):
            mul[i*C+r] = reduction_mod_p(
                multiplication_mod_p(
                    mul[(i-1)*C+r],
                    set_list[i*C+r]
                )
            )

    for r in range(C):
        inv0[r] = inversion_modulo_p(mul[(N-1)*C+r])

    for i in range(N-1, 0, -1):
        for r in range(C):
            tmp = reduction_mod_p(
                multiplication_mod_p(
                    mul[(i-1)*C+r],
                    inv0[r]
                )
            )
            inv0[r] = reduction_mod_p(
                multiplication_mod_p(
                    inv0[r],
                    set_list[i*C+r]
                )
            )
            set_list[i*C+r] = tmp

    for r in range(C):
        set_list[r] = inv0[r]

    return set_list
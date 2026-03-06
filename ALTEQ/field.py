from params import PRIME, LOG_Q

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

def set_inversion_mod_p(set_list: list[int]) -> list[int]:
    n = len(set_list)

    mul = [0] * n
    mul[0] = set_list[0]
    for i in range(1, n):
        mul[i] = reduction_mod_p(multiplication_mod_p(mul[i-1], set_list[i]))

    inv0 = inversion_modulo_p(mul[n-1])

    for i in reversed(range(1, n)):
        set_list[i] = reduction_mod_p(multiplication_mod_p(mul[i-1], inv0))
        inv0 = reduction_mod_p(multiplication_mod_p(inv0, set_list[i]))

    set_list[0] = inv0

    return set_list

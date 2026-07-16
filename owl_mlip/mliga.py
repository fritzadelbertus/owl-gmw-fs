# mLIGA
# A group action based on the module Lattice Isomorphism Problem


from cga import Cga, Group, Set
from owl_mlip.codec import decode_group, decode_set, encode_group, encode_set
from owl_mlip.params import LOGN
from owl_mlip.params import NTT_MODULUS as Q
from owl_mlip.poly import adjoint, poly_mul_ntt
from owl_mlip.sample import sample_group, sample_set
from owl_mlip.typing import GroupElementT, SetElementT


def set_equal(x: SetElementT, y: SetElementT) -> bool:
    return x == y


set_element_length = 2 * 2 * (2**LOGN) * 4

zero_entry = [0 for _ in range(2**LOGN)]
set_placeholder = [[zero_entry, zero_entry], [zero_entry, zero_entry]]

# Initialize the Set
mliga_set = Set(
    sample_set, encode_set, set_element_length, decode_set, set_equal, set_placeholder
)


# Initialize the Group
def poly_add_q(p0: list[int], p1: list[int]) -> list[int]:
    return [(e0 + e1) % Q for (e0, e1) in zip(p0, p1, strict=True)]


def matrix_mult(a: GroupElementT, b: GroupElementT) -> GroupElementT:
    c00 = poly_add_q(
        poly_mul_ntt(a[0][0], b[0][0], Q), poly_mul_ntt(a[0][1], b[1][0], Q)
    )
    c01 = poly_add_q(
        poly_mul_ntt(a[0][0], b[0][1], Q), poly_mul_ntt(a[0][1], b[1][1], Q)
    )
    c10 = poly_add_q(
        poly_mul_ntt(a[1][0], b[0][0], Q), poly_mul_ntt(a[1][1], b[1][0], Q)
    )
    c11 = poly_add_q(
        poly_mul_ntt(a[1][0], b[0][1], Q), poly_mul_ntt(a[1][1], b[1][1], Q)
    )
    return [[c00, c01], [c10, c11]]


def group_inverse(g: GroupElementT) -> GroupElementT:
    a = g[0][0]
    minus_b = [(i * (-1)) % Q for i in g[0][1]]
    minus_c = [(i * (-1)) % Q for i in g[1][0]]
    d = g[1][1]
    return [[d, minus_b], [minus_c, a]]


def operate_group(a: GroupElementT, b: GroupElementT) -> GroupElementT:
    return matrix_mult(b, a)


def group_equal(x: GroupElementT, y: GroupElementT) -> bool:
    return x == y


group_element_length = 2 * 2 * (2**LOGN) * 4
group_placeholder = [[zero_entry, zero_entry], [zero_entry, zero_entry]]


mliga_group = Group(
    sample_group,
    operate_group,
    group_inverse,
    group_element_length,
    encode_group,
    decode_group,
    group_equal,
    group_placeholder,
)

# Initialize the Group Action


def mliga_action(q: SetElementT, b: GroupElementT) -> SetElementT:
    b_star = [
        [adjoint(b[0][0]), adjoint(b[1][0])],
        [adjoint(b[0][1]), adjoint(b[1][1])],
    ]
    return matrix_mult(matrix_mult(b_star, q), b)


mliga = Cga(mliga_set, mliga_group, mliga_action)

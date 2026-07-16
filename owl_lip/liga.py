# LIGA
# A group action based on the Lattice Isomorphism Problem

import numpy as np

from cga import Cga, Group, Set
from owl_lip.codec import decode_group, decode_set, encode_group, encode_set
from owl_lip.field import matmul_fp
from owl_lip.params import INT_NBYTES, MATRIX_ENTRY_TYPE, P
from owl_lip.sample import sample_group, sample_set
from owl_lip.typing import GroupElementT, SetElementT


def set_equal(x: SetElementT, y: SetElementT) -> bool:
    return np.array_equal(x, y)


set_element_length = 4 * INT_NBYTES

set_placeholder = np.zeros((2, 2), dtype=MATRIX_ENTRY_TYPE)

# Initialize the Set
liga_set = Set(
    sample_set, encode_set, set_element_length, decode_set, set_equal, set_placeholder
)


def group_equal(x: GroupElementT, y: GroupElementT) -> bool:
    return np.array_equal(x, y)


# Initialize the Group
def group_inverse(g: GroupElementT) -> GroupElementT:
    if g.shape != (2, 2):
        raise ValueError("Expected 2x2 matrix")

    a = int(g[0, 0])
    b = int(g[0, 1])
    c = int(g[1, 0])
    d = int(g[1, 1])

    det = (a * d - b * c) % P

    if det not in (1, P - 1):
        raise ValueError(f"Matrix is not unimodular, det={det}")

    det_inv = pow(det, -1, P)

    return np.array(
        [
            [(det_inv * d) % P, (det_inv * (-b)) % P],
            [(det_inv * (-c)) % P, (det_inv * a) % P],
        ],
        dtype=MATRIX_ENTRY_TYPE,
    )

def operate_group(a: GroupElementT, b: GroupElementT) -> GroupElementT:
    return matmul_fp(b, a)


group_element_length = 4 * INT_NBYTES

group_placeholder = np.zeros((2, 2), dtype=MATRIX_ENTRY_TYPE)

liga_group = Group(
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
def liga_action(q: SetElementT, b: GroupElementT) -> GroupElementT:
    return matmul_fp(matmul_fp(b.T, q), b)


liga = Cga(liga_set, liga_group, liga_action)

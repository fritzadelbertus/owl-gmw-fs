from cga import Cga, Group, Set
from custom_gmwfs.codec import decode_group, decode_set, encode_group, encode_set
from custom_gmwfs.sample import sample_group, sample_set
from custom_gmwfs.typing import GroupElementT, SetElementT


def set_equal(s1: SetElementT, s2: SetElementT) -> bool:
    return True


set_element_length = 0

set_placeholder = 0

# Initialize the Set
_set = Set(
    sample_set, encode_set, set_element_length, decode_set, set_equal, set_placeholder
)


def group_equal(g1: GroupElementT, g2: GroupElementT) -> bool:
    return True


# Initialize the Group
def group_inverse(g: GroupElementT) -> GroupElementT:
    return g


def operate_group(g1: GroupElementT, g2: GroupElementT) -> GroupElementT:
    return g1


group_element_length = 0

group_placeholder = 0

_group = Group(
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
def action(s: SetElementT, g: GroupElementT) -> SetElementT:
    return s


group_action = Cga(_set, _group, action)

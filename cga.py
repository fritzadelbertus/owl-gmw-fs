# Cryptographic Group Action
from rngcontext import RngContext
from type_defs import (
    Decoder,
    Encoder,
    Equality,
    GroupAction,
    GroupElementT,
    GroupInverse,
    GroupOperator,
    Sample,
    SetElementT,
)


class Set[SetElementT]:
    def __init__(
        self,
        sample: Sample[SetElementT],
        set_encoder: Encoder[SetElementT],
        set_element_length: int,
        set_decoder: Decoder[SetElementT],
        equal_element: Equality[SetElementT],
        placeholder: SetElementT,
    ):
        self._sample = sample
        self._element_length = set_element_length
        self._encoder = set_encoder
        self._decoder = set_decoder
        self._equal = equal_element
        self._placeholder = placeholder

    def sample(self, rng: RngContext) -> SetElementT:
        return self._sample(rng)

    def encode(self, s: SetElementT) -> bytes:
        return self._encoder(s)

    def decode(self, bs: bytes) -> SetElementT:
        return self._decoder(bs)

    def element_length(self) -> int:
        return self._element_length

    def equal(self, s1: SetElementT, s2: SetElementT) -> bool:
        return self._equal(s1, s2)

    def placeholder(self):
        return self._placeholder


class Group[GroupElementT]:
    def __init__(
        self,
        sample: Sample[GroupElementT],
        operator: GroupOperator[GroupElementT],
        group_inverse: GroupInverse[GroupElementT],
        group_element_length: int,
        group_encoder: Encoder[GroupElementT],
        group_decoder: Decoder[GroupElementT],
        equal_element: Equality[GroupElementT],
        placeholder: GroupElementT,
    ):
        self._sample = sample
        self._operator = operator
        self._inverse = group_inverse
        self._encoder = group_encoder
        self._decoder = group_decoder
        self._element_length = group_element_length
        self._equal = equal_element
        self._placeholder = placeholder

    def sample(self, rng: RngContext) -> GroupElementT:
        return self._sample(rng)

    def operate(self, g1: GroupElementT, g2: GroupElementT) -> GroupElementT:
        return self._operator(g1, g2)

    def inverse(self, g: GroupElementT) -> GroupElementT:
        return self._inverse(g)

    def element_length(self) -> int:
        return self._element_length

    def encode(self, g: GroupElementT) -> bytes:
        return self._encoder(g)

    def decode(self, bs: bytes) -> GroupElementT:
        return self._decoder(bs)

    def equal(self, g1: GroupElementT, g2: GroupElementT) -> bool:
        return self._equal(g1, g2)

    def placeholder(self):
        return self._placeholder


class Cga[GroupElementT, SetElementT]:
    def __init__(
        self,
        cga_set: Set[SetElementT],
        cga_group: Group[GroupElementT],
        cga_action: GroupAction[SetElementT, GroupElementT],
    ):
        self.set = cga_set
        self.group = cga_group
        self._action = cga_action

    def action(self, s: SetElementT, g: GroupElementT) -> SetElementT:
        return self._action(s, g)

from collections.abc import Callable
from typing import TypeVar

from rngcontext import RngContext

GroupElementT = TypeVar("GroupElementT")
SetElementT = TypeVar("SetElementT")

PublicKeyT = TypeVar("PublicKeyT")
SecretKeyT = TypeVar("SecretKeyT")
SignatureT = TypeVar("SignatureT")
ContextT = TypeVar("ContextT")

type GmwfsSecretKeyT[GroupElementT] = list[GroupElementT]
type GmwfsPublicKeyT[SetElementT] = tuple[SetElementT, list[SetElementT]]
type GmwfsSignatureT[GroupElementT] = tuple[bytes, list[GroupElementT]]
type GmwfsContextT = None

# cga.py
GroupOperator = Callable[[GroupElementT, GroupElementT], GroupElementT]
GroupInverse = Callable[[GroupElementT], GroupElementT]
GroupAction = Callable[[SetElementT, GroupElementT], SetElementT]
GroupEquality = Callable[[GroupElementT, GroupElementT], bool]
SetEquality = Callable[[SetElementT, SetElementT], bool]


ValueT = TypeVar("ValueT")
Sample = Callable[[RngContext], ValueT]
Equality = Callable[[ValueT, ValueT], bool]
Encoder = Callable[[ValueT], bytes]
Decoder = Callable[[bytes], ValueT]

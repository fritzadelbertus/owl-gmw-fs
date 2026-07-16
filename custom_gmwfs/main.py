import hashlib

from custom_gmwfs.group_action import group_action
from custom_gmwfs.params import CHG_SIZE, ORBIT, C, K, R
from expand import expand_challenge
from gmwfs import Gmwfs


def expand(chg: bytes) -> tuple[list[int], list[int], list[int]]:
    return expand_challenge(chg, C, K, R, CHG_SIZE)


def hashfunc(input: bytes) -> bytes:
    return hashlib.shake_256(input).digest(CHG_SIZE)


gmwfs = Gmwfs("scheme-name", group_action, C, K, R, hashfunc, expand, CHG_SIZE, ORBIT)

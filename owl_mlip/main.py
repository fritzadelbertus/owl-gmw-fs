import hashlib

from expand import expand_challenge
from gmwfs import Gmwfs
from owl_mlip.mliga import mliga
from owl_mlip.params import CHG_SIZE, ORBIT, C, K, R


def expand(chg: bytes) -> tuple[list[int], list[int], list[int]]:
    return expand_challenge(chg, C, K, R, CHG_SIZE)


def hashfunc(input: bytes) -> bytes:
    return hashlib.shake_256(input).digest(CHG_SIZE)


owl_mlip = Gmwfs("owl-mLIP", mliga, C, K, R, hashfunc, expand, CHG_SIZE, ORBIT)

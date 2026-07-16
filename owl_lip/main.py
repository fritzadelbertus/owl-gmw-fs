import hashlib

from expand import expand_challenge
from gmwfs import Gmwfs
from owl_lip.liga import liga
from owl_lip.params import CHG_SIZE, ORBIT, C, K, R


def expand(chg):
    return expand_challenge(chg, C, K, R, CHG_SIZE)


def hashfunc(input):
    return hashlib.shake_256(input).digest(CHG_SIZE)


owl_lip = Gmwfs("owl-LIP", liga, C, K, R, hashfunc, expand, CHG_SIZE, ORBIT)

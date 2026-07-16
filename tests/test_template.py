import os

import numpy as np
import pytest

# == Modify Only this part ==
from custom_gmwfs.group_action import group_action as cga
from custom_gmwfs.main import gmwfs as gmwfs
from rngcontext import RngContext as RngContext

num_of_tests = 5
# == Modify Only this part ==
R = gmwfs.r
K = gmwfs.k
chg_size = gmwfs.chg_size
hash = gmwfs.hash
expand = gmwfs.expand
keygen = gmwfs.keygen
group_equal = cga.group.equal
set_equal = cga.set.equal
sample_set = cga.set.sample
sample_group = cga.group.sample
encode_set = cga.set.encode
encode_group = cga.group.encode
decode_set = cga.set.decode
decode_group = cga.group.decode
group_operator = cga.group.operate
group_inverse = cga.group.inverse
group_placeholder = cga.group.placeholder()
set_placeholder = cga.set.placeholder()
group_action = cga.action


SEEDS = range(num_of_tests)


def make_rng(seed):
    return RngContext(np.array([seed], dtype=np.uint8))


@pytest.mark.parametrize("seed", SEEDS, ids=lambda s: f"seed={s}")
def test_sample_set_is_deterministic(seed):
    x1 = sample_set(make_rng(seed))
    x2 = sample_set(make_rng(seed))
    assert set_equal(x1, x2)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda s: f"seed={s}")
def test_sample_group_is_deterministic(seed):
    g1 = sample_group(make_rng(seed))
    g2 = sample_group(make_rng(seed))
    assert group_equal(g1, g2)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda s: f"seed={s}")
def test_set_codec(seed):
    x = sample_set(make_rng(seed))
    encoded = encode_set(x)

    assert isinstance(encoded, bytes)
    assert set_equal(decode_set(encoded), x)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda s: f"seed={s}")
def test_group_codec(seed):
    g = sample_group(make_rng(seed))
    encoded = encode_group(g)

    assert isinstance(encoded, bytes)
    assert group_equal(decode_group(encoded), g)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda s: f"seed={s}")
def test_group_operator_is_associative(seed):
    rng = make_rng(seed)
    g1 = sample_group(rng)
    g2 = sample_group(rng)
    g3 = sample_group(rng)

    left = group_operator(group_operator(g1, g2), g3)
    right = group_operator(g1, group_operator(g2, g3))

    assert group_equal(left, right)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda s: f"seed={s}")
def test_group_inverse(seed):
    g = sample_group(make_rng(seed))
    g_inv = group_inverse(g)

    left_identity = group_operator(g_inv, g)
    right_identity = group_operator(g, g_inv)

    assert group_equal(left_identity, right_identity)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda s: f"seed={s}")
def test_group_action_compatibility(seed):
    rng = make_rng(seed)
    x = sample_set(rng)
    g1 = sample_group(rng)
    g2 = sample_group(rng)

    combined = group_action(x, group_operator(g1, g2))
    sequential = group_action(group_action(x, g2), g1)

    assert set_equal(combined, sequential)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda s: f"seed={s}")
def test_identity_action(seed):
    rng = make_rng(seed)
    x = sample_set(rng)
    g = sample_group(rng)

    identity = group_operator(g, group_inverse(g))

    assert set_equal(x, group_action(x, identity))


@pytest.mark.parametrize("seed", SEEDS, ids=lambda s: f"seed={s}")
def test_sign_and_verify_element(seed):
    message = os.urandom(32)
    rng = make_rng(seed)
    public_key, secret_key, _ = keygen()
    orbit, public = public_key
    h_i = [sample_group(rng) for i in range(R)]
    t_i = [group_action(orbit, h_i[i]) for i in range(R)]

    hash_input = message
    for s in t_i:
        hash_input += encode_set(s)

    chg = hash(hash_input)
    chg_c, chg_nc, chg_val = expand(chg)

    f_i = [group_placeholder] * R
    for k in range(K):
        f_i[chg_nc[k]] = group_operator(h_i[chg_nc[k]], secret_key[chg_val[k]])

    for k in range(R - K):
        f_i[chg_c[k]] = h_i[chg_c[k]]

    ver_t_i = [set_placeholder] * R
    for k in range(K):
        ver_t_i[chg_nc[k]] = group_action(public[chg_val[k]], f_i[chg_nc[k]])

    for k in range(R - K):
        ver_t_i[chg_c[k]] = group_action(orbit, f_i[chg_c[k]])

    valid = True
    for i in range(len(t_i)):
        if not set_equal(t_i[i], ver_t_i[i]):
            valid = False
            break

    assert valid

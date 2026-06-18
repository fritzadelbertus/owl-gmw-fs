from OWL_LIP.params import L

import random
import numpy as np


def alg_add(x, y):
    return tuple(x[i] + y[i] for i in range(4))

def alg_mul_int(k, x):
    return tuple(k * x[i] for i in range(4))

def alg_zero():
    return (0, 0, 0, 0)



def canonical_pm(B):
    B = np.array(B, dtype=object)
    flat = B.flatten()
    for x in flat:
        if x != 0:
            return B if x > 0 else -B
    return B

def sample_k():
    return random.choice([-2, -1, 1, 2])

Q = [
    [(1,1,0,0), (1,0,0,0)],
    [(1,0,0,0), (1,0,1,0)]
]

# ==============================================================================================================
# The Group Action

def mat_alg_action(B, Q):
    B = np.array(B, dtype=object)

    R = [[alg_zero(), alg_zero()],
         [alg_zero(), alg_zero()]]

    for i in range(2):
        for j in range(2):
            s = alg_zero()
            for p in range(2):
                for q in range(2):
                    term = alg_mul_int(B[p, i] * B[q, j], Q[p][q])
                    s = alg_add(s, term)
            R[i][j] = s

    return R

def action(B,Q):
    return mat_alg_action(B,Q)

# ==============================================================================================================
# Group and Set Sampling
# Sample a random basis matrix
def group_sampler(rng=None): # Samples a random matrix
    B = np.eye(2, dtype=object)

    for _ in range(L):
        move = random.choice(["E12", "E21", "S"])
        if move == "E12":
            k = sample_k()
            E = np.array([[1, k], [0, 1]], dtype=object)
        elif move == "E21":
            k = sample_k()
            E = np.array([[1, 0], [k, 1]], dtype=object)
        else:
            E = np.array([[0, 1], [1, 0]], dtype=object)

        B = B @ E

    return canonical_pm(B)

# Sample a random gram matrix
def set_sampler(rng=None):
    Q = [
        [(1,1,0,0), (1,0,0,0)],
        [(1,0,0,0), (1,0,1,0)]
    ]
    B = group_sampler()
    return mat_alg_action(B,Q)

# ==============================================================================================================
# Group Operations

# Calculate group inverse
def group_inverse(A):
    A = np.array(A, dtype=object)

    if A.shape != (2, 2):
        raise ValueError("Expected 2x2 matrix")

    a, b = A[0, 0], A[0, 1]
    c, d = A[1, 0], A[1, 1]

    det = a*d - b*c

    if det not in (1, -1):
        raise ValueError(f"Matrix is not unimodular, det={det}")

    return np.array([
        [ d // det, -b // det],
        [-c // det,  a // det]
    ], dtype=object)

# Calculate group operation (matrix multiplication)
def group_operator(A, B):
    return np.array(B, dtype=object) @ np.array(A, dtype=object)


group_identity = np.eye(2, dtype=object)
# This section contains OWL Parameters
import numpy as np

LAMBDA = 128
L = 20
C = 7
K = 22
R = 84
CHG_SIZE = LAMBDA // 4
INT_NBYTES = 4
P = 2**32 - 5

MATRIX_ENTRY_TYPE = np.uint32

ORBIT = (
    np.array(
        [[2414213562, 1000000000], [1000000000, 2732050808]],
        dtype=MATRIX_ENTRY_TYPE,
    )
    % P
)

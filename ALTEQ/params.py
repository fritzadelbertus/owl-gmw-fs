LAMBDA = 128
N = 13
LOG_Q = 32
PRIME = 4294967291
C = 7
K = 22
ROUND = 84

LEN = N * (N - 1) * (N - 2) // 6
ALT_SIZE = LEN * LOG_Q // 8 # ATF size (public)
MAT_SIZE = N * N * LOG_Q // 8 # Matrix size (secret)



MAT_SK_SEED_SIZE = 32
SK_SEED_SIZE = LAMBDA//4
PK_SEED_SIZE = LAMBDA//4
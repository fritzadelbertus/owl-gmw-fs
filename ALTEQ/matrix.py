from params import N,K, PRIME
from field import reduction_mod_p, multiplication_mod_p, inversion_modulo_p
def column_mul(mat, column, j, s):

    for i in range(N):
        if (i != j):
            for k in range(s,N):
                for r in range(K):
                    mat[(i*N+k)*K+r] = reduction_mod_p(
                        (mat[(i*N+k)*K+r]& 0xFFFFFFFFFFFFFFFF)+multiplication_mod_p(
                            mat[(j*N+k)*K+r], column[i*K+r]
                        ))
    
    for k in range(s,N):
        for r in range(K):
            mat[(j*N+k)*K+r] = reduction_mod_p(
                multiplication_mod_p(
                    mat[(j*N+k)*K+r], column[(j)*K+r]
                )
            )
        
    return mat

def column_inv(col, j):
    for r in range(K):
        col[(j)*K+r] = inversion_modulo_p(col[(j)*K+r])

    for i in range(N):
        if (i!=j):
            for r in range(K):
                col[(i)*K+r]=reduction_mod_p(
                    multiplication_mod_p(
                        (PRIME-col[(j)*K+r]),
                        col[(i)*K+r]
                    )
                )
    return col

def columns_matrix(colsA, colsB):
    mat = [0]*N*N*K

    for i in range(N):
        for j in range(N):
            for r in range(K):
                mat[(i*N+j)*K+r]=colsA[(j*N+i)*K+r]

    for j in range(N-1, -1, -1):
        mat = column_mul(mat, colsA[j*N*K:(j+1)*N*K], j, j+1)
    for j in range(N-1, -1, -1):
        mat = column_mul(mat, colsB[j*N*K:(j+1)*N*K], j, 0)

    return mat

def columns_decomposition(columns, nb_mats):
    mat0 = columns.copy()

    for i in range(N):
        columns[i*K:(i+1)*K] = mat0[i*N*K:(i*N+1)*K]

    cols = columns[:N*K]

    for r in range(K):
        if (cols[r]==0 or cols[r]==PRIME):
            return 0, []
    
    for j in range(1,N):
        cols = column_inv(cols, j-1)
        mat0 = column_mul(mat0, cols, j-1, j)

        for i in range(N):
            columns[(j*N+i)*K:(j*N+i+1)*K] = mat0[(i*N+j)*K:(i*N+j+1)*K]
        
        cols = columns[j*N*K:(j+1)*N*K]

        for r in range(K):
            if (cols[j*K+r]==0 or cols[j*K+r]==PRIME):
                return 0, []

    return 1, columns



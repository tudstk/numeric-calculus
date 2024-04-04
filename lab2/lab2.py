import copy
import numpy as np

EPS = 10e-15


def LU_decomposition_2(A):
    n = len(A)
    for p in range(n):
        for i in range(p, n):
            sum_LU = 0.0
            for k in range(p):
                sum_LU += A[i][k] * A[k][p]
            A[i][p] -= sum_LU
        for i in range(p + 1, n):
            sum_LU = 0.0
            for k in range(p):
                sum_LU += A[p][k] * A[k][i]
            A[p][i] = (A[p][i] - sum_LU) / A[p][p]
    return A


def determinant_tri(A):
    d = 1
    for i in range(0, len(A)):
        d *= A[i][i]
    return d


def euclidean_norm(vector):
    return np.sqrt(np.sum(vector ** 2))


def solve_LU(A, b):
    n = len(A)
    x = np.zeros(n)
    y = np.zeros(n)

    # substitutia directa Ly = b
    for i in range(n):
        sum_ly = 0.0
        for j in range(i):
            sum_ly += A[i][j] * y[j]
        y[i] = (b[i] - sum_ly) / A[i][i]

    # substitutia inversa Ux = y
    for i in range(n - 1, -1, -1):
        sum_ux = 0.0
        for j in range(i + 1, n):
            sum_ux += A[i][j] * x[j]
        x[i] = (y[i] - sum_ux)

    return x


A = np.array([[2.5, 2, 2], [5, 6, 5], [5, 6, 6.5]])
b = np.array([2, 2, 2])
b_init = copy.deepcopy(b)
A_init = copy.deepcopy(A)

print(LU_decomposition_2(A))

x = solve_LU(A, b)
assert abs(determinant_tri(A) - np.linalg.det(A_init)) < EPS
print(determinant_tri(A))

# norma #
Z = np.dot(A_init, x)
Z -= b
print(euclidean_norm(Z))

print("solutia *ESTIMATIVA* domnule:", x)

xlib = np.linalg.solve(A_init, b_init)
print("solutia calculata de LIBRARIE domnule:", xlib)
inverse_Alib = np.linalg.inv(A_init)
print("norma dintre X si Xlib:", euclidean_norm(x - xlib))
print("norma dintre X si restul...:", euclidean_norm(x - np.dot(inverse_Alib, b_init)))

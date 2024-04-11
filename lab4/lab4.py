import numpy as np
from scipy.sparse import lil_matrix, csr_matrix


def read_A_lil(filename):
    with open(filename, 'r') as file:
        n = int(file.readline())
        A = lil_matrix((n, n), dtype=np.float64)
        for line in file:
            elements = line.split(',')
            i = int(elements[1])
            j = int(elements[2])
            value = float(elements[0])
            A[i, j] = value
    return A


def read_A_csr(filename):
    A = read_A_lil(filename)
    return csr_matrix(A)


def read_B(filename):
    with open(filename, 'r') as file:
        n = int(file.readline())
        b = np.zeros(n)
        for i, line in enumerate(file):
            line = line.strip()
            if line:
                try:
                    b[i] = float(line)
                except ValueError:
                    print("linie ignorata:  ", line)
    return b


def gauss_seidel(A, b, epsilon=1e-8, max_iter=10000):
    n = len(b)
    x = np.zeros(n)
    k = 0
    while True:
        delta_x_max = 0
        for i in range(n):
            sum_ax = A[i].dot(x) - A[i, i] * x[i]
            delta_x = abs((b[i] - sum_ax) / A[i, i] - x[i])
            x[i] = (b[i] - sum_ax) / A[i, i]
            if delta_x > delta_x_max:
                delta_x_max = delta_x
        k += 1
        if delta_x_max < epsilon or k >= max_iter or delta_x_max > 1e8:
            break
    sol_aprox = A.dot(x)
    norma = np.linalg.norm(sol_aprox - b)

    if delta_x_max < epsilon:
        return x, k, norma
    else:
        print("Divergenta")

def read_data(i):
    A_lil = read_A_lil(f'resources/matrix_A/a_{i}.txt')
    A_csr = read_A_csr(f'resources/matrix_A/a_{i}.txt')
    b = read_B(f'resources/matrix_B/b_{i}.txt')
    return A_lil, A_csr, b


if __name__ == '__main__':
    for i in range(1, 6):
        A_lil, A_csr, b = read_data(i)
        diag_lil = np.all(A_lil.diagonal() != 0)
        diag_csr = np.all(A_csr.diagonal() != 0)
        print(f"elementele de pe diagonala matricei {i} pentru lil sunt nenule:", diag_lil)
        print(f"elementele de pe diagonala matricei {i} pentru csr sunt nenule:", diag_csr)
        x, iterations, norm_diff = gauss_seidel(A_lil, b)
        if x is not None:
            print(f"solutia pentru sistemul {i} este: {x}")
            print(f"nr de iteratii: {iterations}")
            print(f"norma: {norm_diff}\n")
import copy

import numpy as np

EPS = 10e-15


def householder(A, n, b):
    Q = np.zeros((n, n))
    for i in range(0, n):
        Q[i][i] = 1
    # print(Q)
    u = np.zeros(n)
    for r in range(0, n - 1):
        sigma = sum(A[i][r] ** 2 for i in range(r, n))
        if sigma <= EPS:
            break
        k = np.sqrt(sigma)
        if A[r][r] > 0:
            k *= -1
        beta = sigma - k * A[r][r]
        u[r] = A[r][r] - k
        for i in range(r + 1, n):
            u[i] = A[i][r]
        for j in range(r + 1, n):
            gamma = sum(u[i] * A[i][j] for i in range(r, n)) / beta
            for i in range(r, n):
                A[i][j] = A[i][j] - gamma * u[i]
        A[r][r] = k
        for i in range(r + 1, n):
            A[i][r] = 0
        gamma = sum(b[i] * u[i] for i in range(r, n)) / beta
        for i in range(r, n):
            b[i] = b[i] - gamma * u[i]
        for j in range(0, n):
            gamma = sum(u[i] * Q[i][j] for i in range(r, n)) / beta
            for i in range(r, n):
                Q[i][j] = Q[i][j] - gamma * u[i]
    return Q, A


def b_vector(s, A, n):
    b = np.zeros(n)
    for i in range(0, n):
        b[i] = sum(s[j] * A[i][j] for j in range(0, n))
    return b


def euclidean_norm(vector):
    return np.sqrt(np.sum(vector ** 2))


def solve(R, Q, b, n, inverse_flag=False):
    x = np.zeros(n)
    if not inverse_flag:
        target = np.dot(Q, b)
    else:
        target = b
    for i in range(n - 1, -1, -1):
        x[i] = (target[i] - sum(R[i][j] * x[j] for j in range(i + 1, n))) / R[i][i]

    return x


def inverse(A, Q, R, n):
    inv = np.zeros((n, n))
    if np.linalg.det(A) < EPS:
        print("Can't compute inverse.")
        return -1

    for j in range(0, n):
        b_inv = [row[j] for row in Q]
        x_star = solve(R, Q, b_inv, n, inverse_flag=True)
        inv[:, j] = x_star

    return inv


if __name__ == "__main__":
    A = [[0, 0, 4],
         [1, 2, 3],
         [0, 1, 2]]

    A_init = copy.deepcopy(A)
    n = np.ndim(A) + 1

    s = np.array([3, 2, 1])

    b = b_vector(s, A, n)
    b_init = copy.deepcopy(b)

    print("b:", b)

    Q_transpose, R = householder(A, n, b)
    print("Q:", Q_transpose)
    print("R:", R)
    x = solve(R, Q_transpose, b_init, n)
    print("Solutia data de catre NOI, domnule:", x)

    Qlib, Rlib = np.linalg.qr(A_init)
    xlib = solve(Rlib, np.transpose(Qlib), b_init, n)
    print("Solutia data de catre LIBRARIE, domnule:", xlib)

    norms = []
    norms.extend([euclidean_norm(xlib - x), euclidean_norm(np.dot(A_init, xlib) - b_init),
                  euclidean_norm(np.dot(A_init, xlib) - b_init), (euclidean_norm((x - s))) / euclidean_norm(s),
                  (euclidean_norm((xlib - s))) / euclidean_norm(s)])

    print("Prima norma:", norms[0])
    print("A doua norma:", norms[1])
    print("A doua norma:", norms[2])
    print("A patra norma:", norms[3])
    print("A cincea norma:", norms[4])

    for i in range(len(norms)):
        assert norms[i] < 10e-6

    A_inverse = inverse(A_init, Q_transpose, R, n)
    A_inverse_lib = np.linalg.inv(A_init)

    print("Ultima norma, domnul meu:", euclidean_norm(A_inverse - A_inverse_lib))

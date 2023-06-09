import numpy as np
import numpy.linalg as LA
import random as rand
import math

# Норма вектора а
def norma(a):
    return np.max(np.sum(np.absolute(a), axis=1))

# Реализация метода Якоби
def Jacobi(A, b, eps):
    #если матрица А обладает свойством диаганального преобладания, то можем строить последовательность векторов x_k
    for i in range(A.shape[0]):
        if (sum([A[i, j] for j in range(A.shape[1])]) - A[i, i] > A[i, i]):
            raise NameError("Якоби не сойдется")

    x = np.transpose(np.matrix([0.0 for i in b]))
    norm = 1
    it = 0
    while (norm > eps):
        it += 1
        x_t = np.copy(x)
        for i in range(len(A)):
            s = 0
            for j in range(len(A)):
                if (i != j): s += A[i, j] * x[j]
            x_t[i] = (b[i] - s) / A[i, i]
        norm = np.sqrt(sum((x_t[i] - x[i]) ** 2 for i in range(len(A)))) #для условия остановки
        x = np.copy(x_t)
    print('>>>> Posteriori', it, '\n')
    return x


# Реализация метода Зейделя
def Seidel(A, b, eps):
    n = len(A)
    x = np.transpose(np.matrix([0.0 for i in b]))

    norm = 1
    it = 0
    while (norm > eps):
        it += 1
        x_t = np.copy(x)
        for i in range(n):
            s1 = sum([A[i, j] * x_t[j] for j in range(i)])
            s2 = sum([A[i, j] * x[j] for j in range(i + 1, n)])
            x_t[i] = (b[i] - s1 - s2) / A[i, i]

        norm = np.sqrt(sum((x_t[i] - x[i]) ** 2 for i in range(n))) #для условия остановки
        x = x_t
    print('>>>> Posteriori', it, '\n')
    return x

# Точность вычислений, генерация матрицы и вектора
eps = 1e-7
size = rand.randint(3, 4)
b = np.transpose(np.matrix([rand.randint(1, 10) for i in range(size)]))
A = np.matrix([[rand.randint(0, 10) + 0.0 for j in range(size)] for i in range(size)])


# Создание нижнетреугольной, диагональной и верхнетреугольной матриц
for i in range(len(A)):
    A[i, i] = sum([A[i, j] + 5 for j in range(A.shape[1])])

L = np.zeros(A.shape)
D = np.zeros(A.shape)
R = np.zeros(A.shape)


for i in range(A.shape[0]):
    if i == 0:
        D[i, i] = A[i, i]
        R[i, i + 1:] = A[i, i + 1:]
    else:
        L[i, :i] = A[i, :i]
        D[i, i] = A[i, i]
        R[i, i + 1:] = A[i, i + 1:]

print("A=\n", A)
print("b=\n", b)
print('L \n', L, '\n')
print('D \n', D, '\n')
print('R \n', R, '\n')



print("\n                                     Jacobi         \n")
c = np.dot(LA.inv(D), b)
try:
    B = -np.dot(LA.inv(D), (L + R))
    x = Jacobi(A, b, eps * (1 - norma(B)) / norma(B))
    print('>>>> Priori', np.ceil(math.log(1e-7 * (1 - norma(B)) / norma(c), norma(B))), '\n')
    print("x (Jacobi) =\n", x)
    print("AX - b =\n", np.subtract(np.dot(A, x), b))
    print("\n||x^*-x|| = ", norma(LA.solve(A, b) - x))
    print("\n                                   Seidel           \n")
    LD = L + D
    q = norma(np.dot(LA.inv(LD), R))

    x = Seidel(A, b, eps * (1 - q) / q)
    print('>>>> Priori', np.ceil(math.log(1e-7 * abs(1 - q) / norma(c), q)), '\n')
    print("x (Seidel) =\n", x)
    print("AX - b =\n", np.subtract(np.dot(A, x), b))
    print("\n||x^*-x|| = ", norma(LA.solve(A, b) - x))
except BaseException as e:
    print(e)
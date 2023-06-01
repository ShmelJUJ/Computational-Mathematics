import random
import numpy

# Переменная для подсчёта кол-ва перестановок при LU-разложении
AMOUNT_PERMUTATIONS = 0

# Функция транспонирования матрицы
def transpose_matrix(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]

# Функция создания матрицы размера n со случайными элементами от 0 до 100
def generating_Matrix(n, matrix):
    for i in range(n):
        for j in range(n):
            number_random = random.randint(0, 100)
            matrix[i][j] = number_random
    print("A = \n", end='')
    print(numpy.matrix(matrix), '\n')


# Функция LUP-разложения матрицы
def LU_decomosition(n, A, P):
    global AMOUNT_PERMUTATIONS
    for i in range(len(A)):
        pivot_value = 0
        pivot = -1
        for row in range(i, len(A)):
            if abs(A[row][i]) > pivot_value:
                pivot_value = A[row][i]
                pivot = row
        if pivot_value != 0:
            AMOUNT_PERMUTATIONS+=1
            P[i], P[pivot] = P[pivot], P[i]
            A[i], A[pivot] = A[pivot], A[i]
            for j in range(i + 1, len(A)):
                A[j][i] /= A[i][i]
                for k in range(i + 1, len(A)):
                    A[j][k] -= A[j][i] * A[i][k]
    print("Матрица после LUP разложения", '\n', numpy.matrix(A), '\n')


# Функция решения уравнения Ax = b
def solving_equation(A, b):
    y = [0 for i in range(len(A))]
    # (Ly=Pb,L-верхняя треугольная)
    for i in range(len(y)):
        y[i] = b[i] - sum([A[i][k] * y[k] for k in range(0, i)])
    x = [0 for i in range(len(A))]
    # (Ux=y,U-нижняя треугольная)
    for i in range(len(x) - 1, -1, -1):
        x[i] = (y[i] - sum([A[i][k] * x[k] for k in range(i + 1, len(y))])) / A[i][i]
    return x


# Функция нахождения обратной матрицы
def inverse_matrix(A, P):
    A_inverse = numpy.eye(len(A))
    E = numpy.eye(len(A))
    PT=P.copy()
    for i in range(len(P)):
        PT[P[i]] = i
    for i in range(len(A)):
        A_inverse[:, i] = solving_equation(A, E[PT[i]])
    return A_inverse


# Функция вычисления нормы матрицы
def norma(matrix) -> float:
    max_sum = 0
    n = len(matrix)
    for i in range(n):
        temp_sum = 0
        for j in range(n):
            temp_sum += abs(matrix[j][i])
        if (temp_sum > max_sum):
            max_sum = temp_sum
    return max_sum


# Функция вычисления числа обусловленности матрицы
def conditional_number(A, A_inverse) -> float:
    number = norma(A) * norma(A_inverse)
    return number


# Функция вычисления определителя
def determinant(A):
    detA = 1
    for i in range(len(A)):
        detA *= A[i][i]
    detA = detA if AMOUNT_PERMUTATIONS % 2 == 0 else -detA
    return detA


def check_solution(A, x, b):
    print("Ax - b =", numpy.subtract(numpy.matmul(A, x), b))


"""def check_inversematrix(A, A_inverse):
    print("\n Левосторонее умножение: \n", numpy.matmul(A_inverse, A))
    print("\n Правосторонее умножение: \n", numpy.matmul(A, A_inverse))"""


def main():
    matrix_size = int(input("Enter matrix size: "))
    A = [[0] * matrix_size for i in range(matrix_size)]
    generating_Matrix(matrix_size, A)
    A_copy = numpy.copy(A)
    P = [i for i in range(matrix_size)]

    LU_decomosition(matrix_size, A, P)
    U=numpy.copy(A)
    L=numpy.eye(matrix_size)
    for i in range(1,matrix_size):
        for j in range(0,i):
            L[i,j]=A[i][j]
            U[i,j]=0
    print('L =\n',L,'\n \n U =\n',U,'\n \n LU= \n',numpy.matmul(L,U),'\nPA=\n',numpy.matmul([[1 if P[i]==j else 0 for j in range(matrix_size)] for i in range(matrix_size)],A_copy))
    b = [random.randint(0, 100) for i in range(len(P))]
    b_copy = b[:]
    for i in range(len(P)):
        b[i] = b_copy[P[i]]
    A_inv = inverse_matrix(A, P)
    x = solving_equation(A, b)
    #print(numpy.matrix(b_copy))
    print("\n Определитель матрицы A:  \n ", determinant(A))
    print(" \n Матрица b: ")
    for i in range(len(b)):
        print("b_{} = ".format(i + 1), b[i])
    print("\n Решение матричного уравнения:  \n ")
    for i in range(len(x)):
        print("x_{} = ".format(i + 1), x[i])
    print(check_solution(A_copy, x, b_copy))
    print("\n Обратная матрица:  \n ", numpy.matrix(A_inv))
    #check_inversematrix(A_copy, A_inv)
    print("\n Число обусловленности: ", conditional_number(A, A_inv))
    print("\n Матрица перестановок: \n", P)


main()
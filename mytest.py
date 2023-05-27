import SLEsolver
import numpy as np

def testJ(A, b, x0, n, m, epsilon, max_iter):
    return SLEsolver.Jacobi_method(A, b, x0, n, m, epsilon, max_iter)


def testS(A, b, x0, n, m, epsilon, max_iter):
    return SLEsolver.Seidel_method(A, b, x0, n, m, epsilon, max_iter)

A1 = [10, 8, -1, -2.8, 4, 1, -1, -0.6, 2]
b1 = [10, 60, 20]
x01 = [0, 0, 0]
n1 = 3
m1 = 3
epsilon1 = 0.01
max_iter1 = 100
print(testJ(A1, b1, x01, n1, m1, epsilon1, max_iter1))
print()
print(testS(A1, b1, x01, n1, m1, epsilon1, max_iter1))
print()


# Module
module with the jacobi and seidel method

Модуль для Python, написанный на языке С (стандарт С90). Предназначен для решения СЛАУ с помощью итерационных методов.
В модуле реализованы два метода (Якоби и Зейделя):
void jacobi_method(double *a, double *b, double *x0, int n, double epsilon, int max_iter)
void seidel(double *a, double *b, double *x0, int n, double epsilon, int max_iter)

На вход каждой функции передается матрица а, записанная как одномерный массив, правый столбец b, начальное приближение, размерность, желаемая точность, максимальное количество итераций.

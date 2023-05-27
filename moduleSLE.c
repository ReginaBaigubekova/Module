#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <python.h>
 
void jacobi_method(double *a, double *b, double *x0, int n, double epsilon, int max_iter) { 
    int i, j, k, iter; 
    double **A, *x, *x_new, diff, sum; 
 
    A = (double**) malloc(n * sizeof(double*)); 
    for (i = 0; i < n; i++) { 
        A[i] = (double*) malloc(n * sizeof(double)); 
        for (j = 0; j < n; j++) { 
            A[i][j] = a[i * n + j]; 
        } 
    } 
 
    x = (double*) malloc(n * sizeof(double)); 
    x_new = (double*) malloc(n * sizeof(double)); 
 
    for (i = 0; i < n; i++) { 
        x[i] = x0[i]; 
        x_new[i] = x0[i]; 
    } 
 
    iter = 0; 
    do { 
        iter++; 
        for (i = 0; i < n; i++) { 
            sum = 0; 
            for (j = 0; j < n; j++) { 
                if (i != j) { 
                    sum += A[i][j] * x_new[j]; 
                } 
            } 
            x_new[i] = (b[i] - sum) / A[i][i]; 
        } 
 
        diff = 0; 
        for (k = 0; k < n; k++) { 
            diff += pow(x_new[k] - x[k], 2); 
        } 
        diff = sqrt(diff); 
 
        for (k = 0; k < n; k++) { 
            x[k] = x_new[k]; 
        } 
    } while (diff > epsilon && iter < max_iter); 

    for (i = 0; i < n; i++) { 
        x0[i] = x[i]; 
    } 
 
    for (i = 0; i < n; i++) { 
        free(A[i]); 
    } 
    free(A); 
    free(x); 
    free(x_new); 
} 
 
static PyObject* Jacobi_method(PyObject* self, PyObject* args){

    PyObject *a;
    PyObject *b;
    PyObject *x0;
    int n;
    int m;
    int max_iter;
    double epsilon;

    if( !PyArg_ParseTuple( args, "OOOdddd", &a,&b,&x0,&n,&m,&epsilon,&max_iter) ) {
        PyErr_SetString(PyExc_TypeError, "The arguments you typed in are wrong !");
        return NULL;}
    int i,j;
    /*Преобразуем входные данные - начальное приближение к x0, в обычный одномерный массив Си*/
    double *x0_matrix = (double**)malloc(n * sizeof(double*));
    for (i=0; i<n; i++){
        PyObject* item = PyList_GetItem(x0, i); 
	    double elem = PyFloat_AsDouble(item); 
	    x0_matrix[i] = elem;
    }
    /*Преобразуем входные данные - матрицу b, в обычный одномерный массив Си*/
    double *b_matrix = (double**)malloc(n * sizeof(double*));
    for (i=0; i<n; i++){
        PyObject* item = PyList_GetItem(b, i); 
	    double elem = PyFloat_AsDouble(item); 
	    b_matrix[i] = elem;
    }
    /*Преобразуем входные данные - матрицу A, в обычный многомерный массив Си*/
    double* A_matrix;
    A_matrix = (double**) malloc(n*n*sizeof(double*));
    for (i = 0; i < n*n; i++){
        A_matrix[i] = PyFloat_AsDouble(PyList_GetItem(a,i));
    }
    /*Запускаем основную функцию simple_iterations, которая преобразует начальное приближение x0 в конечное приближение к решению*/
    jacobi_method(A_matrix,b_matrix,x0_matrix,n,epsilon,max_iter);
    /*Выводим на экран конечно приближение к решению задачи*/
    return Py_BuildValue("ddd",x0_matrix[0],x0_matrix[1],x0_matrix[2]);
}

void seidel_method(double *a, double *b, double* x0, int n, double epsilon) {
    int i, j, k, iter; 
    double **A, *x, *x_new, diff, sum; 
 
    A = (double**) malloc(n * sizeof(double*)); 
    for (i = 0; i < n; i++) { 
        A[i] = (double*) malloc(n * sizeof(double)); 
        for (j = 0; j < n; j++) { 
            A[i][j] = a[i * n + j]; 
        } 
    } 
    x = (double*) malloc(n * sizeof(double));
    x_new = (double*) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        x[i] = 0;
    }
    int converge = 0;
    while (!converge) {
        for (i = 0; i < n; i++) {
            int s1 = 0;
            for (j = 0; j < i; j++) {
                s1 += A[i][j] * x_new[j];
            }
            int s2 = 0;
            for (j = i + 1; j < n; j++) {
                s2 += A[i][j] * x[j];
            }
            x_new[i] = (b[i] - s1 - s2) / A[i][i];
        }
        converge = 1;
        for (i = 0; i < n; i++) {
            if (fabs(x_new[i] - x[i]) > epsilon) {
                converge = 0;
                break;
            }
        }
        for (i = 0; i < n; i++) {
            x[i] = x_new[i];
        }
        k++;
    }
    for (i = 0; i < n; i++) {
            x0[i] = x_new[i];
        }
    for (i = 0; i < n; i++) { 
        free(A[i]); 
    }
    free(x);
    free(x_new);
    free(A);
}

static PyObject* Seidel_method(PyObject* self, PyObject* args){

    PyObject *a;
    PyObject *b;
    PyObject *x0;
    int n;
    int m;
    int max_iter;
    double epsilon;

    if( !PyArg_ParseTuple( args, "OOOdddd", &a,&b,&x0,&n,&m,&epsilon,&max_iter) ) {
        PyErr_SetString(PyExc_TypeError, "The arguments you typed in are wrong !");
        return NULL;}
    int i,j;
    /*Преобразуем входные данные - начальное приближение к x0, в обычный одномерный массив Си*/
    double *x0_matrix = (double**)malloc(n * sizeof(double*));
    for (i=0; i<n; i++){
        PyObject* item = PyList_GetItem(x0, i); 
	    double elem = PyFloat_AsDouble(item); 
	    x0_matrix[i] = elem;
    }
    /*Преобразуем входные данные - матрицу b, в обычный одномерный массив Си*/
    double *b_matrix = (double**)malloc(n * sizeof(double*));
    for (i=0; i<n; i++){
        PyObject* item = PyList_GetItem(b, i); 
	    double elem = PyFloat_AsDouble(item); 
	    b_matrix[i] = elem;
    }
    /*Преобразуем входные данные - матрицу A, в обычный многомерный массив Си*/
    double* A_matrix;
    A_matrix = (double**) malloc(n*n*sizeof(double*));
    for (i = 0; i < n*n; i++){
        A_matrix[i] = PyFloat_AsDouble(PyList_GetItem(a,i));
    }
    seidel_method(A_matrix,b_matrix,x0_matrix,n,epsilon);
    /*Запускаем основную функцию simple_iterations, которая преобразует начальное приближение x0 в конечное приближение к решению*/
    /*Выводим на экран конечно приближение к решению задачи*/
    return Py_BuildValue("ddd",x0_matrix[0],x0_matrix[1],x0_matrix[2]);
}

/* =========================================================== */


static char Jacobi_method_docs[] =
    "Jacobi_method(A, b, x0, n, epsilon, max_iter) returns vector x \
    that solves a system of linear equations A*x = b by means of Jacobi's method, where x0 is the \
    initial value for x, max_iter is the maximum iterations number, \
    epsilon is maximum difference in value between exact solution and \
    approximated solution. \n";


static char Seidel_method_docs[] =
    "Seidel_method(A, b, x0, n, epsilon, max_iter) returns vector x \
    that solves a system of linear equations A*x = b by means of Seidel's method, where x0 is the \
    initial value for x, max_iter is the maximum iterations number, \
    epsilon is maximum difference in value between exact solution and \
    approximated solution. \n";


static PyMethodDef module_methods[] = 
{
    {"Jacobi_method", (PyCFunction)Jacobi_method, METH_VARARGS, Jacobi_method_docs},
    {"Seidel_method", (PyCFunction)Seidel_method, METH_VARARGS, Seidel_method_docs},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC PyInit_SLEsolver(void) {
    PyObject *module;

    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "SLEsolver",
        "Solving systems of linear equations\n",
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
    };

    module = PyModule_Create(&moduledef);
    if (!module) return NULL;

    return module;
}
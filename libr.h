#include <iostream>
#include <cmath>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>
#include <clocale>
#include <random>
using namespace std;

#define n 10

// Решение СЛАУ
double* gauss(double** matrix, double* b);

// Вывод в файл
void to_file(double *vector, int size, string filename);
// Вывод матрицы 
void print_matrix(double **matrix, int rows, int cols);
// Вывод  вектора
void print_vector(double *vector, int size);
// Локальная матрица жесткости
double **K_loc(double E, double J, double l);
//Зануление столбцов и строк
void make_zeros(int m, double ** K);
// Глобальная матрица жесткости
double **K_glob(double E, double J, double *lens, double k_coef);   
// Вектор F
double *F_glob(double P, double Mom, double q12, double q21, double* lens);

// Базисные функции
double N1(double x, double l); 
double N2(double x, double l);
double L1(double x, double l);
double L2(double x, double l);

// Первые производные б. ф-й
double d_N1(double x, double l);
double d_N2(double x, double l);
double d_L1(double x, double l);
double d_L2(double x, double l);

// Вторые производные б. ф-й
double d2_N1(double x, double l);
double d2_N2(double x, double l);
double d2_L1(double x, double l);
double d2_L2(double x, double l);

// Третьи производные б. ф-й
double d3_N1(double x, double l);
double d3_N2(double x, double l);
double d3_L1(double x, double l);
double d3_L2(double x, double l);

// Функции для распределенной нагрузки 
double q_b(double q12, double q21);
double q_a(double q12, double q21);

double fw_0(double l, double q1, double q0);
double ft_0(double l, double q1, double q0);
double fw_1(double l, double q1, double q0);
double ft_1(double l, double q1, double q0);

// Прогибы
double * w_x(double step, double* solution, double* lens, double* points);
// Углы
double * theta_x(double step, double* solution, double* lens, double* points);
// Моменты
double * M_x(double step, double* solution, double* lens, double* points, double E, double J);
// Силы
double * Q_x(double step, double* solution, double* lens, double* points, double E, double J);

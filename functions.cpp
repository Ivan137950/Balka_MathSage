#include "libr.h"

// Вывод в файл
void to_file(double *vector, int size, string filename) {
    ofstream file(filename);
    for (int i = 0; i < size; i++) {
        file << vector[i] << " ";
    }
    file << endl;
}

// Вывод матрицы 
void print_matrix(double **matrix, int rows, int cols) {
    const int column_width = 8;
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << setw(column_width) << matrix[i][j];
        }
        cout << endl;
    }
}

// Вывод  вектора
void print_vector(double *vector, int size) {
    const int column_width = 10;

    for (int i = 0; i < size; i++) {
        cout << setw(column_width) << vector[i];
    }
    cout << endl;
}


// Локальная матрица жесткости
double **K_loc(double E, double J, double l)   //(EI / L³) * | 12 6L -12 6L | | 6L 4L² -6L 2L² | | -12 -6L 12 -6L | | 6L 2L² -6L 4L² |
{
    double **K = new double*[4];
    for (int i = 0; i < 4; i++)
        K[i] = new double[4];
    
    K[0][0] = (12 * E * J) / (l * l * l);
    K[0][1] = (6 * E * J) / (l * l);
    K[0][2] = -(12 * E * J) / (l * l * l);
    K[0][3] = (6 * E * J) / (l * l);
    K[1][0] = (6 * E * J) / (l * l);
    K[1][1] = (4 * E * J) / (l);
    K[1][2] = -(6 * E * J) / (l * l);
    K[1][3] = (2 * E * J) / (l);
    K[2][0] = -(12 * E * J) / (l * l * l);
    K[2][1] = -(6 * E * J) / (l * l);
    K[2][2] = (12 * E * J) / (l * l * l);
    K[2][3] = -(6 * E * J) / (l * l);
    K[3][0] = (6 * E * J) / (l * l);
    K[3][1] = (2 * E * J) / (l);
    K[3][2] = -(6 * E * J) / (l * l);
    K[3][3] = (4 * E * J) / (l);
    return K;
}


//Зануление столбцов и строк
void make_zeros(int m, double ** K){
    for (int i = max(m - 3, 0); i <= m + 3; i++)
        if (i != m){
            K[i][m] = 0;
            K[m][i] = 0;
        }
    K[m][m] = 1;
}


// Глобальная матрица жесткости
double **K_glob(double E, double J, double *lens, double k_coef)   
{
    int N = n * 2;
    double **K = new double*[N];
    for (int i = 0; i < N; i++)
        K[i] = new double[N];
    for (int i = 0; i < n - 1; i++)
    {
        double l = lens[i];
        double **K_local = K_loc(E, J, l);
        int start = i * 2, end = start + 4;
        for (int j = start; j < end; j++){
            for (int k = start; k < end; k++){
                K[k][j] += K_local[j - start][k - start];
            }
        }
    }  
    double zeros[4]  = {2, 6, 8, 12};
    for (int i = 0; i < 4; i++)
        make_zeros(zeros[i], K);
    K[16][16] -= (k_coef);
    return K;
}


// Вектор F
double *F_glob(double P, double Mom, double q12, double q21, double* lens)
{
    double *F = new double[n * 2];
    for (int i = 0; i < n * 2; i++)
        F[i] = 0;
    F[0] = - P;
    F[5] = Mom;
    
    double a = q_a(q12, q21);
    double b = q_b(q12, q21);
    double q15 = a * 15 + b;

    // Добавление распределенной нагрузки
    for (int i = 0; i < 2; i++){
        double q0, q1;
        q0 = (1 - i) * q12 + i* q15;
        q1 = (1 - i) * q15 + i* q21;

        F[10 + 2 * i] -= fw_0(lens[5 + i], q1, q0);
        F[11 + 2 * i] -= ft_0(lens[5 + i], q1, q0);
        F[12 + 2 * i] -= fw_1(lens[5 + i], q1, q0);
        F[13 + 2 * i] -= ft_1(lens[5 + i], q1, q0);
    }
    F[12] = 0;
    return F;
}

// 10 11 12 13 14 15

// Базисные функции
double N1(double x, double l) { return   1 - 3 * (x / l) * (x / l)  + 2 * (x / l) * (x / l) * (x / l);}
double N2(double x, double l) { return   3 * (x / l) * (x / l)  - 2 * (x / l) * (x / l) * (x / l);}
double L1(double x, double l) { return   x - 2 * x * (x / l) + x * (x / l)  * (x / l);}
double L2(double x, double l) { return - x * (x / l) + x * (x / l)  * (x / l);}

// Первые производные б. ф-й
double d_N1(double x, double l) { return - 6 * (x / l) / l  + 6 * (x / l) * (x / l) / l;}
double d_N2(double x, double l) { return   6 * (x / l) / l  - 6 * (x / l) * (x / l) / l;}
double d_L1(double x, double l) { return   1 - 4 * (x / l) + 3 * (x / l)  * (x / l);}
double d_L2(double x, double l) { return - 2 * (x / l) + 3 * (x / l)  * (x / l);}

// Вторые производные б. ф-й
double d2_N1(double x, double l) { return - 6 / l / l  + 12 * (x / l) / l / l;}
double d2_N2(double x, double l) { return   6 / l / l  - 12 * (x / l) / l / l;}
double d2_L1(double x, double l) { return - 4 / l + 6 * (x / l) / l;}
double d2_L2(double x, double l) { return - 2 / l + 6 * (x / l) / l;}

// Третьи производные б. ф-й
double d3_N1(double x, double l) { return   12 / l / l / l;}
double d3_N2(double x, double l) { return - 12 / l / l / l;}
double d3_L1(double x, double l) { return   6 / l / l;}
double d3_L2(double x, double l) { return   6 / l / l;}

// Функции для распределенной нагрузки 
double q_b(double q12, double q21) { return q12 - 4 * (q21 - q12) / 3;}
double q_a(double q12, double q21) { return (q21 - q12) / 9;}

double fw_0(double l, double q1, double q0){ return  l / 20 * (7 * q0 + 3 * q1);}
double ft_0(double l, double q1, double q0){ return  l * l /60 * (3 * q0 + 2 * q1);}
double fw_1(double l, double q1, double q0){ return  l / 20 * (3 * q0 + 7 * q1);}
double ft_1(double l, double q1, double q0){ return -l * l / 60 * (2 * q0 + 3 * q1);}

// Прогибы
double * w_x(double step, double* solution, double* lens, double* points) {
    double x = 0;
    int ind = 0;
    double *w = new double[int(25 / step)];
    for (int i = 0; i < n - 1; i++) {
        while (x < points[i] + 0.00000001) {        
            double l = lens[i];
            double N1_val = N1(x - points[i] + l, l);
            double L1_val = L1(x - points[i] + l, l);
            double N2_val = N2(x - points[i] + l, l);
            double L2_val = L2(x - points[i] + l, l);

            w[ind] = N1_val * solution[i * 2] + L1_val * solution[i * 2 + 1] + N2_val * solution[i * 2 + 2] + L2_val * solution[i * 2 + 3];
            x += step;
            ind ++;
        }
        
    }
    return w;
}


// Углы
double * theta_x(double step, double* solution, double* lens, double* points) {
    double x = 0;
    int ind = 0;
    double *theta = new double[int(25 / step)];
    for (int i = 0; i < n - 1; i++) {
        while (x < points[i] + 0.00000001) {        
            double l = lens[i];
            double N1_val = d_N1(x - points[i] + l, l);
            double L1_val = d_L1(x - points[i] + l, l);
            double N2_val = d_N2(x - points[i] + l, l);
            double L2_val = d_L2(x - points[i] + l, l);

            theta[ind] = N1_val * solution[i * 2] + L1_val * solution[i * 2 + 1] + N2_val * solution[i * 2 + 2] + L2_val * solution[i * 2 + 3];
            x += step;
            ind ++;
        }
        
    }
    return theta;
}


// Моменты
double * M_x(double step, double* solution, double* lens, double* points, double E, double J) {
    double x = 0;
    int ind = 0;
    double *m = new double[int(25 / step)];
    for (int i = 0; i < n - 1; i++) {
        while (x < points[i] + 0.00000001) {        
            double l = lens[i];
            double N1_val = d2_N1(x - points[i] + l, l);
            double L1_val = d2_L1(x - points[i] + l, l);
            double N2_val = d2_N2(x - points[i] + l, l);
            double L2_val = d2_L2(x - points[i] + l, l);

            m[ind] = E * J * ( N1_val * solution[i * 2] + L1_val * solution[i * 2 + 1] + N2_val * solution[i * 2 + 2] + L2_val * solution[i * 2 + 3]);
            x += step;
            ind ++;
        }
        
    }
    return m;
}

// Силы
double * Q_x(double step, double* solution, double* lens, double* points, double E, double J) {
    double x = 0;
    int ind = 0;
    double *Q = new double[int(25 / step)];
    for (int i = 0; i < n - 1; i++) {
        while (x < points[i] + 0.00000001) {        
            double l = lens[i];
            double N1_val = d3_N1(x - points[i] + l, l);
            double L1_val = d3_L1(x - points[i] + l, l);
            double N2_val = d3_N2(x - points[i] + l, l);
            double L2_val = d3_L2(x - points[i] + l, l);

            Q[ind] = E * J * (N1_val * solution[i * 2] + L1_val * solution[i * 2 + 1] + N2_val * solution[i * 2 + 2] + L2_val * solution[i * 2 + 3]);
            x += step;
            ind ++;
        }
        
    }
    return Q;
}


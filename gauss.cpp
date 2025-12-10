#include "libr.h"

double* gauss(double** matrix, double* b) {
    int M = 2 * n;
    double** temp_matrix = new double*[M];
    for (int i = 0; i < M; i++) {
        temp_matrix[i] = new double[M + 1];
        for (int j = 0; j < M; j++) {
            temp_matrix[i][j] = matrix[i][j];
        }
        temp_matrix[i][M] = b[i];
    }

    // Прямой ход
    for (int i = 0; i < M; i++) {
        // Поиск строки с максимальным элементом в i-м столбце
        int max_row = i;
        for (int k = i + 1; k < M; k++) {
            if (fabs(temp_matrix[k][i]) > fabs(temp_matrix[max_row][i])) {
                max_row = k;
            }
        }

        // Обмен строк
        if (max_row != i) {
            for (int j = 0; j <= M; j++) {
                swap(temp_matrix[i][j], temp_matrix[max_row][j]);
            }
        }

        // Нормализация i-й строки
        double tmii = temp_matrix[i][i];
        for (int j = 0; j <= M; j++) {
            temp_matrix[i][j] /= tmii;
        }

        // Обнуление 
        for (int k = i + 1; k < M; k++) {
            double koef = temp_matrix[k][i];
            for (int j = 0; j < M+1; j++) {
                temp_matrix[k][j] -= koef * temp_matrix[i][j];
            }
        }
    }

    // Обратный ход
    double* x = new double[M];
    for (int i = M - 1; i >= 0; i--) {
        x[i] = temp_matrix[i][M];
        for (int j = i + 1; j < M; j++) {
            x[i] -= temp_matrix[i][j] * x[j];
        }
    }

    // for (int i = 0; i < M; i++) {
    //     cout << x[i] << " ";
    // }
    // cout << endl;

    // Освобождение памяти
    for (int i = 0; i < M; i++) {
        delete[] temp_matrix[i];
    }
    delete[] temp_matrix;

    return x;
}
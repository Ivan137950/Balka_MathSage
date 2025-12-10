#include "libr.h"

// MAIN FUNCTION
int main(int argc, char* argv[])
{   
    
    double P = std::atof(argv[1]), Mom = std::atof(argv[2]);
    double k_coef =  std::atof(argv[3]);
    double E =  std::atof(argv[4]), J =  std::atof(argv[5]);
    double q12 =  std::atof(argv[6]), q21 =  std::atof(argv[7]);


    double points[n-1] = {2, 5, 7, 10, 12, 15, 21, 23, 25};
    double lens[n-1] = {2, 3, 2, 3, 2, 3, 6, 2, 2};
    int N = n * 2;
    
    double **K = K_glob(E, J, lens, k_coef);
    double *F = F_glob(P, Mom, q12, q21, lens);
    double *U = gauss(K, F);

    print_matrix(K, N, N);
    cout << "F = " ;
    print_vector(F, N);
        cout << "U = " ;
    print_vector(U, N);

    double step = 0.001;
    double *w     = w_x(step, U, lens, points);
    double *theta = theta_x(step, U, lens, points);
    double *M     = M_x(step, U, lens, points, E, J);
    double *Q     = Q_x(step, U, lens, points, E, J);

    to_file(w, int(25 / step), "W.txt");
    to_file(theta, int(25 / step), "Theta.txt");
    to_file(M, int(25 / step), "M.txt");
    to_file(Q, int(25 / step), "Q.txt");


// ------------------------------
    for (int i = 0; i < N; i++)
        delete[] K[i];
    delete[] K;
    delete[] F;
    delete[] U;
    delete[] w;
    return 0;
}

// 2, 5, 7, 10, 12, 15, 21, 23, 25
// 2, 3, 2, 3,  2,  3,  6,  2,  2
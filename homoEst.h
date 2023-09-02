#include "iterativesvd.h"
#include "normalsvd.h"

void itr_svd(vector<vector<double>> A, vector<vector<double>> &U,vector<vector<double>> &V, vector<double> &S); 
void svd(vector<vector<double>> A, vector<vector<double>>& U, vector<double>& S, vector<vector<double>>& V);
vector<vector<double>> matrix_multiply(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int n = A.size();
    int m = B[0].size();
    int p = B.size();
    std::vector<std::vector<double>> C(n, std::vector<double>(m, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < p; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// function to calculate the normalization matrix
vector<vector<double>> normalization_matrix(vector<vector<double>> points) {
    int n = points.size();
    double sum_x = 0, sum_y = 0;
    double sum_d_x = 0, sum_d_y = 0;
    double mean_d, scale_factor;
    for (int i = 0; i < n; i++) {
        sum_x += points[i][0];
        sum_y += points[i][1];
    }
    double mean_x = sum_x / n;
    double mean_y = sum_y / n;
    for (int i = 0; i < n; i++) {
        sum_d_x += abs(points[i][0] - mean_x);
        sum_d_y += abs(points[i][1] - mean_y);
    }
    mean_d = (sum_d_x + sum_d_y) / (2 * n);
    scale_factor = sqrt(2) / mean_d;
    vector<vector<double>> T = {{scale_factor, 0, -scale_factor * mean_x},
                                {0, scale_factor, -scale_factor * mean_y},
                                {0, 0, 1}};
    return T;
}

// function to normalize feature points
vector<vector<double>> normalize_points(vector<vector<double>> points, 
                                         vector<vector<double>> T) {
    int n = points.size();
    vector<vector<double>> points_normalized(n, vector<double>(3, 1.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 2; j++) {
            points_normalized[i][j] = T[j][0] * points[i][0] + 
                T[j][1] * points[i][1] + T[j][2];
        }
    }
    return points_normalized;
}

// function to find homography by solving linear system of equations
vector<vector<double>> find_homography(vector<vector<double>> points1, 
                                        vector<vector<double>> points2) {
    int n = points1.size();
    vector<vector<double>> A(2 * n, vector<double>(9, 0.0));
    for (int i = 0; i < n/2; i++) {
        A[2 * i][0] = -points1[i][0];
        A[2 * i][1] = -points1[i][1];
        A[2 * i][2] = -1;
        A[2 * i][6] = points1[i][0] * points2[i][1];
        A[2 * i][7] = points1[i][1] * points2[i][1];
        A[2 * i][8] = points2[i][1];
        A[2 * i + 1][3] = -points1[i][0];
        A[2 * i + 1][4] = -points1[i][1];
        A[2 * i + 1][5] = -1;
        A[2 * i + 1][6] = points1[i][0] * points2[i][0];
        A[2 * i + 1][7] = points1[i][1] * points2[i][0];
        A[2 * i + 1][8] = points2[i][0];
    }
    vector<vector<double>> H(3, vector<double>(3, 0.0));
    vector<vector<double>> U(9,vector<double>(9,0)), S(9,vector<double>(9,0)), Vt(9,vector<double>(9,0));
    // svd(A, U, S, Vt);
    H[0][0] = Vt[8][0];
    H[0][1] = Vt[8][1];
    H[0][2] = Vt[8][2];
    H[1][0] = Vt[8][3];
    H[1][1] = Vt[8][4];
    H[1][2] = Vt[8][5];
    H[2][0] = Vt[8][6];
    H[2][1] = Vt[8][7];
    H[2][2] = Vt[8][8];
    return H;
}

vector<vector<double>> transpose(const vector<vector<double>>& A) {
    int m = A.size();
    int n = A[0].size();
    vector<vector<double>> result(n, vector<double>(m, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[i][j] = A[j][i];
        }
    }
    return result;
}

vector<vector<double>> reverese_Normalization(vector<vector<double>>H, 
                                        vector<vector<double>> T,vector<vector<double>> T_prime) {

    // Create B and B_prime as the inverse of T and T_prime
    vector<vector<double>> B(T.size(), vector<double>(T[0].size(), 0.0));
    vector<vector<double>> B_prime(T_prime.size(), vector<double>(T[0].size(), 0.0));

    for (int i = 0; i < 2; i++) {
        B[i][i] = 1.0 / T[i][i];
        B_prime[i][i] = 1.0 / T_prime[i][i];
    }

    // Compute H = T * H' * B_prime * B
    vector<vector<double>> temp1(3, vector<double>(3, 0.0));
    vector<vector<double>> temp2(3, vector<double>(3, 0.0));
    vector<vector<double>> H_new(3, vector<double>(3, 0.0));

    // Compute T * H'
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double val = 0.0;
            for (int k = 0; k < 3; k++) {
                val += T[i][k] * H[k][j];
            }
            temp1[i][j] = val;
        }
    }

    // Compute B_prime * B
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double val = 0.0;
            for (int k = 0; k < 3; k++) {
                val += B_prime[i][k] * B[k][j];
            }
            temp2[i][j] = val;
        }
    }

    // Compute H = T * H' * B_prime * B
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double val = 0.0;
            for (int k = 0; k < 3; k++) {
                val += temp1[i][k] * temp2[k][j];
            }
            H_new[i][j] = val;
        }
    }

    // Update H with the normalized version
    H = H_new;
    return H;

}

float homoEst(vector<vector<double>> points){
    vector<vector<double>> normalize_matrix = normalization_matrix(points);
    vector<vector<double>> normalized_points = normalize_points(points,normalize_matrix);
    vector<vector<double>> H_prime = find_homography(points,normalized_points);
    vector<vector<double>> H =  reverese_Normalization(H_prime,normalize_matrix,transpose(normalize_matrix));
    vector<vector<double>> A= matrix_multiply(H_prime,H);
    vector<vector<double>> U(A.size(),vector<double>(A[0].size(),0));
    vector<vector<double>> V(A[0].size(),vector<double>(A.size(),0));
    vector<double> S(0,A[0].size());


    auto begin = std::chrono::high_resolution_clock::now();
    // svd(A,U,S,V);

    itr_svd(A,U,V,S);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    
    // cout<<("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
    return elapsed.count() * 1e-9;
}


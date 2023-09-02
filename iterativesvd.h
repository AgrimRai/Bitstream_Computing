#include <iostream>
#include <bits/stdc++.h>
using namespace std;

// function to compute the dot product of two vectors
double dot_product(vector<double> a, vector<double> b) {
    double result = 0;
    for (int i = 0; i <(int)a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

// function to compute the norm of a vector
double norm(vector<double> a) {
    double result = 0;
    for (int i = 0; i < (int)a.size(); i++) {
        result += a[i] * a[i];
    }
    return sqrt(result);
}

// function to normalize a vector
vector<double> normalize(vector<double> a) {
    vector<double> result;
    double n = norm(a);
    for (int i = 0; i < (int)a.size(); i++) {
        result.push_back(a[i] / n);
    }
    return result;
}

// function to perform one iteration of the SVD algorithm
void svd_iteration(vector<vector<double>> A, vector<double> v, 
                   vector<double> &u, vector<double> &z, double &sigma) {
    vector<double> w(A.size(), 0.0);
    for (int i = 0; i <(int) A.size(); i++) {
        w[i] = dot_product(A[i], v);
    }
    // double alpha = dot_product(w, w);
    u = normalize(w);
    z.resize(A[0].size(), 0.0);
    for (int j = 0; j < (int)A[0].size(); j++) {
        for (int i = 0; i < (int)A.size(); i++) {
            z[j] += A[i][j] * u[i];
        }
    }
    sigma = norm(z);
    z = normalize(z);
}

// function to compute the SVD of a matrix A
void itr_svd(vector<vector<double>> A, vector<vector<double>> &U, 
         vector<vector<double>> &V, vector<double> &S) {
    int m = A.size();
    int n = A[0].size();
    vector<double> v(n, 1.0);
    vector<double> u(m, 0.0);
    vector<double> z(n, 0.0);
    double sigma = 0.0;
    int max_iter = min(m, n);
    for (int k = 0; k < max_iter; k++) {
        svd_iteration(A, v, u, z, sigma);
        U.push_back(u);
        V.push_back(z);
        S.push_back(sigma);
        vector<vector<double>> u_matrix(m, vector<double>(1, 0.0));
        for (int i = 0; i < m; i++) {
            u_matrix[i][0] = u[i];
        }
        vector<vector<double>> v_matrix(1, z);
        vector<vector<double>> sigma_matrix(1, vector<double>(1, sigma));
        vector<vector<double>> A_prime(m, vector<double>(n, 0.0));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A_prime[i][j] = A[i][j] - sigma_matrix[0][0] * 
                    u_matrix[i][0] * v_matrix[0][j];
            }
        }
        A = A_prime;
        v = z;
    }
}


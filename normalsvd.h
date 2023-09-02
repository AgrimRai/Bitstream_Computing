#include <cmath>
#include <vector>

using namespace std;

void svd(vector<vector<double>> A, vector<vector<double>>& U, vector<double>& S, vector<vector<double>>& V) {
   
        int m = A.size();
        int n = A[0].size();
        int nu = min(m, n);
        U = A;
        S.resize(nu);
        V.resize(n, vector<double>(n, 0.0));
        vector<double> rv1(n, 0.0);
        double scale, h, f, g, x, y, z, norm;
        int i, j, k, l, iter, s;
        for (i = 0; i < nu; i++) {
            l = i + 1;
            rv1[i] = scale * g;
            g = s = scale = 0.0;
            if (i < m) {
                for (k = i; k < m; k++) scale += abs(U[k][i]);
                if (scale != 0.0) {
                    for (k = i; k < m; k++) {
                        U[k][i] /= scale;
                        s += U[k][i] * U[k][i];
                    }
                    f = U[i][i];
                    g = -copysign(sqrt(s), f);
                    h = f * g - s;
                    U[i][i] = f - g;
                    for (j = l; j < n; j++) {
                        for (s = 0.0, k = i; k < m; k++) s += U[k][i] * U[k][j];
                        f = s / h;
                        for (k = i; k < m; k++) U[k][j] += f * U[k][i];
                    }
                    for (k = i; k < m; k++) U[k][i] *= scale;
                }
            }
            S[i] = scale * g;
            g = s = scale = 0.0;
            if (i < m && i != n - 1) {
                for (k = l; k < n; k++) scale += abs(U[i][k]);
                if (scale != 0.0) {
                    for (k = l; k < n; k++) {
                        U[i][k] /= scale;
                        s += U[i][k] * U[i][k];
                    }
                    f = U[i][l];
                    g = -copysign(sqrt(s), f);
                    h = f * g - s;
                    U[i][l] = f - g;
                    for (k = l; k < n; k++) rv1[k] = U[i][k] / h;
                    for (j = l; j < m; j++) {
                        for (s = 0.0, k = l; k < n; k++) s += U[j][k] * U[i][k];
                        for (k = l; k < n; k++) U[j][k] += s * rv1[k];
                    }
                    for (k = l; k < n; k++) U[i][k] *= scale;
                }
            }
            norm = max(norm, abs(S[i]) + abs(rv1[i]));
        
        for (i = nu - 1; i >= 0; i--) {
            if (i < n - 1) {
                if (S[i] != 0.0) {
                    for (j = l; j < n; j++) {
                        double t = 0.0;
                        for (k = l; k < n; k++) t += V[k][j] * U[i][k];
                        t = -t / U[i][i];
                        for (k = l; k < n; k++) V[k][j] += t * U[i][k];
                    }
                }
                for (j = l; j < n; j++) V[i][j] = V[j][i] = 0.0;
            }
            V[i][i] = 1.0;
            l = i;
        }
        zxy--;
    }

}    

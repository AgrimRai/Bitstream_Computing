
#include<bits/stdc++.h>
#include"iterativesvd.h"
#include"normalsvd.h"
using namespace std;

double dot_product(vector<double> a, vector<double> b);
double norm(vector<double> a) ;
vector<double> normalize(vector<double> a);
void svd_iteration(vector<vector<double>> A, vector<double> v,vector<double> &u, vector<double> &z, double &sigma);
void svd(vector<vector<double>> A, vector<vector<double>> &U, vector<vector<double>> &V, vector<double> &S) ;

void homoEst(vector<vector<double>> points);
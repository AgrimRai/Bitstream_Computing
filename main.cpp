#include<bits/stdc++.h>
#include <chrono>
#include "homoEst.h"
using namespace std;

int main(){
    // auto begin = std::chrono::high_resolution_clock::now();
    vector<vector<double>> points(9,vector<double>(9,0));
    points[0][0]=1;
    points[8][8]=2;
    float tt =0;
    for(int i=0;i<1000;i++){
        float time = homoEst(points);
        tt+=time;
    }
    cout<<"Time measured: "<<tt/1000;

    
    
}

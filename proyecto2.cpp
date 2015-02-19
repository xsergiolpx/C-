

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    /*int dimensionx = 9;
    int dimensiony = 5;
    double h=1;
    int CC = h*dimensionx;
    int RR = h*dimensiony;
    
    vector<vector<int> > matrix(RR, vector<int>(CC));
    
    for (int i = 0; i < RR; i++) {
        for (int j = 0; j < CC; j++) {
            matrix[i][j] = 0;
        }
    }*/
    
    double **M;
    new double*


    //contorno
    for (int i=0;i<RR;i++){
        matrix[i][0]=20;
        matrix[i][CC-1]=20;
    }

    
    
        for (int i = 0; i < RR; i++) {
        for (int j = 0; j < CC; j++) {
            cout << matrix[i][j] << " ";

        }
        cout << endl;
    }
    return 0;
}


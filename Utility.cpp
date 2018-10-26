//
// Created by Katarzyna on 06.10.2018.
//

#include "Utility.h"
#include <iostream>

using namespace std;
void Utility::testGrid(GRID A) {

    cout<<"| ------------------------------------------------------------"<<endl;
    cout<<"| GRID TEST CORNERS"<<endl;
    cout<<"| ------------------------------------------------------------"<<endl;
    cout << "(" << A.node[0].x << ", " << A.node[0].y << ")" << endl;
    cout << "(" << A.node[18].x << ", " << A.node[18].y << ")" << endl;
    cout << "(" << A.node[23].x << ", " << A.node[23].y << ")" << endl;
    cout << "(" << A.node[5].x << ", " << A.node[5].y << ")" << endl;
}
void Utility::print(GRID &ptGrid, int nH, int nL) {

    printNODE(ptGrid, nH, nL);
    printELEMENT(ptGrid, nH, nL);

}

void Utility::printNODE(GRID &grid, int nH, int nL) {

    cout<<"| ------------------------------------------------------------"<<endl;
    cout<<"| PRINT NODE"<<endl;
    cout<<"| ------------------------------------------------------------"<<endl;
    for (int i = 0; i < nH * nL; i++) {
        cout << "| idx: " << i << "  x = " << grid.node[i].x << "   y = " << grid.node[i].y << endl;
    }

}

void Utility::printELEMENT(GRID &grid, int nH, int nL) {

    cout<<"| ------------------------------------------------------------"<<endl;
    cout<<"| PRINT ELEMENT"<<endl;
    cout<<"| ------------------------------------------------------------"<<endl;
    for (int i = 0; i < (nH - 1) * (nL - 1); i++) {

        cout << "| idx: " << i << "  id = " << grid.element[i].id[0] << " " << grid.element[i].id[1] << " "
             << grid.element[i].id[2] << " " <<grid.element[i].id[3] << endl;
    }

}
void Utility::readFile(int *H, int *L, int *nH, int *nL, int *K, int *t) {
    ifstream read;
    read.open("../data.txt");
    if (read.is_open()) {
        read >> *H;
        read >> *L;
        read >> *nH;
        read >> *nL;
        read >> *K;
        read >> *t;

    }
    read.close();

    Utility::printData(H, L, nH, nL, K, t);
}
void Utility::printData(int *H, int *L, int *nH, int *nL, int *K, int *t)
{
    cout<<"| ------------------------------------------------------------"<<endl;
    cout<<"| PRINT DATA"<<endl;
    cout<<"| ------------------------------------------------------------"<<endl;
    cout << "| H = " << *H << endl;
    cout << "| L = " << *L << endl;
    cout << "| nH = " << *nH << endl;
    cout << "| nL = " << *nL << endl;
    cout << "| K = " << *K << endl;
    cout << "| t = " << *t << endl;
}
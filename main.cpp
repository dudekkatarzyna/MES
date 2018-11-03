#include <iostream>
#include<fstream>
#include"GRID.h"
#include "Utility.h"
using namespace std;



int main() {

    int H, L, nH, nL, K, t;

    Utility::readFile(&H, &L, &nH, &nL, &K, &t);
    GRID A(H, L, nH, nL, K, t);
    Utility::testGrid(A);

    GRID::createUniversalElement();

    GRID::createJacobian(A,nL, nH);

    return 0;
}
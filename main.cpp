#include <iostream>
#include<fstream>
#include <cmath>
#include"Grid/GRID.h"
#include "Utility/Utility.h"
#include "UniversalElements/UniversalElementSides.h"
#include "Utility/Calculations.h"

using namespace std;

int main() {
    float H, L;
    int nH, nL, K, t, alfa, c, ro;

    Utility::readFile(&H, &L, &nH, &nL, &K, &t, &alfa, &c, &ro);
    GRID A(H, L, nH, nL, K, t);
    //Utility::testGrid(A);

    GRID::setBC(A);

    float **globalH, **globalC, **globalP;
    float **matrixH, **matrixC;
    float **bordCond;

    for (int el = 0; el < ((nL - 1) * (nH - 1)); el++) {

        cout << endl << endl << "ELEMENT " << el << endl;
        matrixH = Calculations::createH(A, el);

        cout << endl << "MATRIX H WITHOUT BOARD CONDITIONS" << endl;
        Utility::printH(matrixH);
        bordCond = Calculations::addBordCondition(A, el, alfa);
        Calculations::sumArrays(matrixH, bordCond);

        Utility::printFinalH(matrixH);

        matrixC = Calculations::matrixC(c, ro);

        Utility::printMatrixC(matrixC);

        //float **vectorP=Calculations::vectorP(A, el, alfa);
    }

    globalH = new float *[nL * nH];
    globalC = new float *[nL * nH];

    for (int k = 0; k < nL * nH; ++k) {
        globalH[k] = new float[nL * nH]();
        globalC[k] = new float[nL * nH]();
    }

    for (int el = 0; el < ((nL - 1) * (nH - 1)); el++) {

        for (int i = 0; i < 4; i++) {

            for (int j = 0; j < 4; j++) {
                globalH[A.element[el].id[i]][A.element[el].id[j]] += matrixH[i][j];
                globalC[A.element[el].id[i]][A.element[el].id[j]] += matrixC[i][j];
            }
        }
    }

    Utility::printGlobalH(globalH, nH, nL);
    Utility::printGlobalC(globalC, nH, nL);

}
//
// Created by Katarzyna on 06.10.2018.
//

#include "Utility.h"
#include <iostream>

using namespace std;

void Utility::testGrid(GRID A) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| GRID TEST CORNERS" << endl;
    cout << "| ------------------------------------------------------------" << endl;
    cout << "(" << A.node[0].x << ", " << A.node[0].y << ")" << endl;
    cout << "(" << A.node[18].x << ", " << A.node[18].y << ")" << endl;
    cout << "(" << A.node[23].x << ", " << A.node[23].y << ")" << endl;
    cout << "(" << A.node[5].x << ", " << A.node[5].y << ")" << endl;
}

void Utility::printGrid(GRID &ptGrid, int nH, int nL) {

    printNODE(ptGrid, nH, nL);
    printELEMENT(ptGrid, nH, nL);

}

void Utility::printNODE(GRID &grid, int nH, int nL) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| PRINT NODE" << endl;
    cout << "| ------------------------------------------------------------" << endl;
    for (int i = 0; i < nH * nL; i++) {
        cout << "| idx: " << i << "  x = " << grid.node[i].x << "   y = " << grid.node[i].y << endl;
    }

}

void Utility::printELEMENT(GRID &grid, int nH, int nL) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| PRINT ELEMENT" << endl;
    cout << "| ------------------------------------------------------------" << endl;
    for (int i = 0; i < (nH - 1) * (nL - 1); i++) {

        cout << "| idx: " << i << "  id = " << grid.element[i].id[0] << " " << grid.element[i].id[1] << " "
             << grid.element[i].id[2] << " " << grid.element[i].id[3] << endl;
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

void Utility::printData(int *H, int *L, int *nH, int *nL, int *K, int *t) {
    cout << "| ------------------------------------------------------------" << endl;
    cout << "| PRINT DATA" << endl;
    cout << "| ------------------------------------------------------------" << endl;
    cout << "| H = " << *H << endl;
    cout << "| L = " << *L << endl;
    cout << "| nH = " << *nH << endl;
    cout << "| nL = " << *nL << endl;
    cout << "| K = " << *K << endl;
    cout << "| t = " << *t << endl;
}

void Utility::printUniversalElement(UniversalElement ue) {


    cout << "| ------------------------------------------------------------" << endl;
    cout << "| CREATE UNIVERSAL ELEMENT" << endl;
    cout << "| ------------------------------------------------------------" << endl;


    for (int i = 0; i < 4; i++) {
        cout << "divKsi punkt nr " << i << endl;
        for (int j = 0; j < 4; j++) {
            cout << ue.divKsi[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl;
    for (int i = 0; i < 4; i++) {
        cout << "divEta punkt nr " << i << endl;
        for (int j = 0; j < 4; j++) {
            cout << ue.divEta[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

}

void Utility::printCreateJacobian(float (*Jacobian)[2][2]) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| CREATE JACOBIAN" << endl;
    cout << "| ------------------------------------------------------------" << endl;

    printJacobian(Jacobian);

}

void Utility::printJacobian(float (*Jacobian)[2][2]) {


    for (int nodeId = 0; nodeId < 4; nodeId++) {
        cout << "node: " << nodeId << ": ";
        cout << Jacobian[nodeId][0][0] << " ";
        cout << Jacobian[nodeId][0][1] << " ";
        cout << Jacobian[nodeId][1][0] << " ";
        cout << Jacobian[nodeId][1][1] << endl;
    }
}

void Utility::printRevertJacobian(float (*Jacobian)[2][2]) {
    cout << "| ------------------------------------------------------------" << endl;
    cout << "| REVERT JACOBIAN" << endl;
    cout << "| ------------------------------------------------------------" << endl;

    printJacobian(Jacobian);
}

void Utility::printDetJ(float detJ[4]) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| DET JACOBIAN" << endl;
    cout << "| ------------------------------------------------------------" << endl;

    for (int nodeId = 0; nodeId < 4; nodeId++) {
        cout << "detJ:" << detJ[nodeId] << endl;
    }

}

void Utility::printMultiplyDetJacobian(float (*Jacobian)[2][2]) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| MULTIPLY JACOBIAN AND 1/DET" << endl;
    cout << "| ------------------------------------------------------------" << endl;

    Utility::printJacobian(Jacobian);

}

void Utility::printdNdXY(float (*divNx)[4], float (*divNy)[4]) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| CALCULATE dN/dx AND dN/dy" << endl;
    cout << "| ------------------------------------------------------------" << endl;


    cout << "dN / dx" << endl;
    for (int nodeId = 0; nodeId < 4; nodeId++) {
        cout << divNx[nodeId][0] << " " << divNx[nodeId][1] << " " << divNx[nodeId][2] << " " << divNx[nodeId][3]
             << endl;

    }
    cout << "dN / dy" << endl;
    for (int nodeId = 0; nodeId < 4; nodeId++) {

        cout << divNy[nodeId][0] << " " << divNy[nodeId][1] << " " << divNy[nodeId][2] << " " << divNy[nodeId][3]
             << endl;
    }

}

void Utility::printdivNxySqr(float (*divNxSqr)[4][4], float (*divNySqr)[4][4]) {


    cout << endl << "{dN/dx}{dN/dx}T" << endl;
    for (int pc = 0; pc < 4; pc++) {
        cout << endl << "PC: " << pc << endl;
        cout << divNxSqr[pc][0][0] << " " << divNxSqr[pc][0][1] << " " << divNxSqr[pc][0][2] << " "
             << divNxSqr[pc][0][3]
             << endl;
        cout << divNxSqr[pc][1][0] << " " << divNxSqr[pc][1][1] << " " << divNxSqr[pc][1][2] << " "
             << divNxSqr[pc][1][3]
             << endl;
        cout << divNxSqr[pc][2][0] << " " << divNxSqr[pc][2][1] << " " << divNxSqr[pc][2][2] << " "
             << divNxSqr[pc][2][3]
             << endl;
        cout << divNxSqr[pc][3][0] << " " << divNxSqr[pc][3][1] << " " << divNxSqr[pc][3][2] << " "
             << divNxSqr[pc][3][3];
        cout << endl;
    }
    cout << endl << "{dN/dy}{dN/dy}T" << endl;
    for (int pc = 0; pc < 4; pc++) {
        cout << endl << "PC: " << pc << endl;
        cout << divNySqr[pc][0][0] << " " << divNySqr[pc][0][1] << " " << divNySqr[pc][0][2] << " "
             << divNySqr[pc][0][3]
             << endl;
        cout << divNySqr[pc][1][0] << " " << divNySqr[pc][1][1] << " " << divNySqr[pc][1][2] << " "
             << divNySqr[pc][1][3]
             << endl;
        cout << divNySqr[pc][2][0] << " " << divNySqr[pc][2][1] << " " << divNySqr[pc][2][2] << " "
             << divNySqr[pc][2][3]
             << endl;
        cout << divNySqr[pc][3][0] << " " << divNySqr[pc][3][1] << " " << divNySqr[pc][3][2] << " "
             << divNySqr[pc][3][3];
        cout << endl;
    }

}

void Utility::printMultiplyT(float (*divNxSqr)[4][4], float (*divNySqr)[4][4]) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| MULTIPLY TRANSPOSITION" << endl;
    cout << "| ------------------------------------------------------------" << endl;

    Utility::printdivNxySqr(divNxSqr, divNySqr);
}

void Utility::printRemoveIntegral(float (*divNxSqr)[4][4], float (*divNySqr)[4][4]) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| REMOVE INTEGRAL" << endl;
    cout << "| ------------------------------------------------------------" << endl;

    Utility::printdivNxySqr(divNxSqr, divNySqr);
}

void Utility::printMultiplyK(float (*arrK)[4][4]) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| MULTIPLY K" << endl;
    cout << "| ------------------------------------------------------------" << endl;

    cout << endl;
    for (int pkt = 0; pkt < 4; pkt++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                cout << arrK[pkt][i][j] << " ";
            }
            cout << endl;
        }
        cout << endl << endl;
    }


}

void Utility::printH(float (*H)[4]) {

    cout << "| ------------------------------------------------------------" << endl;
    cout << "| CALCULATE H" << endl;
    cout << "| ------------------------------------------------------------" << endl;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << H[i][j] << " ";
        }
        cout << endl;
    }

}

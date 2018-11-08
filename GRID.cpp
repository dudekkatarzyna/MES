//
// Created by Katarzyna on 06.10.2018.
//

#include "GRID.h"
#include "Utility.h"
#include "UniversalElement.h"
#include<iostream>
#include <cmath>

using namespace std;

GRID::GRID(int H, int L, int nH, int nL, int K, int t) {

    node = new NODE[nH * nL];
    element = new ELEMENT[(nH - 1) * (nL - 1)];

    int i = 0;
    int deltaX = L / (nL - 1);
    int deltaY = H / (nH - 1);
    for (int l = 0; l < nL; l++) {
        for (int h = 0; h < nH; h++) {
            this->node[i].x = l * deltaX;
            this->node[i].y = h * deltaY;
            this->node[i].t = t;
            i++;
        }
    }

    int j = 0;
    for (int i = 0; i < (nH - 1) * (nL - 1); i++) {
        if (i % (nH - 1) == 0 && j != 0) {
            j++;
        }
        this->element[i].id[0] = j;
        this->element[i].id[1] = j + nH;
        this->element[i].id[2] = j + nH + 1;
        this->element[i].id[3] = j + 1;

        this->element[i].k = K;
        j++;
    }

    Utility::printGrid(*this, nH, nL);
}

float Jacobian[4][2][2];
float detJ[4];
UniversalElement ue;

void GRID::createH(GRID A, int elId) {

    Utility::printUniversalElement(ue);
    GRID::createJacobian(A, elId);

}

void GRID::createJacobian(GRID A, int elId) {


    for (int nodeId = 0; nodeId < 4; nodeId++) {
        Jacobian[nodeId][0][0] =
                ue.divKsi[nodeId][0] * A.node[A.element[elId].id[0]].x +
                ue.divKsi[nodeId][1] * A.node[A.element[elId].id[1]].x +
                ue.divKsi[nodeId][2] * A.node[A.element[elId].id[2]].x +
                ue.divKsi[nodeId][3] * A.node[A.element[elId].id[3]].x;
        Jacobian[nodeId][0][1] =
                ue.divKsi[nodeId][0] * A.node[A.element[elId].id[0]].y +
                ue.divKsi[nodeId][1] * A.node[A.element[elId].id[1]].y +
                ue.divKsi[nodeId][2] * A.node[A.element[elId].id[2]].y +
                ue.divKsi[nodeId][3] * A.node[A.element[elId].id[3]].y;
        Jacobian[nodeId][1][0] =
                ue.divEta[nodeId][0] * A.node[A.element[elId].id[0]].x +
                ue.divEta[nodeId][1] * A.node[A.element[elId].id[1]].x +
                ue.divEta[nodeId][2] * A.node[A.element[elId].id[2]].x +
                ue.divEta[nodeId][3] * A.node[A.element[elId].id[3]].x;
        Jacobian[nodeId][1][1] =
                ue.divEta[nodeId][0] * A.node[A.element[elId].id[0]].y +
                ue.divEta[nodeId][1] * A.node[A.element[elId].id[1]].y +
                ue.divEta[nodeId][2] * A.node[A.element[elId].id[2]].y +
                ue.divEta[nodeId][3] * A.node[A.element[elId].id[3]].y;
    }

    Utility::printCreateJacobian(Jacobian);
    GRID::detJacobian(A);

}

void GRID::detJacobian(GRID A) {

    for (int nodeId = 0; nodeId < 4; nodeId++) {
        detJ[nodeId] =
                Jacobian[nodeId][0][0] * Jacobian[nodeId][1][1] - Jacobian[nodeId][0][1] * Jacobian[nodeId][1][0];
    }
    Utility::printDetJ(detJ);
    GRID::revertJavobian(A);

}

void GRID::revertJavobian(GRID A) {

    for (int nodeId = 0; nodeId < 4; nodeId++) {
        float tmp = Jacobian[nodeId][0][0];
        Jacobian[nodeId][0][0] = Jacobian[nodeId][1][1];
        Jacobian[nodeId][1][1] = tmp;
        Jacobian[nodeId][0][1] *= (-1);
        Jacobian[nodeId][1][0] *= (-1);
    }

    Utility::printRevertJacobian(Jacobian);
    GRID::multiplyDetJacobian(A);

}


void GRID::multiplyDetJacobian(GRID A) {

    for (int nodeId = 0; nodeId < 4; nodeId++) {
        Jacobian[nodeId][0][0] *= (1 / detJ[nodeId]);
        Jacobian[nodeId][0][1] *= (1 / detJ[nodeId]);
        Jacobian[nodeId][1][0] *= (1 / detJ[nodeId]);
        Jacobian[nodeId][1][1] *= (1 / detJ[nodeId]);
    }

    Utility::printMultiplyDetJacobian(Jacobian);
    GRID::dNdXY(A);
}

void GRID::dNdXY(GRID A) {

    float divNx[4][4];
    float divNy[4][4];

    for (int nodeId = 0; nodeId < 4; nodeId++) {

        divNx[nodeId][0] =
                Jacobian[nodeId][0][0] * ue.divKsi[nodeId][0] + Jacobian[nodeId][0][1] * ue.divEta[nodeId][0];
        divNx[nodeId][1] =
                Jacobian[nodeId][0][0] * ue.divKsi[nodeId][1] + Jacobian[nodeId][0][1] * ue.divEta[nodeId][1];
        divNx[nodeId][2] =
                Jacobian[nodeId][0][0] * ue.divKsi[nodeId][2] + Jacobian[nodeId][0][1] * ue.divEta[nodeId][2];
        divNx[nodeId][3] =
                Jacobian[nodeId][0][0] * ue.divKsi[nodeId][3] + Jacobian[nodeId][0][1] * ue.divEta[nodeId][3];

        divNy[nodeId][0] =
                Jacobian[nodeId][1][0] * ue.divKsi[nodeId][0] + Jacobian[nodeId][1][1] * ue.divEta[nodeId][0];
        divNy[nodeId][1] =
                Jacobian[nodeId][1][0] * ue.divKsi[nodeId][1] + Jacobian[nodeId][1][1] * ue.divEta[nodeId][1];
        divNy[nodeId][2] =
                Jacobian[nodeId][1][0] * ue.divKsi[nodeId][2] + Jacobian[nodeId][1][1] * ue.divEta[nodeId][2];
        divNy[nodeId][3] =
                Jacobian[nodeId][1][0] * ue.divKsi[nodeId][3] + Jacobian[nodeId][1][1] * ue.divEta[nodeId][3];

    }

    Utility::printdNdXY(divNx, divNy);

    multiplyT(A, divNx, divNy);
}

void GRID::multiplyT(GRID A, float divNx[4][4], float divNy[4][4]) {


    float divNxSqr[4][4][4];
    float divNySqr[4][4][4];

    for (int pc = 0; pc < 4; pc++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                divNxSqr[pc][i][j] = divNx[pc][i] * divNx[pc][j];

                divNySqr[pc][i][j] = divNy[pc][i] * divNy[pc][j];
            }
        }
    }

    Utility::printMultiplyT(divNxSqr, divNySqr);
    GRID::removeIntegral(A, divNxSqr, divNySqr);
}


void GRID::removeIntegral(GRID A, float divNxSqr[4][4][4], float divNySqr[4][4][4]) {

    for (int arr = 0; arr < 4; arr++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                divNxSqr[arr][i][j] *= detJ[i];
                divNySqr[arr][i][j] *= detJ[i];

            }
        }
    }

    Utility::printRemoveIntegral(divNxSqr, divNySqr);
    GRID::multiplyK(A, divNxSqr, divNySqr);
}

void GRID::multiplyK(GRID A, float divNxSqr[4][4][4], float divNySqr[4][4][4]) {

    float arrK[4][4][4];
    for (int pkt = 0; pkt < 4; pkt++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                arrK[pkt][i][j] = A.element[pkt].k * (divNxSqr[pkt][i][j] + divNySqr[pkt][i][j]);
            }
        }
    }

    Utility::printMultiplyK(arrK);
    GRID::calculateH(arrK);
}

void GRID::calculateH(float arrK[4][4][4]) {

    float H[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            H[i][j] = 0;
            for (int arr = 0; arr < 4; arr++) {
                H[i][j] += arrK[arr][i][j];
            }
        }
    }

    Utility::printH(H);
}
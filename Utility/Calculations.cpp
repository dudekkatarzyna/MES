//
// Created by Katarzyna on 12.11.2018.
//

#include <cmath>
#include <iostream>
#include "Calculations.h"
#include "../UniversalElements/UniversalElement.h"
#include "Utility.h"
#include "../UniversalElements/UniversalElementSides.h"
#include "../UniversalElements/KsiEta.h"

using namespace std;
float Jacobian[4][2][2];
float **H;
float detJ[4];
float divNx[4][4];
float divNy[4][4];
float divNxSqr[4][4][4];
float divNySqr[4][4][4];
float arrK[4][4][4];
float **PC1Sqr;
float **PC2Sqr;
float **sum;
float sideLenght = 0;
float PC1[4], PC2[4];
UniversalElement ue;

float **Calculations::createH(GRID A, int elId) {


    Calculations::createJacobian(A, elId);
    Calculations::detJacobian();
    Calculations::revertJavobian();
    Calculations::multiplyDetJacobian();
    Calculations::dNdXY();
    multiplyT(divNx, divNy);
    Calculations::removeIntegral(divNxSqr, divNySqr);
    Calculations::multiplyK(A, divNxSqr, divNySqr);
    Calculations::calculateH(arrK);

    return H;
}

void Calculations::createJacobian(GRID A, int elId) {


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

    //Utility::printCreateJacobian(Jacobian);
}

void Calculations::detJacobian() {

    for (int nodeId = 0; nodeId < 4; nodeId++) {
        detJ[nodeId] =
                Jacobian[nodeId][0][0] * Jacobian[nodeId][1][1] - Jacobian[nodeId][0][1] * Jacobian[nodeId][1][0];
    }
    //Utility::printDetJ(detJ);
}

void Calculations::revertJavobian() {

    for (auto &nodeId : Jacobian) {
        float tmp = nodeId[0][0];
        nodeId[0][0] = nodeId[1][1];
        nodeId[1][1] = tmp;
        nodeId[0][1] *= (-1);
        nodeId[1][0] *= (-1);
    }

    //Utility::printRevertJacobian(Jacobian);


}


void Calculations::multiplyDetJacobian() {

    for (int nodeId = 0; nodeId < 4; nodeId++) {
        Jacobian[nodeId][0][0] *= (1 / detJ[nodeId]);
        Jacobian[nodeId][0][1] *= (1 / detJ[nodeId]);
        Jacobian[nodeId][1][0] *= (1 / detJ[nodeId]);
        Jacobian[nodeId][1][1] *= (1 / detJ[nodeId]);
    }

    //Utility::printMultiplyDetJacobian(Jacobian);

}

void Calculations::dNdXY() {

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

    //Utility::printdNdXY(divNx, divNy);
}

void Calculations::multiplyT(float divNx[4][4], float divNy[4][4]) {


    for (int pc = 0; pc < 4; pc++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                divNxSqr[pc][i][j] = divNx[pc][i] * divNx[pc][j];
                divNySqr[pc][i][j] = divNy[pc][i] * divNy[pc][j];
            }
        }
    }

    //Utility::printMultiplyT(divNxSqr, divNySqr);
}


void Calculations::removeIntegral(float divNxSqr[4][4][4], float divNySqr[4][4][4]) {

    for (int arr = 0; arr < 4; arr++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                divNxSqr[arr][i][j] *= detJ[i];
                divNySqr[arr][i][j] *= detJ[i];
            }
        }
    }

    //Utility::printRemoveIntegral(divNxSqr, divNySqr);
}

void Calculations::multiplyK(GRID A, float divNxSqr[4][4][4], float divNySqr[4][4][4]) {


    for (int pkt = 0; pkt < 4; pkt++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                arrK[pkt][i][j] = A.element[pkt].k * (divNxSqr[pkt][i][j] + divNySqr[pkt][i][j]);
            }
        }
    }

    //Utility::printMultiplyK(arrK);
}


void Calculations::calculateH(float arrK[4][4][4]) {

    H = new float *[4];

    for (int i = 0; i < 4; i++) {
        H[i] = new float[4];
        for (int j = 0; j < 4; j++) {
            H[i][j] = 0;
            for (int arr = 0; arr < 4; arr++) {
                H[i][j] += arrK[arr][i][j];
            }
        }
    }

    //Utility::printH(H);
}

float **Calculations::addBordCondition(GRID A, int el, int alfa) {

    auto **sumArr = new float *[4];
    for (int i = 0; i < 4; ++i) {
        sumArr[i] = new float[4]{};
    }
    UniversalElementSides a;

    sideLenght = static_cast<float>(sqrt(
            pow((A.node[A.element[el].id[2]].x - A.node[A.element[el].id[1]].x), 2) +
            pow((A.node[A.element[el].id[2]].y - A.node[A.element[el].id[1]].y), 2)));

    if (A.element[el].Q[0]) {
        a.N(PC1, PC2, 0);
        makeArrayFromVector(PC1, PC2, alfa);
        sumArrays(sumArr, sumArraysDet(PC1Sqr, PC2Sqr));
    }
    if (A.element[el].Q[1]) {
        a.N(PC1, PC2, 1);
        makeArrayFromVector(PC1, PC2, alfa);
        sumArrays(sumArr, sumArraysDet(PC1Sqr, PC2Sqr));
    }
    if (A.element[el].Q[2]) {
        a.N(PC1, PC2, 2);
        makeArrayFromVector(PC1, PC2, alfa);
        sumArrays(sumArr, sumArraysDet(PC1Sqr, PC2Sqr));
    }
    if (A.element[el].Q[3]) {
        a.N(PC1, PC2, 3);
        makeArrayFromVector(PC1, PC2, alfa);
        sumArrays(sumArr, sumArraysDet(PC1Sqr, PC2Sqr));
    }


    // Utility::printSumBC(sumArr);
    return sumArr;
}

void Calculations::makeArrayFromVector(float PC1[4], float PC2[4], int alfa) {

    PC1Sqr = new float *[4];
    PC2Sqr = new float *[4];
    for (int i = 0; i < 4; i++) {
        PC1Sqr[i] = new float[4];
        PC2Sqr[i] = new float[4];
        for (int j = 0; j < 4; j++) {
            PC1Sqr[i][j] = PC1[i] * PC1[j] * alfa;
            PC2Sqr[i][j] = PC2[i] * PC2[j] * alfa;
        }
    }

    //  Utility::printNxN(PC1Sqr);
    //  Utility::printNxN(PC2Sqr);
}

float **Calculations::sumArraysDet(float **Arr1, float **Arr2) {

    sum = new float *[4];

    for (int i = 0; i < 4; i++) {
        sum[i] = new float[4];
        for (int j = 0; j < 4; j++) {
            sum[i][j] = (Arr1[i][j] + Arr2[i][j]) * (sideLenght / 2);
        }
    }

//    Utility::printSumPC(sum);
    return sum;
}

void Calculations::sumArrays(float **Arr1, float **Arr2) {

    sum = new float *[4];

    for (int i = 0; i < 4; i++) {
        sum[i] = new float[4];
        for (int j = 0; j < 4; j++) {
            Arr1[i][j] += Arr2[i][j];
        }
    }
}

float **MatrixCN;
float ***MatrixCNSqrt;
float **MatrixC;

float **Calculations::matrixC(int c, int ro) {

    Calculations::createNKsiEta();
    Calculations::createMatrixCNSqrt();

    MatrixC = new float *[4];
    for (int i = 0; i < 4; i++) {
        MatrixC[i] = new float[4];
        for (int j = 0; j < 4; ++j) {
            MatrixC[i][j] =
                    (MatrixCNSqrt[0][i][j] + MatrixCNSqrt[1][i][j] + MatrixCNSqrt[2][i][j] + MatrixCNSqrt[3][i][j]) *
                    c * ro * detJ[i];
        }
    }
    return MatrixC;
}


void Calculations::createNKsiEta() {


    KsiEta elArr[4];
    elArr[0].ksi = static_cast<float>(-1 / sqrt(3));
    elArr[0].eta = static_cast<float>(-1 / sqrt(3));
    elArr[1].ksi = static_cast<float>(1 / sqrt(3));
    elArr[1].eta = static_cast<float>(-1 / sqrt(3));
    elArr[2].ksi = static_cast<float>(1 / sqrt(3));
    elArr[2].eta = static_cast<float>(1 / sqrt(3));
    elArr[3].ksi = static_cast<float>(-1 / sqrt(3));
    elArr[3].eta = static_cast<float>(1 / sqrt(3));


    MatrixCN = new float *[4];
    for (int i = 0; i < 4; i++) {
        MatrixCN[i] = new float[4];
        MatrixCN[i][0] = static_cast<float>(0.25 * (1 - elArr[i].ksi) * (1 - elArr[i].eta));
        MatrixCN[i][1] = static_cast<float>(0.25 * (1 + elArr[i].ksi) * (1 - elArr[i].eta));
        MatrixCN[i][2] = static_cast<float>(0.25 * (1 + elArr[i].ksi) * (1 + elArr[i].eta));
        MatrixCN[i][3] = static_cast<float>(0.25 * (1 - elArr[i].ksi) * (1 + elArr[i].eta));
    }

    //Utility::printMatrix4x4(MatrixCN);
}


void Calculations::createMatrixCNSqrt() {

    MatrixCNSqrt = new float **[4];

    for (int pc = 0; pc < 4; pc++) {
        MatrixCNSqrt[pc] = new float *[4];
        for (int i = 0; i < 4; i++) {
            MatrixCNSqrt[pc][i] = new float[4];
            for (int j = 0; j < 4; j++) {
                MatrixCNSqrt[pc][i][j] = MatrixCN[pc][i] * MatrixCN[pc][j];
            }
        }
        //Utility::printNxNinC(MatrixCNSqrt[pc]);
    }

}

float **Calculations::vectorP(GRID A, int el, int alfa, int tA) {

    auto **vectorP = new float *[4];
    for (int i = 0; i < 4; ++i) {
        vectorP[i] = new float[4]{};
    }

    UniversalElementSides b;

    if (A.element[el].Q[0]) {
        b.N(PC1, PC2, 0);
        for (int i = 0; i < 4; i++) {
            vectorP[i][0] += (PC1[i] + PC2[i]);
        }
    }
    if (A.element[el].Q[1]) {
        b.N(PC1, PC2, 1);
        for (int i = 0; i < 4; i++) {
            vectorP[i][0] += (PC1[i] + PC2[i]);
        }
    }
    if (A.element[el].Q[2]) {
        b.N(PC1, PC2, 2);
        for (int i = 0; i < 4; i++) {
            vectorP[i][0] += (PC1[i] + PC2[i]);
        }
    }
    if (A.element[el].Q[3]) {
        b.N(PC1, PC2, 3);
        for (int i = 0; i < 4; i++) {
            vectorP[i][0] += (PC1[i] + PC2[i]);
        }
    }
    for (int i = 0; i < 4; i++) {
        vectorP[i][0] *= (alfa * tA * 0.5 * sideLenght);
    }

    return vectorP;
}

float **Calculations::solveEqationForT(float **H, float **P, int nH, int nL, float **t1) {

    const double eps = 1e-12;

    auto **HP = new float *[nL * nH];
    for (int k = 0; k < nL * nH; ++k) {
        HP[k] = new float[(nL * nH) + 1]();
    }

    for (int i = 0; i < nH * nL; ++i) {
        for (int j = 0; j < nH * nL; ++j) {
            HP[i][j] = H[i][j];
        }
        HP[i][nH * nL] = P[i][0];
    }

    float m, s;

    // eliminacja współczynników
    for (int i = 0; i < nL * nH; i++) {
        for (int j = i + 1; j < nL * nH; j++) {
            if (fabs(H[i][i]) < eps) break;
            m = -HP[j][i] / HP[i][i];
            for (int k = i + 1; k <= nL * nH; k++)
                HP[j][k] += m * HP[i][k];
        }
    }

    // wyliczanie niewiadomych
    for (int i = (nL * nH) - 1; i >= 0; i--) {
        s = HP[i][nL * nH];
        for (int j = nL * nH - 1; j >= i + 1; j--)
            s -= HP[i][j] * t1[j][0];
        if (fabs(HP[i][i]) < eps) break;
        t1[i][0] = s / HP[i][i];
    }

    return t1;
}

void Calculations::getMinMaxtemp(float **t1Vector, int nH, int nL) {

    float min = t1Vector[0][0], max = t1Vector[0][0];

    for (int i = 1; i < nH * nL; i++) {
        if (t1Vector[i][0] < min) min = t1Vector[i][0];
        else if (t1Vector[i][0] > max) max = t1Vector[i][0];
    }

    Utility::printMinMaxTemp(min, max);
}

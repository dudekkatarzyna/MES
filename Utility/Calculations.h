//
// Created by Katarzyna on 12.11.2018.
//

#ifndef MES_CALCULATIONS_H
#define MES_CALCULATIONS_H


#include "../Grid/GRID.h"

class Calculations {

public:
    static float** createH(GRID A, int elId);
    static void createJacobian(GRID A, int elId);
    static void detJacobian(GRID A);
    static void multiplyDetJacobian(GRID A);
    static void revertJavobian(GRID A);
    static void dNdXY(GRID A);
    static void multiplyT(GRID A, float divNx[4][4], float divNy[4][4]);
    static void multiplyK(GRID A, float divNxSqr[4][4][4], float divNySqr[4][4][4]);
    static void removeIntegral(GRID A, float divNxSqr[4][4][4], float divNySqr[4][4][4]);
    static void calculateH(float arrK[4][4][4]);
    static float **addBordCondition(GRID A, int el, int alfa);
    static void makeArrayFromVector(float PC1[4], float PC2[4], int alfa, float sideLenght);
    static float** sumArraysDet(float **Arr1, float **Arr2);
    static void sumArrays(float **Arr1, float **Arr2);

    static float** vectorP(GRID A, int el, int alfa, int t0);

    static float** matrixC(int c, int ro);
    static void createNKsiEta();
    static void createMatrixCNSqrt();

    static float** solveEqationForT(float **H, float **P, int nH, int nL, float **t1);

    static void getMinMaxtemp(float **t1Vector, int nH, int nL);
};


#endif //MES_CALCULATIONS_H

//
// Created by Katarzyna on 12.11.2018.
//

#ifndef MES_CALCULATIONS_H
#define MES_CALCULATIONS_H


#include "../Grid/GRID.h"

class Calculations {

public:
    static double** createH(GRID A, int elId);
    static void createJacobian(GRID A, int elId);
    static void detJacobian();
    static void multiplyDetJacobian();
    static void revertJavobian();
    static void dNdXY();
    static void multiplyT(double divNx[4][4], double divNy[4][4]);
    static void multiplyK(GRID A, double divNxSqr[4][4][4], double divNySqr[4][4][4], int elId);
    static void removeIntegral(double divNxSqr[4][4][4], double divNySqr[4][4][4]);
    static void calculateH(double arrK[4][4][4]);
    static double **addBordCondition(GRID A, int el, double alfa);
    static void makeArrayFromVector(double PC1[4], double PC2[4], double alfa);
    static double** sumArraysDet(double **Arr1, double **Arr2);
    static void sumArrays(double **Arr1, double **Arr2);
    static double** vectorP(GRID A, int el, double alfa, double tA);
    static double** matrixC(double c, double ro);
    static void createNKsiEta();
    static void createMatrixCNSqrt();
    static double** solveEqationForT(double **H, double **P, int nH, int nL, double **t1);
    static void getMinMaxtemp(double **t1Vector, int nH, int nL);
};


#endif //MES_CALCULATIONS_H

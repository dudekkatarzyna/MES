//
// Created by Katarzyna on 06.10.2018.
//

#ifndef MES_UTILITY_H
#define MES_UTILITY_H

#include "../Grid/GRID.h"
#include "../UniversalElements/UniversalElement.h"
#include<fstream>

class Utility {

public:
    static void testGrid(GRID A);
    static void printNODE(GRID &, int, int);
    static void printELEMENT(GRID &, int, int);
    static void printGrid(GRID &pGRID, int, int);
    static void readFile(int*t0, int*tau, int *stepTau, int *tA, float *H, float *L, int *nH, int *nL, int *K, int *alfa, int *c, int *ro);
    static void printData(int*t0, int*tau, int *stepTau, int *tA, float *H, float *L, int *nH, int *nL, int *K, int *alfa, int *c, int *ro);
    static void printUniversalElement(UniversalElement ue);
    static void printCreateJacobian(float Jacobian[4][2][2]);
    static void printRevertJacobian(float Jacobian[4][2][2]);
    static void printJacobian(float Jacobian[4][2][2]);
    static void printDetJ(float detJ[4]);
    static void printMultiplyDetJacobian(float Jacobian[4][2][2]);
    static void printdNdXY(float divNx[4][4], float divNy[4][4]);
    static void printdivNxySqr(float divNxSqr[4][4][4], float divNySqr[4][4][4]);
    static void printMultiplyT(float divNxSqr[4][4][4], float divNySqr[4][4][4]);
    static void printRemoveIntegral(float divNxSqr[4][4][4], float divNySqr[4][4][4]);
    static void printMultiplyK(float arrK[4][4][4]);
    static void printH(float **H);
    static void printNxN(float **Arr);
    static void printSumPC(float **sum);
    static void printFinalH(float **H);
    static void printMatrix4x4(float **Arr);
    static void printNxNinC(float **MatrixCNSqrt);
    static void printMatrixC(float **MatrixC);
    static void printGlobalH(float **globalH, int nH, int nL);
    static void printGlobalC(float **globalC,  int nH, int nL);
    static void printSumBC(float **Arr);
    static void printP(float **vectorP);
    static void printGlobalP(float **globalP, int nH, int nL);
    static void printTemperature(float **t1Vector, int nH, int nL, int t);
    static void printMinMaxTemp(float min, float max);
};


#endif //MES_UTILITY_H

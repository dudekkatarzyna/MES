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
    static void readFile(double*t0, double*tau, double *stepTau, double *tA, double *H, double *L, int *nH, int *nL, double *K, double *alfa, double *c, double *ro);
    static void printData(double*t0, double*tau, double *stepTau, double *tA, double *H, double *L, int *nH, int *nL, double *K, double *alfa, double *c, double *ro);
    static void printUniversalElement(UniversalElement ue);
    static void printCreateJacobian(double Jacobian[4][2][2]);
    static void printRevertJacobian(double Jacobian[4][2][2]);
    static void printJacobian(double Jacobian[4][2][2]);
    static void printDetJ(double detJ[4]);
    static void printMultiplyDetJacobian(double Jacobian[4][2][2]);
    static void printdNdXY(double divNx[4][4], double divNy[4][4]);
    static void printdivNxySqr(double divNxSqr[4][4][4], double divNySqr[4][4][4]);
    static void printMultiplyT(double divNxSqr[4][4][4], double divNySqr[4][4][4]);
    static void printRemoveIntegral(double divNxSqr[4][4][4], double divNySqr[4][4][4]);
    static void printMultiplyK(double arrK[4][4][4]);
    static void printH(double **H);
    static void printNxN(double **Arr);
    static void printSumPC(double **sum);
    static void printFinalH(double **H);
    static void printMatrix4x4(double **Arr);
    static void printNxNinC(double **MatrixCNSqrt);
    static void printMatrixC(double **MatrixC);
    static void printGlobalH(double **globalH, int nH, int nL);
    static void printGlobalC(double **globalC,  int nH, int nL);
    static void printSumBC(double **Arr);
    static void printP(double **vectorP);
    static void printGlobalP(double **globalP, int nH, int nL);
    static void printTemperature(double **t1Vector, int nH, int nL, int t);
    static void printMinMaxTemp(double min, double max);

    static void
    readProjectFile(double *H, double *L, int *nH, int *nL, double *tau, double *stepTau, double *t0, double *tA, double *alfaPow, double *cPow,
                    double *roPow, double *kPow, double *alfa, double *c, double *ro, double *k);

    static void
    printDataExtended(double *t0, double *tau, double *stepTau, double *a, double *h, double *l, int *nH, int *nL, double *k,
                      double *alfa, double *c, double *ro, double *pow, double *kPow, double *roPow, double *alfaPow);

    static void printGridTemperature(GRID A, int nH, int nL);
};


#endif //MES_UTILITY_H

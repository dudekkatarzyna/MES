//
// Created by Katarzyna on 06.10.2018.
//

#ifndef MES_UTILITY_H
#define MES_UTILITY_H

#include "GRID.h"
#include "UniversalElement.h"
#include<fstream>

class Utility {

public:
    static void testGrid(GRID A);
    static void printNODE(GRID &, int, int);
    static void printELEMENT(GRID &, int, int);

    static void printGrid(GRID &pGRID, int, int);

    static void readFile(int *H, int *L, int *nH, int *nL, int *K, int *t);

    static void printData(int *H, int *L, int *nH, int *nL, int *K, int *t);

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

    static void printH(float H[4][4]);
};


#endif //MES_UTILITY_H

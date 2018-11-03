//
// Created by Katarzyna on 06.10.2018.
//

#ifndef MES_GRID_H
#define MES_GRID_H

#include "NODE.h"
#include "ELEMENT.h"
class KsiEta;
class GRID {

public:
    NODE *node;
    ELEMENT *element;

    GRID(int H, int L, int nH, int nL, int K, int t);

    static void createUniversalElement();

    static void createH(GRID A, int elId);
    static void createJacobian(GRID A, int elId);
    static void detJacobian(GRID A);
    static void multiplyDetJacobian(GRID A);
    static void revertJavobian(GRID A);
    static void dNdXY(GRID A);
    static void multiplyT(GRID A, float divNx[4][4], float divNy[4][4]);
    static void multiplyK(GRID A, float divNxSqr[4][4][4], float divNySqr[4][4][4]);
    static void removeIntegral(GRID A, float divNxSqr[4][4][4], float divNySqr[4][4][4]);
    static void calculateH(float arrK[4][4][4]);



};


#endif //MES_GRID_H

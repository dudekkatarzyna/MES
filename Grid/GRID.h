//
// Created by Katarzyna on 06.10.2018.
//

#ifndef MES_GRID_H
#define MES_GRID_H

#include "Components/NODE.h"
#include "Components/ELEMENT.h"
class GRID {

public:
    NODE *node;
    ELEMENT *element;

    GRID(double H, double L, int nH, int nL, double K, double t0, double alfa, double c, double ro);

    static void setBC(GRID A, int nH, int nL);

    void addSecondMaterial(int nH, int nL, double cPow, double roPow, double kPow, double alfaPow);

};


#endif //MES_GRID_H

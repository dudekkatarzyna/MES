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

    GRID(float H, float L, int nH, int nL, int K, int t0, float alfa, float c, float ro);

    static void setBC(GRID A, int nH, int nL);

};


#endif //MES_GRID_H

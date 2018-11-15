//
// Created by Katarzyna on 06.10.2018.
//

#ifndef MES_GRID_H
#define MES_GRID_H

#include "Components/NODE.h"
#include "Components/ELEMENT.h"
class KsiEta;
class GRID {

public:
    NODE *node;
    ELEMENT *element;

    GRID(float H, float L, int nH, int nL, int K, int t);

    static void setBC(GRID A);

};


#endif //MES_GRID_H

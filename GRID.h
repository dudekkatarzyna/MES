//
// Created by Katarzyna on 06.10.2018.
//

#ifndef MES_GRID_H
#define MES_GRID_H

#include "NODE.h"
#include "ELEMENT.h"

class GRID {

public:
    NODE *node;
    ELEMENT *element;

    GRID(int H, int L, int nH, int nL, int K, int t);

};


#endif //MES_GRID_H

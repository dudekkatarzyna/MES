//
// Created by Katarzyna on 06.10.2018.
//

#include "GRID.h"
#include "Utility.h"
#include<iostream>

using namespace std;

GRID::GRID(int H, int L, int nH, int nL, int K, int t) {

    node = new NODE[nH * nL];
    element = new ELEMENT[(nH - 1) * (nL - 1)];

    int i = 0;
    int deltaX = L / (nL - 1);
    int deltaY = H / (nH - 1);
    for (int l = 0; l < nL; l++) {
        for (int h = 0; h < nH; h++) {
            this->node[i].x = l * deltaX;
            this->node[i].y = h * deltaY;
            this->node[i].t = t;
            i++;
        }
    }

    int j = 0;
    for (int i = 0; i < (nH - 1) * (nL - 1); i++) {
        if (i % (nH - 1) == 0 && j != 0) {
            j++;
        }
        this->element[i].id[0] = j;
        this->element[i].id[1] = j + nH;
        this->element[i].id[2] = j + nH + 1;
        this->element[i].id[3] = j + 1;

        this->element[i].k = K;
        j++;
    }

    Utility::print(*this, nH, nL);
}




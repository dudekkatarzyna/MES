//
// Created by Katarzyna on 06.10.2018.
//

#include "GRID.h"
#include "../Utility/Utility.h"
#include "../UniversalElements/UniversalElement.h"
#include<iostream>
#include <cmath>

using namespace std;

GRID::GRID(float H, float L, int nH, int nL, int K, int t0) {

    node = new NODE[nH * nL];
    element = new ELEMENT[(nH - 1) * (nL - 1)];

    int i = 0;
    float deltaX = L / (nL - 1);
    float deltaY = H / (nH - 1);

    for (int l = 0; l < nL; l++) {
        for (int h = 0; h < nH; h++) {
            this->node[i].x = l * deltaX;
            this->node[i].y = h * deltaY;
            this->node[i].t = t0;
            i++;
        }
    }

    int j = 0;
    for (i = 0; i < (nH - 1) * (nL - 1); i++) {
        if (i % (nH - 1) == 0 && j != 0) {
            j++;
        }
        this->element[i].id[0] = j;
        this->element[i].id[1] = j + nH;
        this->element[i].id[2] = j + nH + 1;
        this->element[i].id[3] = j + 1;

        this->element[i].k = K;
        for (bool &q : this->element[i].Q) {
            q = false;
        }
        j++;
    }

    Utility::printGrid(*this, nH, nL);
}

void GRID::setBC(GRID A) {

    A.element[0].Q[0] = true;
    A.element[0].Q[3] = true;
    A.element[1].Q[3] = true;
    A.element[2].Q[2] = true;
    A.element[2].Q[3] = true;
    A.element[3].Q[0] = true;
    A.element[5].Q[2] = true;
    A.element[6].Q[0] = true;
    A.element[6].Q[1] = true;
    A.element[7].Q[1] = true;
    A.element[8].Q[1] = true;
    A.element[8].Q[2] = true;

}

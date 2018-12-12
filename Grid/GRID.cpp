//
// Created by Katarzyna on 06.10.2018.
//

#include "GRID.h"
#include "../Utility/Utility.h"
#include "../UniversalElements/UniversalElement.h"
#include<iostream>
#include <cmath>

using namespace std;

GRID::GRID(float H, float L, int nH, int nL, int K, int t0, float alfa, float c, float ro) {

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
        this->element[i].c = c;
        this->element[i].ro = ro;
        this->element[i].alfa = alfa;
        for (bool &q : this->element[i].Q) {
            q = false;
        }
        j++;
    }

    // Utility::printGrid(*this, nH, nL);
}

void GRID::setBC(GRID A, int nH, int nL) {

    bool top = true, bottom = true, right = true, left = true;


    if (bottom) {
        for (int i = 0; i < (nL - 1) * (nH - 1); i += (nH - 1)) {

            A.element[i].Q[0] = true;
        }
    }
    if (right) {
        for (int i = (nH - 1) * (nL - 2); i < (nH - 1) * (nL - 1); i++) {
            A.element[i].Q[1] = true;
        }
    }
    if (top) {
        for (int i = nH - 2; i < (nL - 1) * (nH - 1); i += (nH - 1)) {
            A.element[i].Q[2] = true;
        }
    }
    if (left) {
        for (int i = 0; i < nH - 1; i++) {
            A.element[i].Q[3] = true;
        }

    }

}

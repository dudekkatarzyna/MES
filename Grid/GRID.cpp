//
// Created by Katarzyna on 06.10.2018.
//

#include "GRID.h"
#include "../Utility/Utility.h"
#include "../UniversalElements/UniversalElement.h"
#include<iostream>
#include <cmath>

using namespace std;

GRID::GRID(double H, double L, int nH, int nL, double K, double t0, double alfa, double c, double ro) {

    node = new NODE[nH * nL];
    element = new ELEMENT[(nH - 1) * (nL - 1)];

    int i = 0;
    double deltaX = L / (nL - 1);
    double deltaY = H / (nH - 1);

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

     Utility::printGrid(*this, nH, nL);
}

void GRID::setBC(GRID A, int nH, int nL) {

    bool top = true, bottom = false, right = false, left = false;


    if (bottom) {
        for (int i = 0; i < (nL - 1) * (nH - 1); i += (nL - 1)) {

            A.element[i].Q[0] = true;
            cout << i << " ";
        }
        cout << endl;
    }
    if (right) {
        for (int i = (nH - 1) * (nL - 2); i < (nH - 1) * (nL - 1); i++) {
            A.element[i].Q[1] = true;
            cout << i << " ";
        }
        cout << endl;
    }
    if (top) {
        for (int i = nH - 2; i < (nL - 1) * (nH - 1); i += (nL - 1)) {
            A.element[i].Q[2] = true;
            cout << i << " ";
        }
        cout << endl;
    }
    if (left) {
        for (int i = 0; i < nH - 1; i++) {
            A.element[i].Q[3] = true;
            cout << i << " ";
        }
        cout << endl;
    }

}

void
GRID::addSecondMaterial(int nH, int nL,  double cPow, double roPow, double kPow, double alfaPow) {

//              2x2
//    for (int h = 2; h < 10; h = h + 6) {
//        for (int l = 26; l < 185; l = l + 112) {
//
//            cout<<l+h<<" "<<l+14+h<<" "<<l+h+1<<" "<<l+14+h+1<<endl;
//
//            this->element[l + h].k = kPow;
//            this->element[l + h].c = cPow;
//            this->element[l + h].ro = roPow;
//            this->element[l + h].alfa = alfaPow;
//
//            this->element[l + 14 + h].k = kPow;
//            this->element[l + 14 + h].c = cPow;
//            this->element[l + 14 + h].ro = roPow;
//            this->element[l + 14 + h].alfa = alfaPow;
//
//            this->element[l + h + 1].k = kPow;
//            this->element[l + h + 1].c = cPow;
//            this->element[l + h + 1].ro = roPow;
//            this->element[l + h + 1].alfa = alfaPow;
//
//            this->element[l + 14 + h + 1].k = kPow;
//            this->element[l + 14 + h + 1].c = cPow;
//            this->element[l + 14 + h + 1].ro = roPow;
//            this->element[l + 14 + h + 1].alfa = alfaPow;
//
//        }
//    }

//                  3x3 wroking
//    for (int h = 2; h < 10; h = h + 6) {
//        for (int l = 26; l < 107; l = l + 78) {
//
//            for (int k = 0; k < 3; k++) {
//                cout << l + h + k << " " << l + h + k + 13 << " " << l + h + k + 26 << endl;
//
//                this->element[l + h + k].k = kPow;
//                this->element[l + h + k].c = cPow;
//                this->element[l + h + k].ro = roPow;
//                this->element[l + h + k].alfa = alfaPow;
//
//                this->element[l + h + k + 13].k = kPow;
//                this->element[l + h + k + 13].c = cPow;
//                this->element[l + h + k + 13].ro = roPow;
//                this->element[l + h + k + 13].alfa = alfaPow;
//
//                this->element[l + h + k + 26].k = kPow;
//                this->element[l + h + k + 26].c = cPow;
//                this->element[l + h + k + 26].ro = roPow;
//                this->element[l + h + k + 26].alfa = alfaPow;
//
//
//            }
//
//        }
//        cout << endl;
//    }


//    for (int h = 2; h < 14; h = h + 6) {
//        for (int l = 28; l < 185; l = l + 14) {
//
//            if ((h + l) == 86 || (h + l) == 100 || (h + l) == 170 || (h + l) == 184 || (h + l) == 92 ||
//                (h + l) == 106 || (h + l) == 176 || (h + l) == 190)
//                continue;
//
//            for (int i = 0; i < 4; ++i) {
//                cout << l + h + i << " ";
//
//                this->element[l + h + i].k = kPow;
//                this->element[l + h + i].c = cPow;
//                this->element[l + h + i].ro = roPow;
//                this->element[l + h + i].alfa = alfaPow;
//
//            }
//            cout << endl;
//        }
//    }

        cout<<"ZZZZZZZZZZZZZZZZZZZZZZZZ"<<endl;

    for (int h = 3; h < 10; h = h + 1) {
        for (int l = 39; l < 130; l = l + 13) {

            cout << l + h << " ";

            this->element[l + h ].k = kPow;
            this->element[l + h ].c = cPow;
            this->element[l + h ].ro = roPow;
            this->element[l + h ].alfa = alfaPow;

        }
        cout << endl;
    }

    Utility::printGrid(*this, nH, nL);
}

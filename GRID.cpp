//
// Created by Katarzyna on 06.10.2018.
//

#include "GRID.h"
#include "Utility.h"
#include<iostream>
#include <cmath>

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

class KsiEta {
public:
    float ksi;
    float eta;
};

float divEta[4][4]; // dN/dEta
float divKsi[4][4]; // dN/dKsi
void GRID::createUniversalElement() {


    KsiEta elArr[4];
    elArr[0].ksi = static_cast<float>(-1 / sqrt(3));
    elArr[0].eta = static_cast<float>(-1 / sqrt(3));
    elArr[1].ksi = static_cast<float>(1 / sqrt(3));
    elArr[1].eta = static_cast<float>(-1 / sqrt(3));
    elArr[2].ksi = static_cast<float>(1 / sqrt(3));
    elArr[2].eta = static_cast<float>(1 / sqrt(3));
    elArr[3].ksi = static_cast<float>(-1 / sqrt(3));
    elArr[3].eta = static_cast<float>(1 / sqrt(3));

    for (int el = 0; el < 4; el++) {
        divEta[el][0] = static_cast<float>(-(1.0 / 4.0) * (1 - elArr[el].ksi));
        divEta[el][1] = static_cast<float>(-(1.0 / 4.0) * (1 + elArr[el].ksi));
        divEta[el][2] = static_cast<float>((1.0 / 4.0) * (1 + elArr[el].ksi));
        divEta[el][3] = static_cast<float>((1.0 / 4.0) * (1 - elArr[el].ksi));
    }

    for (int el = 0; el < 4; el++) {
        divKsi[el][0] = static_cast<float>(-(1.0 / 4.0) * (1 - elArr[el].eta));
        divKsi[el][1] = static_cast<float>((1.0 / 4.0) * (1 - elArr[el].eta));
        divKsi[el][2] = static_cast<float>((1.0 / 4.0) * (1 + elArr[el].eta));
        divKsi[el][3] = static_cast<float>(-(1.0 / 4.0) * (1 + elArr[el].eta));
    }

    for (int i = 0; i < 4; i++) {
        cout << "divKsi punkt nr " << i << endl;
        cout << elArr[i].ksi << endl;
        for (int j = 0; j < 4; j++) {
            cout << divKsi[i][j] << " ";
        }
        cout << endl;
    }

    for (int i = 0; i < 4; i++) {
        cout << "divEta punkt nr " << i << endl;
        cout << elArr[i].ksi << endl;
        for (int j = 0; j < 4; j++) {
            cout << divEta[i][j] << " ";
        }
        cout << endl;
    }

}

void GRID::createJacobian(GRID A, int nL, int nH) {

    float Jacobian[(nH - 1) * (nL - 1)][4][2][2];
    for (int elId = 0; elId < (nH - 1) * (nL - 1); elId++) {

        for (int nodeId = 0; nodeId < 4; nodeId++) {
            cout << "el: " << elId << " node: " << nodeId << ": ";
            Jacobian[elId][nodeId][0][0] =
                    divKsi[nodeId][0] * A.node[A.element[elId].id[0]].x +
                    divKsi[nodeId][1] * A.node[A.element[elId].id[1]].x +
                    divKsi[nodeId][2] * A.node[A.element[elId].id[2]].x +
                    divKsi[nodeId][3] * A.node[A.element[elId].id[3]].x;
            cout << Jacobian[elId][nodeId][0][0] << " ";
            Jacobian[elId][nodeId][0][1] =
                    divKsi[nodeId][0] * A.node[A.element[elId].id[0]].y +
                    divKsi[nodeId][1] * A.node[A.element[elId].id[1]].y +
                    divKsi[nodeId][2] * A.node[A.element[elId].id[2]].y +
                    divKsi[nodeId][3] * A.node[A.element[elId].id[3]].y;
            cout << Jacobian[elId][nodeId][0][1] << " ";
            Jacobian[elId][nodeId][1][0] =
                    divEta[nodeId][0] * A.node[A.element[elId].id[0]].x +
                    divEta[nodeId][1] * A.node[A.element[elId].id[1]].x +
                    divEta[nodeId][2] * A.node[A.element[elId].id[2]].x +
                    divEta[nodeId][3] * A.node[A.element[elId].id[3]].x;
            cout << Jacobian[elId][nodeId][1][0] << " ";
            Jacobian[elId][nodeId][1][1] =
                    divEta[nodeId][0] * A.node[A.element[elId].id[0]].y +
                    divEta[nodeId][1] * A.node[A.element[elId].id[1]].y +
                    divEta[nodeId][2] * A.node[A.element[elId].id[2]].y +
                    divEta[nodeId][3] * A.node[A.element[elId].id[3]].y;
            cout << Jacobian[elId][nodeId][1][1] << endl;

        }
    }


    //odwracanie jacobianu

    for (int elId = 0; elId < (nH - 1) * (nL - 1); elId++) {
        for (int nodeId = 0; nodeId < 4; nodeId++) {
            float tmp = Jacobian[elId][nodeId][0][0];
            Jacobian[elId][nodeId][0][0] = Jacobian[elId][nodeId][1][1];
            Jacobian[elId][nodeId][1][1] = tmp;
            Jacobian[elId][nodeId][0][1] *= (-1);
            Jacobian[elId][nodeId][1][0] *= (-1);

            cout << Jacobian[elId][nodeId][0][0] << " " << Jacobian[elId][nodeId][0][1] << " "
                 << Jacobian[elId][nodeId][1][0] << " "
                 << Jacobian[elId][nodeId][1][1] << endl;
        }

    }

    //obliczanie wyznacznika
    float detJ[(nH - 1) * (nL - 1)][4];
    for (int elId = 0; elId < (nH - 1) * (nL - 1); elId++) {
        for (int nodeId = 0; nodeId < 4; nodeId++) {

            detJ[elId][nodeId] = Jacobian[elId][nodeId][0][0] * Jacobian[elId][nodeId][1][1] -
                                 Jacobian[elId][nodeId][0][1] * Jacobian[elId][nodeId][1][0];
            cout << "elId:" << elId << " nodeId:" << nodeId << " detJ:" << detJ[elId][nodeId] << endl;

        }
    }

    //wymnożenie Jacobianu przez 1/detJ
    for (int elId = 0; elId < (nH - 1) * (nL - 1); elId++) {
        for (int nodeId = 0; nodeId < 4; nodeId++) {
            Jacobian[elId][nodeId][0][0] *= (1 / detJ[elId][nodeId]);
            Jacobian[elId][nodeId][0][1] *= (1 / detJ[elId][nodeId]);
            Jacobian[elId][nodeId][1][0] *= (1 / detJ[elId][nodeId]);
            Jacobian[elId][nodeId][1][1] *= (1 / detJ[elId][nodeId]);

            cout << Jacobian[elId][nodeId][0][0] << " " << Jacobian[elId][nodeId][0][1] << " "
                 << Jacobian[elId][nodeId][1][0] << " "
                 << Jacobian[elId][nodeId][1][1] << endl;
        }
    }
    // wymnożenie do dN/dx i dN/dy

    float divNx[4][4];
    float divNy[4][4];

    for (int elId = 0; elId < (nH - 1) * (nL - 1); elId++) {
        for (int nodeId = 0; nodeId < 4; nodeId++) {

            divNx[nodeId][0] =
                    Jacobian[elId][nodeId][0][0] * divKsi[0][0] + Jacobian[elId][nodeId][0][1] * divEta[nodeId][0];
            divNx[nodeId][1] =
                    Jacobian[elId][nodeId][0][0] * divKsi[0][1] + Jacobian[elId][nodeId][0][1] * divEta[nodeId][1];
            divNx[nodeId][2] =
                    Jacobian[elId][nodeId][0][0] * divKsi[0][2] + Jacobian[elId][nodeId][0][1] * divEta[nodeId][2];
            divNx[nodeId][3] =
                    Jacobian[elId][nodeId][0][0] * divKsi[0][3] + Jacobian[elId][nodeId][0][1] * divEta[nodeId][3];

            divNy[nodeId][0] =
                    Jacobian[elId][nodeId][1][0] * divKsi[0][0] + Jacobian[elId][nodeId][1][1] * divEta[nodeId][0];
            divNy[nodeId][1] =
                    Jacobian[elId][nodeId][1][0] * divKsi[0][1] + Jacobian[elId][nodeId][1][1] * divEta[nodeId][1];
            divNy[nodeId][2] =
                    Jacobian[elId][nodeId][1][0] * divKsi[0][2] + Jacobian[elId][nodeId][1][1] * divEta[nodeId][2];
            divNy[nodeId][3] =
                    Jacobian[elId][nodeId][1][0] * divKsi[0][3] + Jacobian[elId][nodeId][1][1] * divEta[nodeId][3];
        }
    }
    cout << divNx[0][0] << " " << divNx[0][1] << " " << divNx[0][2] << " " << divNx[0][3] << endl;
    cout << divNx[1][0] << " " << divNx[1][1] << " " << divNx[1][2] << " " << divNx[1][3] << endl;
    cout << divNx[2][0] << " " << divNx[2][1] << " " << divNx[2][2] << " " << divNx[2][3] << endl;
    cout << divNx[3][0] << " " << divNx[3][1] << " " << divNx[3][2] << " " << divNx[3][3] << endl;
    cout << endl;
    cout << divNy[0][0] << " " << divNy[0][1] << " " << divNy[0][2] << " " << divNy[0][3] << endl;
    cout << divNy[1][0] << " " << divNy[1][1] << " " << divNy[1][2] << " " << divNy[1][3] << endl;
    cout << divNy[2][0] << " " << divNy[2][1] << " " << divNy[2][2] << " " << divNy[2][3] << endl;
    cout << divNy[3][0] << " " << divNy[3][1] << " " << divNy[3][2] << " " << divNy[3][3] << endl;

    //dN/dx * (dN/dx)T  mnożenie razy transpozycja TODO: błędne obliczanie pkt 3 i 4

    cout << endl << endl;
    float divNxSqr[8][4][4];
    float divNySqr[8][4][4];

    for (int arr = 0; arr < 8; arr++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                divNxSqr[arr][i][j] = divNx[arr][i] * divNx[arr][j];

                divNySqr[arr][i][j] = divNy[arr][i] * divNy[arr][j];
            }
        }
    }
    cout<<"MNOZENIE PRZEZ TRANSPOZYCJE PKT 3"<<endl;
    cout << divNxSqr[2][0][0] << " " << divNxSqr[2][0][1] << " " << divNxSqr[2][0][2] << " " << divNxSqr[2][0][3]
         << endl;
    cout << divNxSqr[2][1][0] << " " << divNxSqr[2][1][1] << " " << divNxSqr[2][1][2] << " " << divNxSqr[2][1][3]
         << endl;
    cout << divNxSqr[2][2][0] << " " << divNxSqr[2][2][1] << " " << divNxSqr[2][2][2] << " " << divNxSqr[2][2][3]
         << endl;
    cout << divNxSqr[2][3][0] << " " << divNxSqr[2][3][1] << " " << divNxSqr[2][3][2] << " " << divNxSqr[2][3][3]
         << endl;
    cout << endl;
    cout << divNySqr[0][0][0] << " " << divNySqr[0][0][1] << " " << divNySqr[0][0][2] << " " << divNySqr[0][0][3]
         << endl;
    cout << divNySqr[0][1][0] << " " << divNySqr[0][1][1] << " " << divNySqr[0][1][2] << " " << divNySqr[0][1][3]
         << endl;
    cout << divNySqr[0][2][0] << " " << divNySqr[0][2][1] << " " << divNySqr[0][2][2] << " " << divNySqr[0][2][3]
         << endl;
    cout << divNySqr[0][3][0] << " " << divNySqr[0][3][1] << " " << divNySqr[0][3][2] << " " << divNySqr[0][3][3]
         << endl;

    //każda macierz przez swój wyznacznik

    for (int arr = 0; arr < 4; arr++) {
        float detNx = divNxSqr[arr][0][0] * divNxSqr[arr][1][1] * divNxSqr[arr][2][2] * divNxSqr[arr][3][3] +
                      divNxSqr[arr][1][0] * divNxSqr[arr][2][1] * divNxSqr[arr][3][2] * divNxSqr[arr][0][4] +
                      divNxSqr[arr][2][0] * divNxSqr[arr][3][1] * divNxSqr[arr][0][2] * divNxSqr[arr][1][3] +
                      divNxSqr[arr][3][0] * divNxSqr[arr][0][1] * divNxSqr[arr][1][2] * divNxSqr[arr][2][3] -
                      divNxSqr[arr][0][3] * divNxSqr[arr][1][2] * divNxSqr[arr][2][1] * divNxSqr[arr][3][0] -
                      divNxSqr[arr][1][3] * divNxSqr[arr][2][2] * divNxSqr[arr][3][1] * divNxSqr[arr][0][0] -
                      divNxSqr[arr][2][3] * divNxSqr[arr][3][2] * divNxSqr[arr][0][1] * divNxSqr[arr][1][0] -
                      divNxSqr[arr][3][3] * divNxSqr[arr][0][2] * divNxSqr[arr][1][1] * divNxSqr[arr][2][0];
        float detNy = divNySqr[arr][0][0] * divNySqr[arr][1][1] * divNySqr[arr][2][2] * divNySqr[arr][3][3] +
                      divNySqr[arr][1][0] * divNySqr[arr][2][1] * divNySqr[arr][3][2] * divNySqr[arr][0][4] +
                      divNySqr[arr][2][0] * divNySqr[arr][3][1] * divNySqr[arr][0][2] * divNySqr[arr][1][3] +
                      divNySqr[arr][3][0] * divNySqr[arr][0][1] * divNySqr[arr][1][2] * divNySqr[arr][2][3] -
                      divNySqr[arr][0][3] * divNySqr[arr][1][2] * divNySqr[arr][2][1] * divNySqr[arr][3][0] -
                      divNySqr[arr][1][3] * divNySqr[arr][2][2] * divNySqr[arr][3][1] * divNySqr[arr][0][0] -
                      divNySqr[arr][2][3] * divNySqr[arr][3][2] * divNySqr[arr][0][1] * divNySqr[arr][1][0] -
                      divNySqr[arr][3][3] * divNySqr[arr][0][2] * divNySqr[arr][1][1] * divNySqr[arr][2][0];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {


                divNxSqr[arr][i][j] *= detJ[0][0]; //swój czy jacobianu?

                divNySqr[arr][i][j] *= detJ[0][0];

                cout<<endl<<divNxSqr[arr][i][j]<<" ";
            }
            cout<<endl;
        }
    }
    cout << endl;
    cout << divNxSqr[0][0][0] << " " << divNxSqr[0][0][1] << " " << divNxSqr[0][0][2] << " " << divNxSqr[0][0][3]
         << endl;
    cout << divNxSqr[0][1][0] << " " << divNxSqr[0][1][1] << " " << divNxSqr[0][1][2] << " " << divNxSqr[0][1][3]
         << endl;
    cout << divNxSqr[0][2][0] << " " << divNxSqr[0][2][1] << " " << divNxSqr[0][2][2] << " " << divNxSqr[0][2][3]
         << endl;
    cout << divNxSqr[0][3][0] << " " << divNxSqr[0][3][1] << " " << divNxSqr[0][3][2] << " " << divNxSqr[0][3][3]
         << endl;
    cout << endl;
    cout << divNySqr[0][0][0] << " " << divNySqr[0][0][1] << " " << divNySqr[0][0][2] << " " << divNySqr[0][0][3]
         << endl;
    cout << divNySqr[0][1][0] << " " << divNySqr[0][1][1] << " " << divNySqr[0][1][2] << " " << divNySqr[0][1][3]
         << endl;
    cout << divNySqr[0][2][0] << " " << divNySqr[0][2][1] << " " << divNySqr[0][2][2] << " " << divNySqr[0][2][3]
         << endl;
    cout << divNySqr[0][3][0] << " " << divNySqr[0][3][1] << " " << divNySqr[0][3][2] << " " << divNySqr[0][3][3]
         << endl;

    //mnożenie razy współczynnik k

    cout << endl;
    float arrK[4][4][4];
    for (int pkt = 0; pkt < 4; pkt++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                arrK[pkt][i][j] = A.element[pkt].k * (divNxSqr[pkt][i][j] + divNySqr[pkt][i][j]);
                cout << arrK[pkt][i][j] << " ";
            }
            cout << endl;
        }
        cout << endl << endl;
    }

    //macierz H

    float H[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            H[i][j]=0;
            for(int arr=0;arr<4;arr++) {
                H[i][j] += arrK[arr][i][j];
            }
            cout<<H[i][j]<<" ";
        }
        cout<<endl;
    }
}


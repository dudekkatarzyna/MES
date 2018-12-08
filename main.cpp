#include <iostream>
#include<fstream>
#include <cmath>
#include"Grid/GRID.h"
#include "Utility/Utility.h"
#include "UniversalElements/UniversalElementSides.h"
#include "Utility/Calculations.h"

using namespace std;

int main() {
    float H, L;
    int nH, nL, K, alfa, c, ro, t0, tau, stepTau, tA;

    Utility::readFile(&t0, &tau, &stepTau, &tA, &H, &L, &nH, &nL, &K, &alfa, &c, &ro);
    GRID A(H, L, nH, nL, K, t0);
    //Utility::testGrid(A);

    GRID::setBC(A);

    float **globalH, **globalC, **globalP;
    float **matrixH, **matrixC, **vectorP;
    float **bordCond;
    float **t1Vector;

    globalH = new float *[nL * nH];
    globalC = new float *[nL * nH];
    globalP = new float *[nL * nH];
    t1Vector = new float *[nL * nH];


    for (int k = 0; k < nL * nH; ++k) {
        globalH[k] = new float[nL * nH]();
        globalC[k] = new float[nL * nH]();
        globalP[k] = new float[nL * nH]();
        t1Vector[k] = new float[nL * nH]();
    }

    for (int t = 0; t < tau; t += stepTau) {



        //clean arrays from next iteration
        for (int k = 0; k < nH * nL; ++k) {
            for (int i = 0; i < nH * nL; ++i) {
                globalH[k][i] = 0;
                globalC[k][i] = 0;
                globalP[k][i] = 0;
                t1Vector[k][i] = 0;
            }
        }

        for (int el = 0; el < ((nL - 1) * (nH - 1)); el++) {

            matrixH = Calculations::createH(A, el);
            bordCond = Calculations::addBordCondition(A, el, alfa);
            Calculations::sumArrays(matrixH, bordCond);
            matrixC = Calculations::matrixC(c, ro);
            vectorP = Calculations::vectorP(A, el, alfa, tA);

            //agregation to nH*nL x nH*nL array
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    globalH[A.element[el].id[i]][A.element[el].id[j]] += matrixH[i][j];
                    globalC[A.element[el].id[i]][A.element[el].id[j]] += matrixC[i][j];
                }
                globalP[A.element[el].id[i]][0] += vectorP[i][0];
            }

        }

        // [C]/dTau
        for (int i = 0; i < nL * nH; i++) {
            for (int j = 0; j < nL * nH; j++) {
                globalC[i][j] /= stepTau;
            }
        }

        //[H] + [C]/dTau
        for (int i = 0; i < nL * nH; i++) {
            for (int j = 0; j < nL * nH; j++) {
                globalH[i][j] += globalC[i][j];
            }
        }

        // {P} + [C]/dTau * t0
        for (int i = 0; i < nL * nH; i++) {
            for (int j = 0; j < nL * nH; j++) {
                globalP[i][0] += (globalC[i][j] * A.node[j].t);
            }
        }

        //obliczenie równania  t1= [H]^{-1}*{P}
        t1Vector = Calculations::solveEqationForT(globalH, globalP, nL, nH, t1Vector);
        Utility::printTemperature(t1Vector, nH, nL, t);
        Calculations::getMinMaxtemp(t1Vector, nH, nL);

        //switch temperatures for next iteration
        for (int i = 0; i < nL * nH; i++) {
            A.node[i].t = t1Vector[i][0];
        }
    }
}
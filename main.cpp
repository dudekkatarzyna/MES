#include <iostream>
#include<fstream>
#include <cmath>
#include"Grid/GRID.h"
#include "Utility/Utility.h"
#include "UniversalElements/UniversalElementSides.h"
#include "Utility/Calculations.h"

using namespace std;

int main() {
    double H, L, k, alfa, c, ro, t0, tau, stepTau, tA;
    int nH, nL;

    double alfaPow, cPow, roPow, kPow;
    Utility::readProjectFile(&H, &L, &nH, &nL, &tau, &stepTau, &t0, &tA, &alfaPow, &cPow, &roPow, &kPow, &alfa, &c, &ro,
                             &k);

    //  Utility::readFile(&t0, &tau, &stepTau, &tA, &H, &L, &nH, &nL, &K, &alfa, &c, &ro);

    GRID A(H, L, nH, nL, k, t0, alfa, c, ro);
    //Utility::testGrid(A);

    A.addSecondMaterial(nH, nL, cPow, roPow, kPow, alfaPow);

    //  Utility::printGrid(A, nH, nL);

    GRID::setBC(A, nH, nL);

    double **globalH, **globalC, **globalP;
    double **matrixH, **matrixC, **vectorP;
    double **bordCond;
    double **t1Vector;

    globalH = new double *[nL * nH];
    globalC = new double *[nL * nH];
    globalP = new double *[nL * nH];
    t1Vector = new double *[nL * nH];


    for (int m = 0; m < nL * nH; ++m) {
        globalH[m] = new double[nL * nH]();
        globalC[m] = new double[nL * nH]();
        globalP[m] = new double[nL * nH]();
        t1Vector[m] = new double[nL * nH]();
    }

    for (int t = 0; t <= tau; t += stepTau) {

        //clean arrays from prev iteration
        for (int m = 0; m < nH * nL; ++m) {
            for (int i = 0; i < nH * nL; ++i) {
                globalH[m][i] = 0;
                globalC[m][i] = 0;
                globalP[m][i] = 0;
                t1Vector[m][i] = 0;
            }
        }

        for (int el = 0; el < ((nL - 1) * (nH - 1)); el++) {

            matrixH = Calculations::createH(A, el);
            bordCond = Calculations::addBordCondition(A, el, A.element[el].alfa);
            Calculations::sumArrays(matrixH, bordCond);
            matrixC = Calculations::matrixC(A.element[el].c, A.element[el].ro);
            vectorP = Calculations::vectorP(A, el, A.element[el].alfa, tA);

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

        //   Utility::printGlobalH(globalH,nH,nL);

        t1Vector = Calculations::solveEqationForT(globalH, globalP, nL, nH, t1Vector);
        // Utility::printTemperature(t1Vector, nH, nL, t);
        Calculations::getMinMaxtemp(t1Vector, nH, nL, t);

        //switch temperatures for next iteration
        for (int i = 0; i < nL * nH; i++) {
            A.node[i].t = t1Vector[i][0];
        }

        //Utility::printGridTemperature(A, nH, nL);
    }

}
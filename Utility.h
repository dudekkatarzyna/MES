//
// Created by Katarzyna on 06.10.2018.
//

#ifndef MES_UTILITY_H
#define MES_UTILITY_H

#include "GRID.h"
#include<fstream>

class Utility {

public:
    static void testGrid(GRID A);
    static void printNODE(GRID &, int, int);
    static void printELEMENT(GRID &, int, int);

    static void print(GRID &pGRID, int, int);

    static void readFile(int *H, int *L, int *nH, int *nL, int *K, int *t);

    static void printData(int *H, int *L, int *nH, int *nL, int *K, int *t);
};


#endif //MES_UTILITY_H

//
// Created by Katarzyna on 12.11.2018.
//

#include <cmath>
#include <iostream>
#include "UniversalElementSides.h"

UniversalElementSides::UniversalElementSides() {


    side[0][0].ksi = -1 / sqrt(3);
    side[0][0].eta = -1;
    side[0][1].ksi = 1 / sqrt(3);
    side[0][1].eta = -1;
    side[1][0].ksi = 1;
    side[1][0].eta = -1 / sqrt(3);
    side[1][1].ksi = 1;
    side[1][1].eta = 1 / sqrt(3);
    side[2][0].ksi = 1 / sqrt(3);
    side[2][0].eta = 1;
    side[2][1].ksi = -1 / sqrt(3);
    side[2][1].eta = 1;
    side[3][0].ksi = -1;
    side[3][0].eta = 1 / sqrt(3);
    side[3][1].ksi = -1;
    side[3][1].eta = -1 / sqrt(3);

}

void UniversalElementSides::N(double PC1[4], double PC2[4], int s) {

    PC1[0] = 0.25 * (1 - side[s][0].ksi) * (1 - side[s][0].eta);
    PC1[1] = 0.25 * (1 + side[s][0].ksi) * (1 - side[s][0].eta);
    PC1[2] = 0.25 * (1 + side[s][0].ksi) * (1 + side[s][0].eta);
    PC1[3] = 0.25 * (1 - side[s][0].ksi) * (1 + side[s][0].eta);

    PC2[0] = 0.25 * (1 - side[s][1].ksi) * (1 - side[s][1].eta);
    PC2[1] = 0.25 * (1 + side[s][1].ksi) * (1 - side[s][1].eta);
    PC2[2] = 0.25 * (1 + side[s][1].ksi) * (1 + side[s][1].eta);
    PC2[3] = 0.25 * (1 - side[s][1].ksi) * (1 + side[s][1].eta);


}

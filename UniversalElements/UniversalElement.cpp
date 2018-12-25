//
// Created by Katarzyna on 08.11.2018.
//

#include <cmath>
#include "UniversalElement.h"
#include "../Utility/Utility.h"
#include "KsiEta.h"

UniversalElement::UniversalElement() {

    KsiEta elArr[4];
    elArr[0].ksi = -1 / sqrt(3);
    elArr[0].eta = -1 / sqrt(3);
    elArr[1].ksi = 1 / sqrt(3);
    elArr[1].eta = -1 / sqrt(3);
    elArr[2].ksi = 1 / sqrt(3);
    elArr[2].eta = 1 / sqrt(3);
    elArr[3].ksi = -1 / sqrt(3);
    elArr[3].eta = 1 / sqrt(3);

    for (int el = 0; el < 4; el++) {
        divEta[el][0] = -(1.0 / 4.0) * (1 - elArr[el].ksi);
        divEta[el][1] = -(1.0 / 4.0) * (1 + elArr[el].ksi);
        divEta[el][2] = (1.0 / 4.0) * (1 + elArr[el].ksi);
        divEta[el][3] = (1.0 / 4.0) * (1 - elArr[el].ksi);
    }

    for (int el = 0; el < 4; el++) {
        divKsi[el][0] = -(1.0 / 4.0) * (1 - elArr[el].eta);
        divKsi[el][1] = (1.0 / 4.0) * (1 - elArr[el].eta);
        divKsi[el][2] = (1.0 / 4.0) * (1 + elArr[el].eta);
        divKsi[el][3] = -(1.0 / 4.0) * (1 + elArr[el].eta);
    }

}

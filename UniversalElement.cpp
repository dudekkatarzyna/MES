//
// Created by Katarzyna on 08.11.2018.
//

#include <cmath>
#include "UniversalElement.h"
#include "Utility.h"

UniversalElement::UniversalElement() {


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


}

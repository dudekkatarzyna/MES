//
// Created by Katarzyna on 08.11.2018.
//

#ifndef MES_UNIVERSALELEMENT_H
#define MES_UNIVERSALELEMENT_H


class UniversalElement {
public:
    UniversalElement();

private:

public:
    double divEta[4][4]; // dN/dEta
    double divKsi[4][4]; // dN/dKsi
};


#endif //MES_UNIVERSALELEMENT_H

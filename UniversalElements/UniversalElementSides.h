//
// Created by Katarzyna on 12.11.2018.
//

#ifndef MES_UNIVERSALELEMENTSIDES_H
#define MES_UNIVERSALELEMENTSIDES_H


class UniversalElementSides {

public:
    UniversalElementSides();

private:

    class KsiEta {
    public:
        double ksi;
        double eta;
    };

public:
    KsiEta side[4][2];
    void N(double PC1[4], double PC2[4], int el);

};


#endif //MES_UNIVERSALELEMENTSIDES_H

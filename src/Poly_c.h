//
// Created by kore on 2/25/22.
//

#ifndef KNN_SERVER_POLY_C_H
#define KNN_SERVER_POLY_C_H
#include "gmpxx.h"

class Poly_c {
public:
    int num;
    mpz_t co[];

    Poly_c(mpz_t coe[], int num_out){
        this->num = num_out;
        for(int i=0; i<num; i++){
            mpz_init_set(this->co[i], coe[i]);
        }
    }

    int evaluate(mpz_t result ,mpz_t x);

};


#endif //KNN_SERVER_POLY_C_H

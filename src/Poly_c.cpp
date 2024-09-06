//
// Created by kore on 2/25/22.
//

#include "Poly_c.h"

int Poly_c::evaluate(mpz_t result ,mpz_t x){
    mpz_t sum, tmp;
    mpz_init(tmp);
    mpz_init_set_ui(sum, 0);
    for(int i=0; i< this->num; i++){
        mpz_pow_ui(tmp, x, i);
        mpz_mul(tmp, this->co[i], tmp);
        mpz_add(sum, sum, tmp);
    }
    mpz_set(result, sum);

    mpz_clear(sum);
    mpz_clear(tmp);

    return 0;
}
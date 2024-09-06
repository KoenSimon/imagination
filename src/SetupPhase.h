//
// Created by kore on 9/8/21.
//

#ifndef kNN_Client_SETUPPHASE_H
#define kNN_Client_SETUPPHASE_H


#include "crypto.h"

class SetupPhase {
public:
    paillier_public_key pks;
    hss_pk_p hpk;
    hss_ek_p hek0, hek1;
    hss_men_p mof11;
    hss_men_p mof12;

    int KeySize;

    mpz_t base_sk, base_msg;


    SetupPhase(){
        int result = 0;

        mpz_t randTmp;
        mpz_init(randTmp);

        mpz_init(pks.n);
        hpk.pk = &pks;

        //int KeySize;
        mpz_init(base_sk);
        mpz_init(base_msg);
        setup_PKI_out(&hpk, &hek0, &hek1, &KeySize, base_sk, base_msg);

        mpz_init(mof11.sharex);
        mpz_init(mof12.sharex);

        gmp_randstate_t rstate;
        gmp_randinit_mt(rstate);

        mpz_urandomm(randTmp, rstate, (hpk.pk)->n);
        mpz_set(mof11.sharex, randTmp);
        mpz_add_ui(randTmp, randTmp, 1);
        mpz_set(mof12.sharex, randTmp);


        for (int i = 0; i < KeySize; i++) {
            //cout << hek0.ds[i] << endl;
            mpz_init_set(mof11.sharedx[i], hek0.ds[i]);
            mpz_init_set(mof12.sharedx[i], hek1.ds[i]);
        }


    }

    SetupPhase(char const * file_path, int server_id){
        int result = 0;

        mpz_t randTmp;
        mpz_init(randTmp);

        mpz_init(pks.n);
        hpk.pk = &pks;

        mpz_init(base_sk);
        mpz_init(base_msg);
        setup_PKI_out(&hpk, &hek0, &hek1, &KeySize, base_sk, base_msg);
        std::cout << "simulated initialize complete." << std::endl;
        //Not suppose to do this, but for convinient.

        mpz_init(mof11.sharex);
        mpz_init(mof12.sharex);
/*
        gmp_randstate_t rstate;
        gmp_randinit_mt(rstate);

        mpz_urandomm(randTmp, rstate, (hpk.pk)->n);
        mpz_set(mof11.sharex, randTmp);
        mpz_add_ui(randTmp, randTmp, 1);
        mpz_set(mof12.sharex, randTmp);
*/

        for (int i = 0; i < KeySize; i++) {
            //cout << hek0.ds[i] << endl;
            mpz_init(mof11.sharedx[i]);
            mpz_init(mof12.sharedx[i]);
        }
        importKeys(file_path, server_id);
        //changeKeySize(KeySize);
    }

    int exprotKeys(char const * file_path);
    int importKeys(char const * file_path, int server_id);
};



#endif //kNN_Client_SETUPPHASE_H

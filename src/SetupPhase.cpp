//
// Created by kore on 9/8/21.
//

#include "SetupPhase.h"
using namespace std;

int SetupPhase::exprotKeys(char const * file_path) {
    unsigned int result_flag, sum_out = 0;

    FILE * fp_out = NULL;
    fp_out = fopen(file_path, "w+");

    cout << "Begin!" << endl;

    //    paillier_public_key pks;
    mpz_t sizeOut;
    mpz_init(sizeOut);
    mpz_set_ui(sizeOut, KeySize);
    result_flag = mpz_out_raw(fp_out, sizeOut);
    if(result_flag == 0){
        cerr << "Unable to write!" << endl;
        return -1;
    }else{
        sum_out += result_flag;
    }


    int len_of_pk = this->pks.len;
    mpz_t mpz_len_of_pk;
    mpz_init(mpz_len_of_pk);
    mpz_set_ui(mpz_len_of_pk, len_of_pk);

    
    result_flag = mpz_out_raw(fp_out, mpz_len_of_pk);
    if(result_flag == 0){
        cerr << "Unable to write!" << endl;
        return -1;
    }else{
        sum_out += result_flag;
    }

    result_flag = mpz_out_raw(fp_out, this->pks.n);
    if(result_flag == 0){
        cerr << "Unable to write!" << endl;
        return -1;
    }else{
        sum_out += result_flag;
    }


    //    hss_pk_p hpk;
    //    omit the pk

    for(int keyind=0; keyind < KeySize; keyind++){
        result_flag = mpz_out_raw(fp_out, this->hpk.Ds[keyind]);
        if(result_flag == 0){
            cerr << "Unable to write!" << endl;
            return -1;
        }else{
            sum_out += result_flag;
        }
    }

    cout << "hpk complete!" << endl;


    //    hss_ek_p hek0, hek1;
    result_flag = mpz_out_raw(fp_out, this->hek0.kprf);
    if(result_flag == 0){
        cerr << "Unable to write!" << endl;
        return -1;
    }else{
        sum_out += result_flag;
    }




    for(int keyind=0; keyind < KeySize; keyind++){
        result_flag = mpz_out_raw(fp_out, this->hek0.ds[keyind]);
        if(result_flag == 0){
            cerr << "Unable to write!" << endl;
            return -1;
        }else{
            sum_out += result_flag;
        }
    }

    cout << "hek0 complete!" << endl;

    result_flag = mpz_out_raw(fp_out, this->hek1.kprf);
    if(result_flag == 0){
        cerr << "Unable to write!" << endl;
        return -1;
    }else{
        sum_out += result_flag;
    }


    for(int keyind=0; keyind < KeySize; keyind++){
        result_flag = mpz_out_raw(fp_out, this->hek1.ds[keyind]);
        if(result_flag == 0){
            cerr << "Unable to write!" << endl;
            return -1;
        }else{
            sum_out += result_flag;
        }
    }

    cout << "hek1 complete!" << endl;

    //    hss_men_p mof11;
    result_flag = mpz_out_raw(fp_out, this->mof11.sharex);
    if(result_flag == 0){
        cerr << "Unable to write!" << endl;
        return -1;
    }else{
        sum_out += result_flag;
    }


    for(int keyind=0; keyind < KeySize; keyind++){
        result_flag = mpz_out_raw(fp_out, this->mof11.sharedx[keyind]);
        if(result_flag == 0){
            cerr << "Unable to write!" << endl;
            return -1;
        }else{
            sum_out += result_flag;
        }
    }


    //    hss_men_p mof12;
    result_flag = mpz_out_raw(fp_out, this->mof12.sharex);
    if(result_flag == 0){
        cerr << "Unable to write!" << endl;
        return -1;
    }else{
        sum_out += result_flag;
    }


    for(int keyind=0; keyind < KeySize; keyind++){
        result_flag = mpz_out_raw(fp_out, this->mof12.sharedx[keyind]);
        if(result_flag == 0){
            cerr << "Unable to write!" << endl;
            return -1;
        }else{
            sum_out += result_flag;
        }
    }

    cout << "Result: Write " << sum_out << " bytes." << endl;
    fclose(fp_out);

    return 0;
}

int SetupPhase::importKeys(char const * file_path, int server_id) {
    unsigned int result_flag, sum_in = 0;

    FILE *fp_in = NULL;
    fp_in = fopen(file_path, "r");

    mpz_t tmpNum;
    mpz_init(tmpNum);


    mpz_t inSize;
    mpz_init(inSize);

    result_flag = mpz_inp_raw(inSize, fp_in);
    if(result_flag == 0){
        cerr << "Unable to read!" << endl;
        return -1;
    }else{
        sum_in += result_flag;
    }
    KeySize = mpz_get_ui(inSize);
    cout << "Import keysize:" << KeySize << endl;


    //    paillier_public_key pks;

    mpz_t mpz_len_of_pk;
    mpz_init(mpz_len_of_pk);

    result_flag = mpz_inp_raw(tmpNum, fp_in);
    if(result_flag == 0){
        cerr << "Unable to read!" << endl;
        return -1;
    }else{
        sum_in += result_flag;
    }

    mpz_set(mpz_len_of_pk, tmpNum);
    //int len_of_pk = this->pks.len;
    this->pks.len = mpz_get_ui(mpz_len_of_pk);

    result_flag = mpz_inp_raw(this->pks.n, fp_in);
    if(result_flag == 0){
        cerr << "Unable to read!" << endl;
        return -1;
    }else{
        sum_in += result_flag;
    }

    cout << "report n is " << this->pks.n << endl;

    //    hss_pk_p hpk;
    //    omit the pk

    for(int keyind=0; keyind < KeySize; keyind++){
        result_flag = mpz_inp_raw(this->hpk.Ds[keyind], fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }
    }

    if(server_id == 0){
        result_flag = mpz_inp_raw(this->hek0.kprf, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }



        for(int keyind=0; keyind < KeySize; keyind++){
            result_flag = mpz_inp_raw(this->hek0.ds[keyind], fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }
        }


        result_flag = mpz_inp_raw(this->hek1.kprf, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }


        for(int keyind=0; keyind < KeySize; keyind++){
            result_flag = mpz_inp_raw(this->hek1.ds[keyind], fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }
        }

        //    hss_men_p mof11;
        result_flag = mpz_inp_raw(this->mof11.sharex, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }


        for(int keyind=0; keyind < KeySize; keyind++){
            result_flag = mpz_inp_raw(this->mof11.sharedx[keyind], fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }
        }


        //    hss_men_p mof12;
        result_flag = mpz_inp_raw(this->mof12.sharex, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }


        for(int keyind=0; keyind < KeySize; keyind++){
            result_flag = mpz_inp_raw(this->mof12.sharedx[keyind], fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }
        }

    }else{
        result_flag = mpz_inp_raw(this->hek1.kprf, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }



        for(int keyind=0; keyind < KeySize; keyind++){
            result_flag = mpz_inp_raw(this->hek1.ds[keyind], fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }
        }


        result_flag = mpz_inp_raw(this->hek0.kprf, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }


        for(int keyind=0; keyind < KeySize; keyind++){
            result_flag = mpz_inp_raw(this->hek0.ds[keyind], fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }
        }

        //    hss_men_p mof11;
        result_flag = mpz_inp_raw(this->mof12.sharex, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }


        for(int keyind=0; keyind < KeySize; keyind++){
            result_flag = mpz_inp_raw(this->mof12.sharedx[keyind], fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }
        }


        //    hss_men_p mof12;
        result_flag = mpz_inp_raw(this->mof11.sharex, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }


        for(int keyind=0; keyind < KeySize; keyind++){
            result_flag = mpz_inp_raw(this->mof11.sharedx[keyind], fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }
        }
    }




    cout << "Result: Read " << sum_in << " bytes." << endl;
    fclose(fp_in);

    return 0;
}
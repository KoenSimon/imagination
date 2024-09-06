//
// Created by kore on 9/12/21.
//

#ifndef PPKNN_KNN_SINGLE_H
#define PPKNN_KNN_SINGLE_H


#include <iostream>
#include <vector>
#include <algorithm>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "SetupPhase.h"
#include<sys/time.h>







class KNN_single {
public:

    int server_id;

    //SetupPhase* setup;
    hss_men_p * mof11;
    hss_pk_p *hpk;
    hss_ek_p* hek0;

    int * KeySize;
    mpz_t base_sk;

    mpz_t * X2;

    mpz_t Qs[pthread_max_number][attribute_number];
    mpz_t Q2s[pthread_max_number];
    mpz_t randoms[pthread_max_number];

    mpz_t refreshed_label[pthread_max_number];
    mpz_t final_label[pthread_max_number];

    //std::vector<std::vector<hss_men_p> > Xi1M;
    std::array<std::array<hss_men_p, attribute_number> , train_number > * Xi1M;// = new std::array<std::array<hss_input_p, attribute_number> , train_number >(); //input

    std::array<std::array<mpz_t, attribute_number> ,train_number > *  Xi1Md;

    std::array<mpz_t, train_number> * Yw1;

    mpz_t * distance_perm;

    struct node_paillier * qX1;

    int sind[train_number];

    KNN_single(//std::array<std::array<hss_input_p, attribute_number> , train_number >* training_data_s,   std::array<mpz_t, train_number> * input_Yw1,
               hss_men_p * input_Mof1, hss_pk_p *input_hpk, hss_ek_p* input_ek, TrainingParams params,
               int * input_KeySize, mpz_t input_base_sk, int server_id){
        //hss_men_p* Mof1, hss_input_p* Ix, hss_ek_p* ek, paillier_public_key* pub, int serverId
        struct timeval t1,t2;
        double timeuse;
        gettimeofday(&t1,NULL);

        this->server_id = server_id;

        std::array<std::array<hss_input_p, attribute_number> , train_number > * Xi1 = new std::array<std::array<hss_input_p, attribute_number> , train_number >(); //input

        this->mof11 = input_Mof1;
        this->hpk = input_hpk;
        this->hek0 = input_ek;
        this->qX1 = new node_paillier[online_train_number+100];

        X2 = new mpz_t[train_number];
        distance_perm = new mpz_t [const_num];


        this->KeySize = input_KeySize;
        std::cout << "Setup key size: " << (*this->KeySize)  << std::endl;
        mpz_init_set(base_sk, input_base_sk);


        std::cout << "Loading input data..." << std::endl;
        char * const DO_file_name_out = (char*) "DO_HSS_Input_Data.data";
        Yw1 = new std::array<mpz_t, train_number>();
        Xi1Md = new std::array<std::array<mpz_t, attribute_number> ,train_number >();
        int read_result = server_read_input_DO_array(Xi1, Yw1, DO_file_name_out, KeySize);

        mpz_t test_tmp;
        mpz_init(test_tmp);
        mpz_set_ui(test_tmp, 150);

        std::cout << "convert begin..." << std::endl;
        Xi1M = new std::array<std::array<hss_men_p , attribute_number> , train_number >(); //input
        for(int i=0;i<train_number;i++)  //for each train sample
        {
            if(i % 100 == 0){
                gettimeofday(&t2,NULL);
                timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
                std::cout<<"Offline runtime on server for " << i << " train data is = "<< timeuse << " seconds. " << std::endl;  //seconds
            }
            //std::vector<hss_men_p> Xi1m; //men
            for(int t=0;t<attribute_number-1;t++){
                hss_men_p ynum1;
                men_share_init(&ynum1);  //men
                //convert_shares_paillier(&ynum1, mof11, &(*Xi1)[0][t], hek0, hpk->pk, server_id, KeySize, base_sk);
                //input_paillier_value(&(*Xi1)[i][t], hpk, test_tmp, KeySize);
                convert_shares_paillier(&(*Xi1M)[i][t], mof11, &(*Xi1)[i][t], hek0, hpk->pk, server_id, KeySize, base_sk);

                //if(t == attribute_number -1){
                    //check:
                    //std::cout << "Keysize: " << sizeOfDigits << std::endl;
                    //std::cout << "Xi1M: "  << (*Xi1M)[i][t].sharex << std::endl;
                    /*std::cout << "> " << (*Xi1)[i][t].encx << " <" << std::endl;

                    std::cout << "encdx:***" << std::endl;
                    for (int j = 0; j < (*KeySize); j++) {
                        std::cout << (*Xi1)[i][t].encdx[j] << ", ";
                    }
                    std::cout << "OK" << std::endl;
*/
                //}


                //convert_shares_paillier(&men_11, &setup->mof11, &input_1, &setup->hek0, setup->hpk.pk, 0, &setup->KeySize, setup->base_sk);
                //Xi1m.push_back(ynum1);
                //Xi2m.push_back(ynum2);
                //std::cout << std::endl;
            }
            //std::cout << std::endl;
            //if(i > 3) return;
            //Xi1M.push_back(Xi1m);
        }
        std::cout << "convert end" << std::endl;
        std::cout << "Check setup params:" << std::endl;
        std::cout << mof11->sharex << ", " << hek0->ds[3] << ", " << hpk->pk->n << ", " << server_id << ", " << *KeySize << ", " << base_sk << std::endl;

        mpz_t op1;
        mpz_init_set_ui(op1, 3);

        //Pre calculate d^2 for all Ix
        mpz_t tmp_X2;
        mpz_init(tmp_X2);

        //std::cout << "Check setup params:" << std::endl;
        //std::cout << mof11->sharex << ", " << hek0->ds[3] << ", " << hpk->pk->n << ", " << server_id << ", " << (*KeySize) << ", " << base_sk << std::endl;
        //std::cout << mof12->sharex << ", " << hek1->ds[3] << ", " << hpk->pk->n << ", " << 1 << ", " << (KeySize) << ", " << base_sk << std::endl;


        for(int i=0;i<train_number;i++)  //for each train sample
        {
            mpz_init(X2[i]);
            mpz_set_ui(X2[i], 0);
            if(i >= 5 and i < 10){
                std::cout << "-----------: " << std::endl;
            }
            for (int j = 0; j < attribute_number-1; j++) {  // calculate the distance between the t train sample and the prediction sample
                mult_shares_paillier_simpler(tmp_X2, &((*Xi1M)[i][j]), (*Xi1)[i][j].encx, hek0, hpk->pk, server_id,(KeySize), base_sk);///????
                //mult_shares_paillier_simpler(re2, &men_12, input_2.encx, &setup->hek1, setup->hpk.pk, 1, &setup->KeySize, setup->base_sk);


                //try to dec (*Xi1)[i][j].encx
                //tmp_X2 << ", " << Xi1M[i][j].sharex << ", " << ((*Xi1)[i][j]).encx
                /*
                if(i >= 5 and i < 10 and j == attribute_number-1){
                    std::cout << tmp_X2 << ", " << (*Xi1M)[i][j].sharex << ", " << ((*Xi1)[i][j]).encx << std::endl;
                    //std::cout << "sharedx: " << std::endl;


                    for(int p=0; p< *KeySize; p++){
                        std::cout << ((*Xi1M)[i][j]).sharedx[p] << ", ";
                    }
                    std::cout << std::endl;
                    //return;
                }*/


                mpz_add(X2[i], X2[i], tmp_X2);

                mpz_init((*Xi1Md)[i][j]);
                mult_shares_paillier_prepare(&((*Xi1M)[i][j]), ((*Xi1Md)[i][j]), hek0, hpk->pk, server_id,
                                             (KeySize), base_sk);///????
            }
            if(i >= 5 and i < 10) {
                std::cout << "Additive share of the square sum of train_data[" << i << "] : " << X2[i] << std::endl;
            }
        }
        mpz_clear(tmp_X2);

        for(int i=0; i<pthread_max_number; i++){
            for(int j=0; j<attribute_number; j++){
                mpz_init(this->Qs[i][j]);
            }
            mpz_init(this->Q2s[i]);
            mpz_init(this->randoms[i]);

            mpz_init(this->refreshed_label[i]);
            mpz_init(this->final_label[i]);
        }

        //clear Xi1
        for (int i = 0; i < train_number; i++)  //for each train sample
        {
            for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi
                mpz_clear((*Xi1)[i][t].encx);
                for(int j=0; j<(*KeySize); j++){
                    mpz_clear((*Xi1)[i][t].encdx[j]);
                }
            }
        }
        delete [] Xi1;

    }

    //void start_knn();
    //void start_knn_origin();
    int newKNN(int pint);
    int set_query_values(mpz_t input_Q[], mpz_t input_Q2, mpz_t random, int i);
    int query(int pint);
    int oblivious_transfer(int pint);
    int refresh(int pint);
    int cheat_select_top_k(int pint, int select_i[], char * k_message);
        ~KNN_single(){
        delete [] X2;
        delete [] distance_perm;
        delete [] qX1;
        delete [] Yw1;
        delete [] Xi1Md;
        delete [] Xi1M;
        delete [] qX1;
    }
    //void Time_result();
    //int hey(int myb);

    };


#endif //PPKNN_KNN_SINGLE_H

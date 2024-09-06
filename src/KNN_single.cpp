//
// Created by kore on 9/12/21.
//
#include <string>
#include "KNN_single.h"
#include <random>


using namespace std;

int KNN_single::newKNN(int pint)//预测函数
{
    return 0;
}

int KNN_single::set_query_values(mpz_t input_Q[attribute_number], mpz_t input_Q2, mpz_t random, int i){
    //cout << " here" << endl;
    for(int j=0; j<attribute_number; j++){
        //cout << i << endl;
        //cout << j << endl;
        mpz_set(this->Qs[i][j], input_Q[j]);
    }
    //cout << " here2" << endl;
    mpz_set(this->Q2s[i], input_Q2);
    //cout << " here3" << endl;
    mpz_set(this->randoms[i], random);
    return 0;
}

int KNN_single::query(int pint){
    mpz_t * distance_origin = new mpz_t[online_const_num];
    mpz_t * label_origin = new mpz_t[online_const_num];

    cout << "Clear Complet." << endl;
    struct timeval t1,t2;
    double timeuse;
    gettimeofday(&t1,NULL);
    int dis=0;
    int midint1, midint2;
    double middou1, middou2;


    mpz_t tmpmpz;
    mpz_init(tmpmpz);


    mpz_t numw;
    mpz_init(numw);

    mpz_t distance, tmp_distance, add_Q2;
    mpz_init(distance);
    mpz_init(tmp_distance);
    mpz_init(add_Q2);

    //Simulating
    //for(int s_i=0;s_i<20;s_i++)  //for each train sample
    //{
    /*
    mpz_t tmpvalue1, tmpvalue2, tmpvalue3;
    mpz_init(tmpvalue1);
    mpz_init(tmpvalue2);
    mpz_init(tmpvalue3);

    mpz_set_ui(tmpvalue1, 11111);
    mpz_set_ui(tmpvalue2, 12222);
    mpz_set_ui(tmpvalue3, 33331);*/

    cout << "First check the values of query: " << endl;
    for (int j = 0; j < attribute_number-1; j++) {
        cout << Qs[pint][j] << "; ";
    }
    cout << endl;
    cout << Q2s[pint] << endl;

    cout << "Calculating distances... " << endl;
    mult_shares_paillier_simpler(add_Q2, mof11, Q2s[pint], hek0, hpk->pk,
                                 server_id, KeySize, base_sk);
    for (int i = 0; i < online_train_number; i++)  //for each train sample
    {
        /*
        if(i % 100 == 0){
            gettimeofday(&t2,NULL);
            timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
            std::cout<<"Offline runtime on server for " << i << " train data is = "<< timeuse << " seconds. " << std::endl;  //seconds
            //gettimeofday(&t1,NULL);
        }
         */
        mpz_set_ui(numw, 0);


        //cout << "Distance: " << distance <<  endl;
/*
        if(i>= 5 and i < 10){
            cout << "X2[ " << i << "] is " << X2[i] << endl;
            cout << "Q2[ " << i << "] is " << add_Q2 << endl;
        }
*/
        mpz_add(distance, add_Q2, X2[i]);
        //cout << "point 2" << endl;
        for (int j = 0; j < attribute_number-1; j++) {  // calculate the distance between the t train sample and the prediction sample
            mult_shares_paillier_efficient(tmp_distance, (*Xi1Md)[i][j], Qs[pint][j], hek0, hpk->pk,
                                         server_id, KeySize, base_sk);

            //if(i <2 and j <2) {
            //    cout << tmp_distance << ", " << Xi1M[i][j].sharex << ", " << Qs[pint][j] << endl;
            //}

            mpz_mul_ui(tmp_distance, tmp_distance, 2);
            mpz_sub(distance, distance, tmp_distance);
/*
            if(i>= 5 and i < 10){
                cout << "Distances[ " << i << "][" << j << "] is " << tmp_distance << endl;
            }*/
        }
        //cout << endl;
        mpz_init_set((qX1[i].dis), distance);
        //mpz_init_set((qX1[i].labels), (*Yw1)[i]); //origin !!
        mpz_init_set_ui((qX1[i].labels), i);

    }
    //}
    cout << "Distances for " << online_train_number << " records is over." << endl;

    //cout << "Value for [5]~[9] is :"<< endl;
    for (int i = 0; i < train_number; i++)  //for each train sample
    {
        cout << qX1[i].dis << ",  ";
    }
    cout << endl;




    //PRP
    /*
    int perm_len;
    if(online_train_number == 100) perm_len = 128;
    else if(online_train_number == 1000) perm_len = 1024;
    else if(online_train_number == 10000) perm_len = 16384;

    mpz_t prp_key;
    mpz_init_set_ui(prp_key, 1001);
*/
    //since we take const = 8192, not enough for 10000, we short it...
    /*
    for(int i=0; i<online_const_num; i++){
        if(i<online_train_number) mpz_init_set(distance_origin[i], qX1[i].dis);
        else mpz_init_set_ui(distance_origin[i], 999);
        mpz_init(distance_perm[i]);
    }*/
    //prand_perm(distance_perm, distance_origin, prp_key, online_const_num);

    cout << "Psuedorandom Shuffling..." << endl;
    for(int i=0; i<online_const_num; i++){
        mpz_init_set(distance_origin[i], qX1[i].dis);
        mpz_init_set(label_origin[i], qX1[i].labels);
    }

    int train_len = online_const_num;

    gmp_randstate_t state_shuffle;
    gmp_randinit_mt(state_shuffle);
    mpz_t key_shuffle;
    mpz_init_set_ui(key_shuffle, 20220402);
    gmp_randseed(state_shuffle, key_shuffle);
    //PRShuffle(distance_origin, train_len, state_shuffle);
    //PRShuffle(label_origin, train_len, state_shuffle);
    PRShuffle(distance_origin, label_origin, train_len, state_shuffle);

    cout << "Shuffled Labels :"<< endl;
    for (int i = 0; i < train_number; i++)  //for each train sample
    {
        cout << i << ":" << label_origin[i] << ",  ";
        sind[i] = mpz_get_ui(label_origin[i]);
    }
    cout << endl;

/*
 * delta[(const_num)*(const_num - 1)/2],
            rand_c[(const_num)*(const_num-1)/2];//delta[train_number][train_number-1];//delta[(train_number*(train_number-1))/2];



    //cout << "Simulate k-out-of-n OT." << endl;
    //delta calculate
    //mpz_t rand_c;
    //mpz_init_set_ui(rand_c, 37);
    int index = 0;
    for(int i=0; i<(const_num); i++){
        for(int j=0; j<(const_num)-1-i; j++){
            mpz_init(delta[index]);
            mpz_sub(delta[index], distance_perm[i], distance_perm[j]);
            mpz_init_set_ui(rand_c[index], 37+i+j);
            mpz_mul(delta[index], delta[index], rand_c[index]);
            mpz_mod(delta[index], delta[index], hpk->pk->n);
            index++;
        }
        //d
    }




    //output delta...



    */
    cout << "Calculating delta." << endl;
    int index = 0;
    array <mpz_t, ( online_const_num * (online_const_num-1) / 2 ) > * delta = new array <mpz_t, ( online_const_num * (online_const_num-1) / 2 ) >();
    //array <unsigned long, ( online_const_num * (online_const_num-1) / 2 ) > * rand_c = new array<unsigned long, ( online_const_num * (online_const_num-1) / 2 ) >();

    cout << "Mem ok" << endl;

    mpz_t randseed, randnum;
    mpz_init(randseed);
    mpz_init(randnum);

    mpz_set_ui(randseed, 20220407);
    index = 0;
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed(state, randseed);

    //const char * output_file1_name = "Compare_tree.data";
    char output_file1_name[128];
    sprintf(output_file1_name,"Compare_tree_%d.data",server_id);
    FILE * fp_out = NULL;
    fp_out = fopen(output_file1_name, "w+");
    unsigned int result_flag, sum_out = 0;

    for(int i=0; i<online_const_num; i++){
        for(int j=online_const_num-1; j>i; j--){
            mpz_init((*delta)[index]);
            mpz_sub((*delta)[index], distance_origin[i], distance_origin[j]);
            mpz_urandomb(randnum, state, message_bit_len);
            mpz_mul((*delta)[index], (*delta)[index], randnum);


            result_flag = mpz_out_raw(fp_out, (*delta)[index]);
            if(result_flag == 0){
                cerr << "Unable to write!" << endl;
                return -1;
            }else{
                sum_out += result_flag;
            }



            index++;
        }
        //cout << i << endl;
    }
    /*
    char end_char[2];
    end_char[0] = '#';
    end_char[1] = '\0';
    write(fp_out, end_char, 1);*/
    int ch;
    ch = '#';
    fputc(ch, fp_out);

    cout << "Result: Write " << sum_out+1 << " bytes." << endl;
    fclose(fp_out);


    for(int i=0; i<online_const_num; i++){
        mpz_clear(distance_origin[i]);
        mpz_clear(distance_perm[i]);
    }
/*
    for(int i=0; i< online_train_number; i++){
        mpz_clear((*delta)[i]);
        mpz_clear((*rand_c)[i]);
    }
*/
/*
    mpz_clear(tmpvalue1);
    mpz_clear(tmpvalue2);
    mpz_clear(tmpvalue3);
*/

    delete [] delta;
    //delete [] rand_c;

    delete [] distance_origin;
    delete [] label_origin;

    cout << "Clear." << endl;

    return 0;
}

int KNN_single::oblivious_transfer(int pint) {
    //calculate comparison tree
    //vector version:
    /*
     *
     * int index = 0;
    for(int i=0; i<(const_num); i++){
        for(int j=0; j<(const_num)-1-i; j++){
            mpz_init(delta[index]);
            mpz_sub(delta[index], distance_perm[i], distance_perm[j]);
            mpz_init_set_ui(rand_c[index], 37+i+j);
            mpz_mul(delta[index], delta[index], rand_c[index]);
            mpz_mod(delta[index], delta[index], hpk->pk->n);
            index++;
        }
        //d
    }

    array< array <mpz_t, const_num> , const_num> delta_distance, rand_c;
    for(int i=0; i<(const_num); i++){
        for(int j=0; j<(const_num); j++){
            mpz_init(delta_distance[i][j]);
            mpz_sub(delta_distance[i][j], distance_perm[i], distance_perm[j]);
            mpz_init_set_ui(rand_c[i][j], 37+i+j);
            mpz_mul(delta_distance[i][j], delta_distance[i][j], rand_c[i][j]);
            mpz_mod(delta_distance[i][j], delta_distance[i][j], hpk->pk->n);
            //index++;
        }
        //d
    }*/


    cout << "Simulate k-out-of-n OT." << endl;

    mpz_t p, q, g, h;
    mpz_init(p);
    mpz_init(q);
    mpz_init(g);
    mpz_init(h);

    //should be initialized by setup..
    mpz_set_ui(q, 5);
    mpz_set_ui(p, 11);
    mpz_set_ui(g, 3);
    mpz_set_ui(h, 7);

    mpz_t tmp, tmp2, tmp3, rand;
    mpz_init(tmp);
    mpz_init(tmp2);
    mpz_init(tmp3);
    mpz_init(rand);

    int choice[k_value];
    //for(int i=0;i)

    //mpz_t messages[const_num];

    //init
    //for(int i=0; i< const_num; i++){
    //    mpz_init(delta[i]);
    //}

    cout << "Init As." << endl;
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    mpz_t As[k_value];
    for(int i=0; i< k_value; i++){
        mpz_init(As[i]);
        mpz_set_ui(As[i], i+k_value);
    }


    cout << "Init array" << endl;
    array< array<mpz_t, 2>, online_train_number> * ci = new array< array<mpz_t, 2>, online_train_number>();
    double expint;
    //init
    for(int i=0; i< online_train_number; i++){
        //cout << i << endl;
        mpz_init((*ci)[i][0]);
        mpz_init((*ci)[i][1]);

        mpz_urandomm(rand, state, q);
        mpz_powm((*ci)[i][0], g, rand, p);
        //cout << i << endl;
        mpz_set_ui(tmp2, 1);
        for(int j=0; j< k_value; j++){
            expint = pow(i+1, j+1);
            mpz_powm_ui(tmp3, As[j], expint, p);
            mpz_mul(tmp2, tmp2, tmp3);
        }
        //cout << i << endl;
        expint = pow(i, k_value);
        mpz_mul(tmp3, g, h);
        mpz_powm_ui(tmp3, tmp3, expint, p);
        mpz_mul(tmp2, tmp2, tmp3);
        //cout << i << endl;
        mpz_powm(tmp3, tmp2, rand, p);
        mpz_mul((*ci)[i][1], tmp2, tmp3);   //in fact should be permed labels
    }



    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(g);
    mpz_clear(h);

    mpz_clear(tmp);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
    mpz_clear(rand);

    for(int i=0; i< online_train_number; i++){
        mpz_clear((*ci)[i][0]);
        mpz_clear((*ci)[i][1]);
    }

    for(int i=0; i< k_value; i++){
        mpz_clear(As[i]);
    }

    delete [] ci;

    return 0;
}

int KNN_single::refresh(int pint){

    //mult_shares_paillier_simpler(tmp_X2, &(Xi1M[i][j]), final_label[pint], hek0, hpk->pk, server_id,(KeySize), base_sk);///????
    mult_shares_paillier_simpler(final_label[pint], mof11, refreshed_label[pint], hek0, hpk->pk, server_id, KeySize, base_sk);

    return 0;
}

int KNN_single::cheat_select_top_k(int pint, int select_i[], char * k_message){
        mpz_t top_k_result[k_value];
        char tmp[1024];
        cout << "Top k labels:" << endl;
        const char * sp = "-";

        //k_message[0] = '\0';
        cout << "[0]:" << sind[select_i[0]] << endl;
        mpz_init(top_k_result[0]);
        mult_shares_paillier_simpler(top_k_result[0], mof11, (*Yw1)[sind[select_i[0]]], hek0, hpk->pk, server_id, KeySize, base_sk);
        mpz_get_str(k_message, 10, top_k_result[0]);

        for (int i = 1; i < k_value; i++) {
            //cout << select_i[i];// << ": ";
            //cout << .sharex;
            cout << "[" << i <<  "]:" << sind[select_i[i]] << endl;
            mpz_init(top_k_result[i]);
            mult_shares_paillier_simpler(top_k_result[i], mof11, (*Yw1)[sind[select_i[i]]], hek0, hpk->pk, server_id, KeySize, base_sk);
            mpz_get_str(tmp, 10, top_k_result[i]);
            strcat(k_message, sp);
            strcat(k_message, tmp);
        }
        //cout << endl;
        printf("Result k messages: %s", k_message);


        return 0;
    }
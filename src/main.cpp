#include <iostream>
#include <vector>
#include <gmpxx.h>
#include <cmath>
#include<sys/time.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <thread>
#include<arpa/inet.h>

#include <queue>
#include <functional>
using namespace std;


#include "KNN_single.h"
#include "readData.h"
#include "Poly.h"
#include "polynomial.h"
#include "spec_sort.h"

#define CLIENT 3

#define MAXLINE 1024
#define SERVER_PORT 1234

#define MAXDATASIZE 1024    //缓冲区大小

using namespace std;

int NUM_IMAGES = 10240;
int PARTY;

//std::mutex mtx;

int pp_cmp_tree[10][10];


typedef struct SockInfo{
    int fd;
    struct sockaddr_in addr;
    socklen_t addr_len;
    pthread_t id;
}SockInfo;


struct knnInfo{
    int fd;
    struct sockaddr_in addr;
    socklen_t addr_len;
    pthread_t id;
    KNN_single * server = NULL;
};

int Simulate_case = 0; //1 if S1-S2 ; 2 if C-S1 ; 3 if C-S2;

int comm_times_SS = 0, comm_times_CS = 0; //count call times
int comm_overhead_SS = 0, comm_overhead_CS = 0;  //count mpz_t
int n=100, m=8, k_num=5, l=32, w=20; //l=bits of encryption
int enc_times = 0, dec_times = 0;

int server_id;
int server_port;

bool cmpt(node a,node b)
{
    return a.dis<b.dis;
}

template <typename T>
void div_m2c(vector<vector<T> > & Mat2,vector<vector<T> > & Mat1,\
double constNum, int IMAX, int JMAX){

    //Mat2.clear();
    vector<vector<T> > tempm;
    for(int i = 0; i < IMAX; ++i){
        vector<T> tempv;
        for(int k = 0; k < JMAX; ++k){

            tempv.push_back( (Mat1[i][k] / constNum));

        }
        //Mat2[i].clear();
        tempm.push_back(tempv);
    }
    Mat2.clear();
    for(int i = 0; i < IMAX; ++i){
        Mat2.push_back(tempm[i]);
    }
}

template <typename T>
void div_v2c(vector<T> &Vec2,vector<T> &Vec1,\
double constNum, int IMAX){

    vector<T> tempv;
//cout << "Who are you!"<< constNum << endl;
    for(int i = 0; i < IMAX; ++i){
        //if(Vec1[i] != 0) cout << Vec1[i]<< " " ;
        tempv.push_back((T) (Vec1[i] / constNum));
        //if(tempv[i] != 0) cout << tempv[i]<< " " ;
        //cout << Vec1[i]<<" ";
    }
    Vec2.clear();
    for(int i = 0; i < IMAX; ++i){
        Vec2.push_back(tempv[i]);
        //if(tempv[i] != 0) cout << tempv[i]<< " " ;
    }

}

template <typename T>
void mul_v2c(vector<T> &Vec2,vector<T> &Vec1,\
double constNum, int IMAX){

    vector<T> tempv;
//cout << "Who are you!"<< constNum << endl;
    for(int i = 0; i < IMAX; ++i){
        //if(Vec1[i] != 0) cout << Vec1[i]<< " " ;
        tempv.push_back((T) (Vec1[i] * constNum));
        //if(tempv[i] != 0) cout << tempv[i]<< " " ;
        //cout << Vec1[i]<<" ";
    }
    Vec2.clear();
    for(int i = 0; i < IMAX; ++i){
        Vec2.push_back(tempv[i]);
        //if(tempv[i] != 0) cout << tempv[i]<< " " ;
    }

}
template <typename T>
void v_expand(vector<T> &Vec1,\
double constNum, int IMAX){
    for(int i = 0; i < IMAX; ++i){
        Vec1[i] = (T)((int) (Vec1[i] * constNum));
    }

}

template <typename T>
void m_expand(vector<vector<T> > & Mat1,\
double constNum, int IMAX, int JMAX){
    for(int i = 0; i < IMAX; ++i){
        for(int k = 0; k < JMAX; ++k){
            Mat1[i][k] = (T)((int) (Mat1[i][k] * constNum));
        }
    }

}
/*
void exec_a_query(KNN_single &server, int pint, int fd){//, mpz_t (*input_Q)[attribute_number], mpz_t input_Q2){
    // , Q, Q2
    //mpz_t input_Q[attribute_number], input_Q2;
    server.query(pint);


    char message[1024];
    memset(message,0,sizeof(message));
    message[1024-1] = '\0';
    strcat(message, "Cal OK!");
    message[1024-1] = '\0';

    mtx.lock();
    if (-1 != send(fd, message ,1024,0)){
        cout<<"send msg: \""<< message <<" \" success"<<endl;
        close(fd);
    }
    mtx.unlock();


}*/



int server_read_input(mpz_t Q[attribute_number], mpz_t Q2, mpz_t random, char const * filepath){
    unsigned int result_flag, sum_in = 0;

    FILE *fp_in = NULL;
    fp_in = fopen(filepath, "r");


    mpz_t tmpNum;
    mpz_init(tmpNum);

    for (int t = 0; t < attribute_number-1; t++) { //read a train samples t to Xi



        result_flag = mpz_inp_raw(tmpNum, fp_in);
        if (result_flag == 0) {
            cerr << "Unable to read!" << endl;
            return -1;
        } else {
            sum_in += result_flag;
        }

        mpz_set(Q[t], tmpNum);
    }

    result_flag = mpz_inp_raw(tmpNum, fp_in);
    if (result_flag == 0) {
        cerr << "Unable to read!" << endl;
        return -1;
    } else {
        sum_in += result_flag;
    }

    mpz_set(Q2, tmpNum);


    result_flag = mpz_inp_raw(tmpNum, fp_in);
    if (result_flag == 0) {
        cerr << "Unable to read!" << endl;
        return -1;
    } else {
        sum_in += result_flag;
    }

    mpz_set(random, tmpNum);

    cout << "Result: Read " << sum_in << " bytes." << endl;
    fclose(fp_in);

    return 0;
}

int testIntoFunc(){
    int mydream;
    mpz_t ptime;

    return 0;
}

void* pthread_Callback(void* arg){
    knnInfo *info = (knnInfo*)arg;

    while(1)
    {
        char ip[64]={0};
        /*
        mpz_t Q[attribute_number];
        mpz_t Q2;
        mpz_t result_values[train_number+100];

        for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi
            mpz_init(Q[t]);
            mpz_set_ui(Q[t], t);
        }
        mpz_init(Q2);
        mpz_set_ui(Q2, 1000);
        for(int k=0; k<train_number; k++){
            mpz_init(result_values[k]);
            mpz_set_ui(result_values[k], 0);
        }
        */
        /*
        char buf[overhead];
        memset(buf, '*', overhead);
        buf[overhead] = '\0';

        if(PARTY == ALICE){
            if( send(fd, buf, strlen(buf), 0) < 0){
                printf("send msg error: %s(errno: %d)\n", strerror(errno), errno);
                return 0;
            }
        }else{
            int n = recv(fd, buf, MAXLINE, 0);
            buf[n] = '\0';
            printf("recv msg from client: %s\n", buf);
        }
        */

        char buf[300*(attribute_number+1)+1024] = {0};

        char * delta = new char[4*(train_number*(train_number-1)/2)+1024];
        memset(delta,'$',4*(train_number*(train_number-1)/2));
        delta[4*(train_number*(train_number-1)/2)-1] = '\0';

        char * ot1 = new char[16*train_number+1024];
        char * ot2 = new char[16*train_number+1024];
        char * refreshed_recv = new char[16+1];
        char * refreshed_return = new char[16+1];


        //testIntoFunc();

        printf("New Client IP: %s, Port: %d\n",
               inet_ntop(AF_INET,&info->addr.sin_addr.s_addr,ip,sizeof(ip)),
               ntohs(info->addr.sin_port));


/*
        int recv_num;
        char recvbuf[MAXDATASIZE], sendbuf[MAXDATASIZE];

        while (true) {
            recv_num = recv(info->fd, recvbuf, MAXDATASIZE,0);
            if(recv_num == -1)
            {
                perror("read err");
                pthread_exit(NULL);
            }else if(recv_num == 0)
            {
                printf("客户端断开连接\n");
                close(info->fd);
                break;
            }else{
                recvbuf[recv_num] = '\0';
                printf("Received client( %s ) message: %s",inet_ntoa(info->addr.sin_addr) , recvbuf);
            }
        }
        */



        int len = read(info->fd,buf,300*(attribute_number+1+1));
        if(len == -1)
        {
            perror("read err");
            pthread_exit(NULL);
        }else if(len == 0)
        {
            printf("客户端断开连接\n");
            close(info->fd);
            break;
        }else{
            //buf[len] = '\0';
            //printf("recv1: %s\n",buf);
            printf("[1] Receive %d Bytes (query) from client.\n",len);
            //=============distance calculaate==============

            mpz_t preQ[attribute_number];
            mpz_t preQ2;
            mpz_t random;

            for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi
                mpz_init(preQ[t]);
            }
            mpz_init(preQ2);
            mpz_init(random);

            const char sep[2] = "-";
            char * token;
            int lab_ind = 1;

            /* 获取第一个子字符串 */
            token = strtok(buf, sep);
            if(token != NULL){
                mpz_set_str(preQ[0], token, 36);
            }else{
                cout << "Sep Error!" << endl;
                cout << buf << endl;
            }

            /* 继续获取其他的子字符串 */
            while( token != NULL ) {
                if(lab_ind >= attribute_number-1){
                    //cout << "Too many selections!" << endl;
                    break;
                }
                cout  << "Extract " << lab_ind << " labels..." << endl;
                printf( "%s\n", token );
                token = strtok(NULL, sep);
                mpz_set_str(preQ[lab_ind], token, 36);
                lab_ind++;
            }

            token = strtok(NULL, sep);
            mpz_set_str(preQ2, token, 36);

            token = strtok(NULL, sep);
            mpz_set_str(random, token, 36);

            //*******

            //testIntoFunc();
            //char const * query_file_path = "client_query.data";
            //server_read_input(preQ, preQ2, random, query_file_path);
            cout << "Receive query from clients, q[0] is " << preQ[0] << endl;
            info->server->set_query_values(preQ, preQ2, random, info->id);
            cout << info->id << endl;
            info->server->query(info->id);

            //info->server->hey(info->id);
            //testIntoFunc();

            //S randomize the distances
            /*
            for(int k=0; k<train_number; k++){
                mpz_add_ui(result_values[k], result_values[k], 60);
            }
            */
            //===============================
            //memset(delta, '#', 16*(train_number*(train_number-1)/2));
            //delta[16*(train_number*(train_number-1)/2 - 1] = '\0';
            //write(info->fd,delta,4*(train_number*(train_number-1)/2));  //AVAILABLE

            char  buffer[MAXLINE];
            FILE *fq;
            char output_file1_name[128];
            sprintf(output_file1_name,"/home/kore/CLionProjects/KNN_Server/cmake-build-debug/Compare_tree_%d.data",server_id);

            if( ( fq = fopen(output_file1_name,"rb") ) == NULL ){
                printf("File open error.\n");
                //return -1;
                //close(sockfd);
                exit(1);
            }

            bzero(buffer,sizeof(buffer));
            while(!feof(fq)){
                len = fread(buffer, 1, sizeof(buffer), fq);
                if(len != write(info->fd, buffer, len)){
                    printf("write error.\n");
                    break;
                }
            }
            //sleep(2000);
            //close(sockfd);
            fclose(fq);

            for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi
                mpz_clear(preQ[t]);
            }
            mpz_clear(preQ2);
            delete [] delta;
        }

        int select_i[k_value];
        len = read(info->fd,ot1,(16+1)*(k_value));
        if(len == -1)
        {
            perror("read err");
            pthread_exit(NULL);
        }else if(len == 0)
        {
            printf("客户端断开连接\n");
            close(info->fd);
            break;
        }else{
            //ot1[len] = '\0';
            //printf("recv2:%d, %s\n", len, ot1);
            printf("[2] Receive %d Bytes (oblivious transfer) from client.\n",len);

            const char sep[2] = "-";
            char *token;
            mpz_t select[k_value];

            int lab_ind = 1;

            /* 获取第一个子字符串 */
            token = strtok(ot1, sep);
            if(token != NULL){
                mpz_init(select[0]);
                mpz_set_str(select[0], token, 10);
                select_i[0] = atoi(token);
            }else{
                cout << "Sep Error!" << endl;
                cout << ot1 << endl;
            }


            /* 继续获取其他的子字符串 */
            while( token != NULL ) {
                if(lab_ind >= k_value){
                    cout << "Too many selections!" << endl;
                    break;
                }
                cout  << "Extract " << lab_ind << " labels..." << endl;
                printf( "%s\n", token );
                token = strtok(NULL, sep);
                mpz_init(select[lab_ind]);
                mpz_set_str(select[lab_ind], token, 10);
                select_i[lab_ind] = atoi(token);
                lab_ind++;
            }
            cout << "Show selected indexes: " << endl;
            for(int j=0; j<k_value; j++){
                cout << j << ": " << select[j] << " ";
            }
            cout  << endl;
            char k_message[16*k_value];
            info->server->cheat_select_top_k(info->id, select_i, k_message);


            //=========proceed the returned labels=============
            /*
            for(int i=0;i<train_number;i++){
                mpz_add_ui(result_values[i], result_values[i], 11);
            }
             */
            //==================================================
            info->server->oblivious_transfer(info->id);

            //memset(ot2, '*', 16*train_number);
            //ot2[16*train_number-1] = '\0';
            strncpy(ot2, k_message, 16*train_number);
            //sprintf(ot,"%1600.10s \n","Good bye!");
            write(info->fd,ot2,16*train_number);

        }

        len = read(info->fd,refreshed_return,16);
        if(len == -1)
        {
            perror("read err");
            pthread_exit(NULL);
        }else if(len == 0)
        {
            printf("客户端断开连接\n");
            close(info->fd);
            break;
        }else {
            refreshed_recv[len] = '\0';
            //printf("recv2:%d, %s\n", len, ot1);
            printf("[3] Receive %d Bytes (refresh label) from client.\n", len);
            mpz_set_ui(info->server->refreshed_label[info->id], info->id);
            info->server->refresh(info->id);




            memset(refreshed_return, '/', 16);
            refreshed_return[16] = '\0';
            //sprintf(ot,"%1600.10s \n","Good bye!");
            write(info->fd,refreshed_return,16);

            close(info->fd);
            printf("task finished.\n");
            break;
        }

        delete [] ot1;
        delete [] ot2;
    }

    pthread_exit(NULL);
}


int client_enc_input(mpz_t Q[], mpz_t Q2, int q[], SetupPhase * setup){
    int result, sum=0;
    mpz_t tmp;
    mpz_init(tmp);
    for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi
        mpz_set_ui(tmp, q[t]);
        result = paillier_encrypt(Q[t], tmp, setup->hpk.pk);
        sum += q[t] * q[t];
        //mpz_set_ui(Q[t], t);
    }
    mpz_set_ui(tmp, sum);
    result = paillier_encrypt(Q2, tmp, setup->hpk.pk);
    //if(server_id == 0) mpz_set_ui(Q2, 33333);
    //else mpz_add_ui(Q2, tmp, 33333);

    mpz_clear(tmp);

    return 0;
}

int client_write_input(mpz_t Q[], mpz_t Q2, char * filepath){
    unsigned int result_flag, sum_out = 0;

    FILE * fp_out = NULL;
    fp_out = fopen(filepath, "w+");

    for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi

        result_flag = mpz_out_raw(fp_out, Q[t]);
        if(result_flag == 0){
            cerr << "Unable to write!" << endl;
            return -1;
        }else{
            sum_out += result_flag;
        }
    }
    result_flag = mpz_out_raw(fp_out, Q2);
    if(result_flag == 0){
        cerr << "Unable to write!" << endl;
        return -1;
    }else{
        sum_out += result_flag;
    }

    cout << "Result: Write " << sum_out << " bytes." << endl;
    fclose(fp_out);

    return 0;
}


/*
int client_write_input(mpz_t Q[], mpz_t Q2, int q[], SetupPhase * setup, int server_id){
    int result, sum=0;
    mpz_t tmp;
    mpz_init(tmp);
    for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi
        mpz_set_ui(tmp, q[t]);
        result = paillier_encrypt(Q[t], tmp, setup->hpk.pk);
        sum += q[t] * q[t];
        //mpz_set_ui(Q[t], t);
    }
    mpz_set_ui(tmp, sum);
    //result = paillier_encrypt(Q2, tmp, setup->hpk.pk);
    if(server_id == 0) mpz_set_ui(Q2, 33333);
    else mpz_add_ui(Q2, tmp, 33333);

    mpz_clear(tmp);

    return 0;
}*/



int server_evaluate(){

    cout << "========================================" << endl;
    cout << "Server " << server_id << " proceeding..." << endl;
    cout << "========================================" << endl;
    char const * file_path_his = (char *) "history_record.log";
    char const * KeyFilePath = (char*) "KeyFilePath.data";

    char * const file_name_server_out_0 = (char*) "Server_Result_Data_0.data";
    char * const file_name_server_out_1 = (char*) "Server_Result_Data_1.data";

    char * const C_file_name_out_0 = (char*) "C_HSS_Input_Data_0.data";
    char * const C_file_name_out_1 = (char*) "C_HSS_Input_Data_1.data";

    struct timeval t1,t2;
    double timeuse;

    gettimeofday(&t1,NULL);

    TrainingParams params;
    SetupPhase * setup = new SetupPhase(KeyFilePath, server_id);

    cout << "========================================" << endl;
    cout << "KNN Predicting" << endl;
    cout << "========================================" << endl;

    //struct node_paillier * qX1 = new node_paillier[train_number+100];//, qX2[train_number+100];

    KNN_single * kp = new KNN_single(&setup->mof11,  &setup->hpk, &setup->hek0, params, &setup->KeySize, setup->base_sk, server_id);

    //delete [] qX1;

    gettimeofday(&t2,NULL);
    timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
    cout<<"Offline runtime of server is = "<< timeuse << " seconds. " << endl;  //seconds

    cout << "========================================" << endl;
    cout << "Load ok, get online..." << endl;


    int port = server_port;
    // 创建套接字
    int listenfd = socket(AF_INET,SOCK_STREAM,0);
    struct sockaddr_in serv_addr;
    memset(&serv_addr,0,sizeof(serv_addr));
    serv_addr.sin_port = htons(port+server_id);
    serv_addr.sin_family = AF_INET;
    socklen_t serv_len = sizeof(serv_addr);
    int ret;
    cout << "Try to bind" << port+server_id << endl;
    if((ret =bind(listenfd,(struct sockaddr*)&serv_addr,serv_len)) != 0){
        perror("err: bind");
        exit(1);
    }
    if((ret =listen(listenfd,36))!=0)
    {
        perror("err: listen");
        exit(2);
    }



    printf("Start accept\n");


    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);

    knnInfo infos[pthread_max_number];  //////

    for (int j=0;j<sizeof(infos)/sizeof(knnInfo);++j)
        infos[j].fd = -1;


    int i;
    while(1)
    {
        for (i=0;i<sizeof(infos)/sizeof(knnInfo);++i)
        {
            if(infos[i].fd == -1)
                break;
        }
        if(i == pthread_max_number)
            break;

        struct sockaddr_in cli_addr;
        socklen_t cli_len = sizeof(cli_addr);

        // 主线程等待接收连接请求
        infos[i].fd =accept(listenfd,(struct sockaddr*)&cli_addr,&cli_len);

        infos[i].addr = cli_addr;
        infos[i].addr_len = cli_len;

        infos[i].server = kp;

        infos[i].id = i;


        // 创建子线程通讯
        pthread_t tid;
        pthread_create(&tid,&attr,pthread_Callback,&(infos[i]));

    }

    close(listenfd);
    pthread_attr_destroy(&attr);
    pthread_exit(NULL);


    return 0;

}

int HSS_MULT_TEST(){



    //setup
    SetupPhase * setup = new SetupPhase("KeyFilePath.data", 0);

    //load system setup? params


    //input two number
    mpz_t op1, op2, re1, re2;
    mpz_init_set_ui(op1, 100);
    mpz_init_set_ui(op2, 250);
    mpz_init_set_ui(re1, 0);
    mpz_init_set_ui(re2, 0);

    //encrypt
    hss_input_p input_1;
    input_paillier(&input_1, &setup->hpk, op1, &setup->KeySize);

    hss_input_p input_2;
    input_paillier(&input_2, &setup->hpk, op2, &setup->KeySize);

    //input_paillier_value(&input_2, &setup->hpk, op2, &setup->KeySize);

    //convert
    hss_men_p men_11;
    men_share_init(&men_11);  //men
    convert_shares_paillier(&men_11, &setup->mof11, &input_1, &setup->hek0, setup->hpk.pk, 0, &setup->KeySize, setup->base_sk);

    hss_men_p men_12;
    men_share_init(&men_12);  //men
    convert_shares_paillier(&men_12, &setup->mof12, &input_1, &setup->hek1, setup->hpk.pk, 1, &setup->KeySize, setup->base_sk);

    hss_men_p men_21;
    men_share_init(&men_21);  //men
    convert_shares_paillier(&men_21, &setup->mof11, &input_2, &setup->hek0, setup->hpk.pk, 0, &setup->KeySize, setup->base_sk);

    hss_men_p men_22;
    men_share_init(&men_22);  //men
    convert_shares_paillier(&men_22, &setup->mof12, &input_2, &setup->hek1, setup->hpk.pk, 1, &setup->KeySize, setup->base_sk);

    cout << "Check men of 10" << endl;
    cout << men_11.sharex << ", " << men_12.sharex << endl;

    std::cout << "Check setup params:" << std::endl;
    std::cout << setup->mof11.sharex << ", " << setup->hek0.ds[3] << ", " << setup->hpk.pk->n << ", " << 0 << ", " << (setup->KeySize) << ", " << setup->base_sk << std::endl;
    std::cout << setup->mof12.sharex << ", " << setup->hek1.ds[3] << ", " << setup->hpk.pk->n << ", " << 1 << ", " << (setup->KeySize) << ", " << setup->base_sk << std::endl;


    for(int k=0; k<1; k++){
        mpz_set_ui(op1, k+100);
        mpz_set_ui(op2, k+250);

        input_paillier_value(&input_1, &setup->hpk, op1, &setup->KeySize);
        input_paillier_value(&input_2, &setup->hpk, op2, &setup->KeySize);


        convert_shares_paillier(&men_11, &setup->mof11, &input_1, &setup->hek0, setup->hpk.pk, 0, &setup->KeySize, setup->base_sk);
        convert_shares_paillier(&men_12, &setup->mof12, &input_1, &setup->hek1, setup->hpk.pk, 1, &setup->KeySize, setup->base_sk);

        mult_shares_paillier_simpler(re1, &men_11, input_2.encx, &setup->hek0, setup->hpk.pk, 0, &setup->KeySize, setup->base_sk);
        mult_shares_paillier_simpler(re2, &men_12, input_2.encx, &setup->hek1, setup->hpk.pk, 1, &setup->KeySize, setup->base_sk);

        for(int p=0; p< (setup->KeySize); p++){
            std::cout << (men_11).sharedx[p] << ", " << (men_12).sharedx[p] << ";" << endl;
        }
        std::cout << "simple mult result:" << std::endl;
        std::cout << re1 << std::endl;
        std::cout << re2 << std::endl;
    }


    return 0;
}


int HSS_test(){
    //setup
    SetupPhase * setup = new SetupPhase();

    //load system setup? params


    //input two number
    mpz_t op1, op2, re1, re2;
    mpz_init_set_ui(op1, 10);
    mpz_init_set_ui(op2, 25);
    mpz_init_set_ui(re1, 0);
    mpz_init_set_ui(re2, 0);

    //encrypt
    hss_input_p input_1;
    input_paillier(&input_1, &setup->hpk, op1, &setup->KeySize);

    hss_input_p input_2;
    input_paillier(&input_2, &setup->hpk, op2, &setup->KeySize);

    //convert
    hss_men_p men_11;
    men_share_init(&men_11);  //men
    convert_shares_paillier(&men_11, &setup->mof11, &input_1, &setup->hek0, setup->hpk.pk, 0, &setup->KeySize, setup->base_sk);

    hss_men_p men_12;
    men_share_init(&men_12);  //men
    convert_shares_paillier(&men_12, &setup->mof12, &input_1, &setup->hek1, setup->hpk.pk, 1, &setup->KeySize, setup->base_sk);

    hss_men_p men_21;
    men_share_init(&men_21);  //men
    convert_shares_paillier(&men_21, &setup->mof11, &input_2, &setup->hek0, setup->hpk.pk, 0, &setup->KeySize, setup->base_sk);

    hss_men_p men_22;
    men_share_init(&men_22);  //men
    convert_shares_paillier(&men_22, &setup->mof12, &input_2, &setup->hek1, setup->hpk.pk, 1, &setup->KeySize, setup->base_sk);

    cout << "Check men of 10" << endl;
    cout << men_11.sharex << ", " << men_12.sharex << endl;

    std::cout << "Check setup params:" << std::endl;
    std::cout << setup->mof11.sharex << ", " << setup->hek0.ds[3] << ", " << setup->hpk.pk->n << ", " << 0 << ", " << (setup->KeySize) << ", " << setup->base_sk << std::endl;
    std::cout << setup->mof12.sharex << ", " << setup->hek1.ds[3] << ", " << setup->hpk.pk->n << ", " << 1 << ", " << (setup->KeySize) << ", " << setup->base_sk << std::endl;

    //add

    //mul
    //mult_shares_paillier_simpler
    mult_shares_paillier_simpler(re1, &men_11, input_1.encx, &setup->hek0, setup->hpk.pk, 0, &setup->KeySize, setup->base_sk);
    mult_shares_paillier_simpler(re2, &men_12, input_1.encx, &setup->hek1, setup->hpk.pk, 1, &setup->KeySize, setup->base_sk);
    std::cout << "simple mult result:" << std::endl;
    std::cout << re1 << std::endl;
    std::cout << re2 << std::endl;




    return 0;

    //PRP
    int train_1024 = 1024;
    cout << "Before perm: " << endl;
    mpz_t ds[train_1024], d[train_1024];
    for(int i=0; i< train_1024; i++){
        mpz_init_set_ui(d[i], i);
        mpz_init_set_ui(ds[i], 0);
        cout << d[i] << ", ";
    }
    cout << endl;
    prand_perm(ds, d, op1, train_1024);
    cout << "After perm: " << endl;
    for(int i=0; i< train_1024; i++){
        cout << ds[i] << ", ";
    }
    cout << endl;
    //output


    //
    int index;
    mpz_t cmp_tree[train_number][train_number], delta[(train_number*(train_number-1))], delta1[(train_number*(train_number-1))], delta2[(train_number*(train_number-1))];
    for(int i=0; i<train_number; i++){
        for(int j=0; j<train_number-1-i; j++){
            mpz_init(cmp_tree[i][i+j+1]);
            mpz_init_set_ui(delta[index], index+1);
            mpz_init_set_ui(delta[index], index+2);

            mpz_sub(cmp_tree[i][i+j+1], delta2[index], delta1[index]);
        }
    }

    int mins_index[k_value];
    int is_min[train_number] = {0};

    index = 0;

    for(int s=0; s<k_value; s++){
        mpz_t min;
        mpz_init(min);
        int mini = 0, tmp = 0;
        while(is_min[mini] == -1 && mini < train_number) mini++;

        //mpz_set(min, result_values[mini]);

        for(int i=mini;i<train_number;i++){
            for (int j = i+1; j < train_number; j++){
                if(mpz_cmp_ui(cmp_tree[i][j] , 0) > 0){
                    mini = j;
                    break;
                }
            }
            if( i == mini) break;
            else i = mini;
        }

        is_min[mini] = -1;
        mins_index[s] = mini;
    }

    //ot retrieve corresponding k labels
    mpz_t labels[k_value];
    for(int s=0; s<k_value; s++){
        mpz_init_set_ui(labels[s], s);
    }

    int max_w = 5, most_ind = 0;
    int pop_ind[5] = {0};
    mpz_t pop_sam[max_w];
    for(int p=0; p<max_w; p++){
        mpz_init_set_ui(pop_sam[p], p);
    }
    //mpz_init();
    for(int s=0; s<k_value; s++){
        for(int p=0; p<max_w; p++){
            if(mpz_cmp(labels[s], pop_sam[p]) == 0){
                pop_ind[p]++;
                if(pop_ind[p] > pop_ind[most_ind]){
                    most_ind = p;
                }
                break;
            }
        }
    }

    //HC send pop_sam[most_ind] to User

    //User
    //mpz_sub(label, labels, rc);
    //cout << "predict label is " << label << endl;



    return 0;
}

int k_out_of_n_OT_int_test(){
    Poly  x;
    Poly  y;

    char *p_path =  (char*) "poly1.data";
    char *q_path =  (char*) "poly2.data";

    x.readPoly(p_path);
    printf("p=");
    x.printPoly();

    y.readPoly(q_path);
    printf("q=");
    y.printPoly();

    char op = '+';
    //char op = '-';
    //char op = '*';

    Poly *res = x.algorithm(y, op);
    printf("The result is=");
    res->printPoly();
    delete res;


    return 0;
}

int poly_evaluate(mpz_t result, int type, int x){
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_ui(tmp, 0);

    mpz_set_ui(result, 0);

    int val;
    for(int i=0; i<k_value; i++){
        val = (int) pow(x, i);
        mpz_mul_ui(tmp, tmp, val); //simul.
        //mpz_mul_ui(tmp, a[i], val);
        mpz_add(result, result, tmp);
    }

    mpz_clear(tmp);
    return 0;
}

int find_gen(mpz_t p, mpz_t gen1, mpz_t gen2){
    mpz_t tmp, tmp2, tmp3, tmp4;
    mpz_init(tmp);
    mpz_init(tmp2);
    mpz_init(tmp3);
    mpz_init(tmp4);

    mpz_set_ui(tmp2, 2);
    mpz_sub_ui(tmp4, p, 1);
    while(mpz_cmp(tmp2, p) <= 0){
        mpz_add_ui(tmp2, tmp2, 1);
        mpz_powm(tmp, tmp2, tmp3, p);
        if(mpz_cmp(tmp, tmp4) == 0){
            mpz_set(gen1, tmp2);
            break;
        }
    }
    while(mpz_cmp(tmp2, p) <= 0){
        mpz_add_ui(tmp2, tmp2, 1);
        mpz_powm(tmp, tmp2, tmp3, p);
        if(mpz_cmp(tmp, tmp4) == 0){
            mpz_set(gen2, tmp2);
            break;
        }
    }

    mpz_clear(tmp);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
    mpz_clear(tmp4);
    return 0;
}

int k_out_of_n_OT_test(){
    //System Params: g, h, G
    mpz_t p, q, n, g, h;
    mpz_init(p);
    mpz_init(q);
    mpz_init(g);
    mpz_init(h);
    mpz_init(n);
/*
    gmp_randstate_t gmp_state;
    int msize=128;//长度

    gmp_randinit_mt(gmp_state);
    mpz_urandomb(q,gmp_state,msize);//生成128bit的随机数
    mpz_nextprime(q,q);//生成素数q
    //p=2*q+1,构造一个素数p，使得p-1的因子只有1,2,q，p-1四个
    mpz_mul_ui(p,q,2);
    mpz_add_ui(p,p,1);
    //求生成元，1到p-1之间的随机数a不是单位元、阶不为2和q，
    //则其阶为p-1，进而证明得a为p的生成元
    mpz_t a_2,a_q;
    mpz_init(a_2);
    mpz_init(a_q);
    do{
        mpz_urandomm(g,gmp_state,p);//生成1到p-1的随机数
        mpz_powm_ui(a_2,g,2,p);//a_2=a^2 % p
        mpz_powm(a_q,g,q,p);//a_q=a^q % p
        if((mpz_cmp_ui(g,1)!=0)&&(mpz_cmp_ui(a_2,1)!=0)&&(mpz_cmp_ui(a_q,1)!=0))
            break;
    }while(1);

    do{
        mpz_urandomm(h,gmp_state,p);//生成1到p-1的随机数
        mpz_powm_ui(a_2,h,2,p);//a_2=a^2 % p
        mpz_powm(a_q,h,q,p);//a_q=a^q % p
        if((mpz_cmp_ui(h,1)!=0)&&(mpz_cmp_ui(a_2,1)!=0)&&(mpz_cmp_ui(a_q,1)!=0))
            break;
    }while(1);

    //mpz_urandomm(Xa,gmp_state,p);//生成私钥Xa
    //mpz_powm(Ya,a,Xa,p);//生成公钥Ya,Ya=a^Xa % p
*/

    //gen_prime(q, 512/2);
    mpz_set_ui(q, 13);
    mpz_mul_ui(p, q, 2);
    mpz_add_ui(p, p, 1);
    mpz_mul(n, p, q);
    //set g = 1+n
    //mpz_add_ui(g, q, 1);
    //mpz_mod(g, g, q);
    //mpz_sub_ui(h, q, 1);
    mpz_set_ui(g, 7);
    mpz_set_ui(h, 11);

    //???
    //mpz_set(p, q);

    cout << "q: " << q << endl;
    cout << "p: " << p << endl;
    //cout << "n: " << n << endl;
    cout << "g: " << g << endl;
    cout << "h: " << h << endl;

    //Sender's messages: m1, ..., mn
    mpz_t messages[train_number];
    for(int i=0; i< train_number; i++){
        mpz_init(messages[i]);
        mpz_set_ui(messages[i], i+3);
        mpz_mod(messages[i], messages[i], q);
    }

    //Receiver's choices: sigma1, ..., sigmak
    int choice[k_value];
    mpz_t recv_messages[k_value];
    for(int i=0; i< k_value; i++) {
        //mpz_init(choice[i]);
        choice[i] = i+1;
        mpz_init(recv_messages[i]);
    }

    //Init//
    mpz_t tmp, tmp2, tmp3, rand;
    mpz_init(tmp);
    mpz_init(tmp2);
    mpz_init(tmp3);
    mpz_init(rand);
    mpz_t co_a[k_value], co_b[k_value];
    for(int i=0; i< k_value; i++){
        mpz_init(co_a[i]);
        mpz_init(co_b[i]);
    }


    //Step1: R chooses two poly
    //choose random a_i in Zq
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    for(int i=0; i< k_value; i++){
        mpz_urandomm(co_a[i], state, q);
    }

    vector<int> ra;
    vector<mpz_class> coa;
    for(int i=0; i<k_value; i++){
        ra.push_back(i);
        coa.push_back(mpz_class(co_a[i]));
    }
    ra.push_back(k_value);
    coa.push_back(mpz_class(1));

    polynomial * a_poly = new polynomial(k_value+1, ra, coa);;
    mpz_class result_params_a[k_value+1] = {0};
    cout << *a_poly << endl;
    polynomial a_result = a_poly->sort(result_params_a);
    cout << "f(x): ";
    cout << a_result << endl;
    for(int i=0; i<k_value; i++){
        cout << "a" << i << "=" << co_a[i] << endl;
    }
    //input choice two construct b_i:
    vector<int> r;
    vector<mpz_class> co;

    r.push_back(0);
    co.push_back(0);
    r.push_back(1);
    co.push_back(1);

    vector<polynomial*> choice_poly;

    for(int i=0; i<k_value; i++){
        co[0] = -choice[i];
        //co[0] = (co[0] + 13) % 13; //?


        polynomial * polyi = new polynomial(2, r, co);
        choice_poly.push_back(polyi);
    }

    polynomial product_poly = (*choice_poly[0]);
    for(int i=1; i<k_value; i++){
        product_poly = (product_poly) * (*choice_poly[i]);
    }

    mpz_class result_params[k_value+1] = {0};
    cout << product_poly << endl;
    polynomial result = product_poly.sort(result_params);
    cout << "f'(x): ";
    cout << result << endl;

    for(int i=0; i<k_value; i++){
        mpz_set(co_b[i], result_params[i].get_mpz_t());
        mpz_mod(co_b[i], co_b[i], q);
        cout << "b" << i << "=" << co_b[i] << endl;
    }
    /*
    mpz_set(co_b[k_value], result_params[k_value].get_mpz_t());
    cout << "b" << k_value << "=" << co_b[k_value] << endl;
     */

    //polynomial(int K, std::vector<int> r, std::vector<double> co)



    //Step2: R cal and send A0, ..., Ak-1
    mpz_t As[k_value];
    //init

    cout << "Masking..." << endl;
    for(int i=0; i< k_value; i++){
        mpz_init(As[i]);
        mpz_powm(tmp, g, co_a[i], p);
        //cout << h << ", " <<  co_b[i] << ", " <<  p << endl;
        mpz_powm(tmp2, h, co_b[i], p); ///
        //cout << tmp2 << ", " <<  h << ", " <<  co_b[i] << endl;
        mpz_mul(As[i], tmp, tmp2);
        //cout << As[i] << ", " <<  tmp << ", " <<  tmp2 << endl;
    }
    cout << endl;


    //Step3: S cal c1, ... ,cn
    mpz_t ci[train_number][2];
    mpz_t Bs[train_number][2];
    long expint;
    //init
    cout << "Check if is correct B" << endl;
    for(int i=0; i< train_number; i++){
        mpz_init(ci[i][0]);
        mpz_init(ci[i][1]);

        mpz_urandomm(rand, state, q);
        mpz_powm(ci[i][0], g, rand, p);
        //cout << ci[i][0] << ", " << g << ", " << rand << endl;

        mpz_set_ui(tmp2, 1);
        for(int j=0; j< k_value; j++){
            expint = pow(i+1, j);

            mpz_powm_ui(tmp3, As[i], expint, p);
            mpz_mul(tmp2, tmp2, tmp3);
        }
        //cout << tmp2 << endl;
        expint = pow(i+1, k_value);
        mpz_mul(tmp3, g, h);
        //cout << tmp3 << endl;
        mpz_powm_ui(tmp3, tmp3, expint, p);
        mpz_mul(tmp2, tmp2, tmp3);
        mpz_mod(tmp2, tmp2, p);
        //cout << tmp3 << endl;

        mpz_init(Bs[i][0]);
        //mpz_init(Bs[i][1]);
        mpz_set(Bs[i][0], tmp2);
        //mpz_set(Bs[i][1], tmp2);
        mpz_powm(tmp3, tmp2, rand, p);
        mpz_mul(ci[i][1], messages[i], tmp3);
        //cout << ci[i][1] << ", " << messages[i] << ", " << tmp3 << ", " << tmp2 << ", " << rand << endl;
    }


    cout << "Poly result comparison:" << endl;
    for(int i=0; i< 5; i++){
        mpz_set_ui(tmp, i+1);
        a_result.evaluate(tmp2, tmp, q);
        //mpz_mod(tmp2, tmp2, q);
        mpz_powm(tmp3, g, tmp2, p);
        cout << "[" << i << "]: " << tmp3 << " -- ";
        result.evaluate(tmp2, tmp, q);
        //mpz_mod(tmp2, tmp2, q);
        //cout << tmp2 << ", ";
        mpz_powm(tmp, h, tmp2, p);
        //cout << tmp << endl;
        mpz_mul(tmp, tmp, tmp3);
        mpz_mod(tmp, tmp, p);
        cout  << Bs[i][0] << endl;
    }



/*
    for(int i=0; i< 20; i++){
        cout << ci[i][0] << ", " << ci[i][1] << endl;
    }
*/


    //Step4: S send c1, ... ,cn



    //Step4: R cal m_{delta_1}, ..., m_{delta_k}
    for(int i=0; i< k_value; i++){
        mpz_set_ui(tmp2, choice[i]);
        a_result.evaluate(tmp3, tmp2, q);
        //result.evaluate(tmp, tmp2, q);
        //cout << "are you 0 :" << tmp << endl;
        //poly_evaluate(tmp3, 0, choice[i]);

        mpz_powm(tmp, ci[choice[i]][0], tmp3, p);
        mpz_invert(tmp, tmp, p);
        mpz_mul(recv_messages[i], ci[choice[i]][1], tmp);
        mpz_mod(recv_messages[i], recv_messages[i], q);

        cout << recv_messages[i] << ", " << choice[i] << ", " << ci[choice[i]][1] << ", " << tmp << ", " << endl;
        cout << "Result messages[" << choice[i] << "] = " << recv_messages[i] << " while origin messages is " << messages[choice[i]] << endl;
    }
    //return 0;


    //for() mul;->get b[]







    //poly_evaluate(tmp3, 0, i);//f(i)
    //        mpz_pow(tmp2, g, tmp3);
    //        poly_evaluate(tmp3, 1, i);//f'(i)
    //        mpz_pow(tmp2, h, fs[i]);
    //        mpz_mul(As[i], tmp, tmp2);



    //get random ci[][0] = g^ki

    //calculate ci[][1] = mul(pow()).....

    //evaluate(tmp, choice[i])
    //mpz_powm(ci[choice[i]][0], tmp)
    //mpz_inv()


    //mpz_mul()
    mpz_clear(tmp);
    mpz_clear(tmp2);
    mpz_clear(tmp3);

    return 0;
}

int k_out_of_n_OT_test2(){
    //System Params: g, h, G
    mpz_t p, q, n, g, h;
    mpz_init(p);
    mpz_init(q);
    mpz_init(g);
    mpz_init(h);
    mpz_init(n);
/*
    gmp_randstate_t gmp_state;
    int msize=128;//长度

    gmp_randinit_mt(gmp_state);
    mpz_urandomb(q,gmp_state,msize);//生成128bit的随机数
    mpz_nextprime(q,q);//生成素数q
    //p=2*q+1,构造一个素数p，使得p-1的因子只有1,2,q，p-1四个
    mpz_mul_ui(p,q,2);
    mpz_add_ui(p,p,1);
    //求生成元，1到p-1之间的随机数a不是单位元、阶不为2和q，
    //则其阶为p-1，进而证明得a为p的生成元
    mpz_t a_2,a_q;
    mpz_init(a_2);
    mpz_init(a_q);
    do{
        mpz_urandomm(g,gmp_state,p);//生成1到p-1的随机数
        mpz_powm_ui(a_2,g,2,p);//a_2=a^2 % p
        mpz_powm(a_q,g,q,p);//a_q=a^q % p
        if((mpz_cmp_ui(g,1)!=0)&&(mpz_cmp_ui(a_2,1)!=0)&&(mpz_cmp_ui(a_q,1)!=0))
            break;
    }while(1);

    do{
        mpz_urandomm(h,gmp_state,p);//生成1到p-1的随机数
        mpz_powm_ui(a_2,h,2,p);//a_2=a^2 % p
        mpz_powm(a_q,h,q,p);//a_q=a^q % p
        if((mpz_cmp_ui(h,1)!=0)&&(mpz_cmp_ui(a_2,1)!=0)&&(mpz_cmp_ui(a_q,1)!=0))
            break;
    }while(1);

    //mpz_urandomm(Xa,gmp_state,p);//生成私钥Xa
    //mpz_powm(Ya,a,Xa,p);//生成公钥Ya,Ya=a^Xa % p
*/

    //gen_prime(q, 512/2);
    mpz_set_ui(q, 13);
    mpz_mul_ui(p, q, 2);
    mpz_add_ui(p, p, 1);
    mpz_mul(n, p, q);
    //set g = 1+n
    //mpz_add_ui(g, q, 1);
    //mpz_mod(g, g, q);
    //mpz_sub_ui(h, q, 1);
    mpz_set_ui(g, 7);
    mpz_set_ui(h, 11);

    //???
    mpz_set(p, q);

    cout << "q: " << q << endl;
    cout << "p: " << p << endl;
    //cout << "n: " << n << endl;
    cout << "g: " << g << endl;
    cout << "h: " << h << endl;

    //Sender's messages: m1, ..., mn
    mpz_t messages[train_number];
    for(int i=0; i< train_number; i++){
        mpz_init(messages[i]);
        mpz_set_ui(messages[i], i+3);
        mpz_mod(messages[i], messages[i], q);
    }

    //Receiver's choices: sigma1, ..., sigmak
    int choice[k_value];
    mpz_t recv_messages[k_value];
    for(int i=0; i< k_value; i++) {
        //mpz_init(choice[i]);
        choice[i] = i+1;
        mpz_init(recv_messages[i]);
    }

    //Init//
    mpz_t tmp, tmp2, tmp3, rand;
    mpz_init(tmp);
    mpz_init(tmp2);
    mpz_init(tmp3);
    mpz_init(rand);

    mpz_t NA[3], NB[3];
    for(int i=0; i< 3; i++){
        mpz_init(NA[i]);
        mpz_init(NB[i]);
    }


    //Step1: R chooses two poly
    //choose random a_i in Zq
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    mpz_urandomm(NB[0], state, p);
    mpz_mul_ui(NB[0], NB[0], k_value);
    mpz_urandomm(NB[1], state, p);
    mpz_mul_ui(NB[2], NB[1], k_value);

    for(int i=0; i< 2; i++){
        mpz_urandomm(NA[i], state, p);
    }

    mpz_t MA;
    mpz_init(MA);

    return 0;
}
/*
int newPRShuffle_test(){

    int train_len = 10;
    cout << "Before perm: " << endl;
    mpz_t ds[train_len], d[train_len];
    for(int i=0; i< train_len; i++){
        mpz_init_set_ui(d[i], i);
        //mpz_init_set_ui(ds[i], 0);
        cout << d[i] << ", ";
    }
    cout << endl;
    //prand_perm(ds, d, op1, train_1024);
    gmp_randstate_t state_shuffle;
    gmp_randinit_mt(state_shuffle);
    mpz_t key_shuffle;
    mpz_init_set_ui(key_shuffle, 20220402);
    gmp_randseed(state_shuffle, key_shuffle);
    PRShuffle(d, train_len, state_shuffle);
    cout << "After perm1: " << endl;
    for(int i=0; i< train_len; i++){
        cout << d[i] << ", ";
    }
    cout << endl;

    for(int i=0; i< train_len; i++){
        mpz_set_ui(d[i], i);
    }
    mpz_init_set_ui(key_shuffle, 20220402);
    gmp_randseed(state_shuffle, key_shuffle);
    PRShuffle(d, train_len, state_shuffle);
    cout << "After perm2: " << endl;
    for(int i=0; i< train_len; i++){
        cout << d[i] << ", ";
    }
    cout << endl;

    for(int i=0; i< train_len; i++){
        mpz_set_ui(d[i], i);
    }
    mpz_init_set_ui(key_shuffle, 20220402);
    gmp_randseed(state_shuffle, key_shuffle);
    PRShuffle(d, train_len, state_shuffle);
    cout << "After perm3: " << endl;
    for(int i=0; i< train_len; i++){
        cout << d[i] << ", ";
    }
    cout << endl;

    return 0;
}
*/


int random_mult_test(){
    int distances_test[10] = {0, 8, 5, 4, 6, 7, 2, 9, 3, 1};
    mpz_t origin_messages[10];
    mpz_t origin_messages0[10], origin_messages1[10];
    mpz_t delta[10*9/2], delta0[10*9/2], delta1[10*9/2];
    mpz_t randseed, randnum;
    mpz_init(randseed);
    mpz_init(randnum);

    gmp_randstate_t state;
    gmp_randinit_mt(state);

    int index = 0;
    for(int i=0; i<10; i++){
        mpz_init(origin_messages[i]);
        mpz_set_ui(origin_messages[i], distances_test[i]);
        //sub share
        mpz_init(origin_messages0[i]);
        mpz_init(origin_messages1[i]);
        //rand
        mpz_urandomb(origin_messages0[i], state, message_bit_len);
        mpz_add(origin_messages1[i], origin_messages[i], origin_messages0[i]);

        for(int j=10-1; j>i; j--){
            mpz_init(delta0[index]);
            mpz_init(delta1[index]);
            index++;
        }
        //cout << i << endl;
    }

    cout << "Origin: " << endl;
    for(int i=0; i<10; i++){
        cout << origin_messages[i] << endl;
    }
    cout << endl;

    mpz_set_ui(randseed, 20220407);
    index = 0;

    gmp_randseed(state, randseed);

    for(int i=0; i<10; i++){
        for(int j=10-1; j>i; j--){
            mpz_sub(delta0[index], origin_messages0[i], origin_messages0[j]);
            mpz_sub(delta1[index], origin_messages1[i], origin_messages1[j]);
            mpz_urandomb(randnum, state, message_bit_len);
            mpz_mul(delta0[index], delta0[index], randnum);
            mpz_mul(delta1[index], delta1[index], randnum);
            index++;
        }
        //cout << i << endl;
    }
/*
    index = 0;
    for(int i=0; i<10; i++){
        for(int j=10-1; j>i; j--){
            cout << "[" << i << "][" << j << "]: " << delta[index] << "; ";
            index ++;
        }
        cout  << endl;
    }
*/

    signed long tmp_delta[10*9/2], tmp_delta0[10*9/2], tmp_delta1[10*9/2];
    bool tmp_cmp_result[10][10];
    index = 0;
    for(int i=0; i<10; i++){
        for(int j=10-1; j>i; j--){
            tmp_delta0[index] = mpz_get_si(delta0[index]);
            tmp_delta1[index] = mpz_get_si(delta1[index]);
            tmp_delta[index] = tmp_delta1[index] - tmp_delta0[index];


            if(tmp_delta[index] > 0){
                tmp_cmp_result[i][j] = true;
                //cmp_result[j][i] = false;
            }else if(tmp_delta[index] == 0){
                tmp_cmp_result[i][j] = false;
                //cmp_result[j][i] = false;
            }else{
                tmp_cmp_result[i][j] = false;
                //cmp_result[j][i] = true;
            }
            index++;
        }
        //cout << i << endl;
    }

    cout << "Delta0: " << endl;
    index = 0;
    for(int i=0; i<10; i++){
        for(int j=10-1; j>i; j--){
            cout << "[" << i << "][" << j << "]: " << tmp_delta0[index] << ": " << delta0[index] << "; ";
            index ++;
        }
        cout  << endl;
    }

    cout << "Delta1: " << endl;
    index = 0;
    for(int i=0; i<10; i++){
        for(int j=10-1; j>i; j--){
            cout << "[" << i << "][" << j << "]: " << tmp_delta1[index] << ": " << delta1[index] << "; ";
            index ++;
        }
        cout  << endl;
    }

//test tmp rm.
/*
    bool cmp_result[10][10];

    index = 0;
    for(int i=0; i<10; i++){
        for(int j=10-1; j>i; j--){
            mpz_init(delta[index]);
            mpz_sub(delta[index], delta1[index], delta0[index]);
            if(mpz_cmp_ui(delta[index], 0) > 0){
                cmp_result[i][j] = true;
                //cmp_result[j][i] = false;
            }else if(mpz_cmp_ui(delta[index], 0) == 0){
                cmp_result[i][j] = false;
                //cmp_result[j][i] = false;
            }else{
                cmp_result[i][j] = false;
                //cmp_result[j][i] = true;
            }
            index++;
        }
        //cout << i << endl;
    }

    index = 0;
    for(int i=0; i<10; i++){
        for(int j=10-1; j>i; j--){
            cout << "[" << i << "][" << j << "]: " << cmp_result[i][j] << "; ";
            index ++;
        }
        cout  << endl;
    }

*/


    int topk_result[3];
    int min=0, flag = 1;
    for (int k = 0; k < 3; k++){
        //while(flag){
            for(int j=min+1; j<10; j++){
                if(tmp_cmp_result[min][j]){
                    min = j;
                    j=min+1;
                    continue;
                }else{
                    continue;
                }
            }
            //cout << i << endl;
        //}
        //flag = 1;
        topk_result[k] = min;
        for(int j=min+1; j<10; j++){
            tmp_cmp_result[min][j] = true;
            //cmp_result[j][min] = true;
        }
        for(int j=0; j<min; j++){
            //cmp_result[min][j] = false;
            tmp_cmp_result[j][min] = false;
        }

        min = 0;

    }

    for(int i=0; i<3; i++){
        cout << topk_result[i] << endl;
    }


    //sorting...
    /*
    long result_values[train_number+100];
    long min = result_values[0];
    long mins[k_value];
    //mpz_init(min);
    //mpz_set(min, result_values[0]);
    result_values[0] = min;

    for(int j=0; j<k_value; j++){
        //long min;
        //mpz_init(min);
        int mini = 0;
        while((result_values[mini] == -1) && mini < train_number) mini++;
        min = result_values[mini];

        for(int i=0;i<train_number;i++){
            if((result_values[i]== -1)) continue;
            if((result_values[i]< min)){
                (min= result_values[i]);
                (result_values[i]= -1);
                continue;
            }
        }
        (mins[j]= min);
    }*/

/*
    index = 0;
    for(int i=0; i<10; i++){
        for(int j=10-1; j>i; j--){
            cout << "[" << i << "][" << j << "]: " << cmp_result[i][j] << "; ";
            index ++;
        }
        cout  << endl;
    }
*/
/*
    auto cmp = [](const int x, const int y) {
        if(x < y){
            return cmp_result[x][y];
        }else if(x == y){
            return false;
        }else{
            return cmp_result[y][x];
        }
    };
*/
/*
    vector<int> ind_of_labels;
    for(int i=0; i<10; i++){
        ind_of_labels.push_back(i);
    }

    vector<int> topk_ind_of_labels;

    if (ind_of_labels.empty() || 3 > ind_of_labels.size())
        cout << "Error!" << endl;

    for (int i = 0; i < 3; i++)
        topk_ind_of_labels.push_back(ind_of_labels[i]);

    //spec_sort * srt = new spec_sort();
    //srt->set_tree(cmp_result);
    make_heap(topk_ind_of_labels.begin(), topk_ind_of_labels.end(), spec_sort(cmp_result));


    for (int i = 3; i < ind_of_labels.size(); i++)
    {
        if (ind_of_labels[i] < topk_ind_of_labels[0])
        {
            pop_heap(topk_ind_of_labels.begin(), topk_ind_of_labels.end(), spec_sort(cmp_result));
            topk_ind_of_labels.pop_back();
            topk_ind_of_labels.push_back(ind_of_labels[i]);
            push_heap(topk_ind_of_labels.begin(), topk_ind_of_labels.end(), spec_sort(cmp_result));
        }
    }*/

    //make_heap(testv.begin(),testv.end(),F(help));

    /*

    for (int i = 0; i < 3; i++)
        topk_ind_of_labels.push_back(ind_of_labels[i]);

    make_heap(topk_ind_of_labels.begin(), topk_ind_of_labels.end(), cmp);


    for (int i = 3; i < ind_of_labels.size(); i++)
    {
        if (ind_of_labels[i] < topk_ind_of_labels[0])
        {
            pop_heap(topk_ind_of_labels.begin(), topk_ind_of_labels.end(), cmp);
            topk_ind_of_labels.pop_back();
            topk_ind_of_labels.push_back(ind_of_labels[i]);
            push_heap(topk_ind_of_labels.begin(), topk_ind_of_labels.end(), cmp);
        }
    }
    */


/*
    for (int i = 0; i < 3; i++)
    {
        pop_heap(ind_of_labels.begin(), ind_of_labels.end(), cmp);
        topk_ind_of_labels.push_back(ind_of_labels.back());
        ind_of_labels.pop_back();
    }
*//*
    for(int i=0; i<3; i++){
        cout << topk_ind_of_labels[i] << endl;
    }

*/

    //cout << "Ok" << endl;
    return 0;
}

int main(int argc, char** argv) {

    //k_out_of_n_OT_test();
    //newPRShuffle_test();
    //return 0;

    //random_mult_test();
    //return 0;
    //HSS_MULT_TEST();
    //return 0;
    //HSS_test();
    //return 0;
    //test_OT();
    //return 0;

    //input
    //int port, num_iters;
    string address;

    PARTY = atoi(argv[1]);
    //port = atoi(argv[2]); //unuse
    server_port = atoi(argv[2]);


    //PARTY = 1;
    //server_port = 10240;


    try{
        int x = -1;
        if(argc <= 3)
            throw x;
        address = argv[3];
    }catch(int x) {
        address = "127.0.0.1";
    }

    int mnist;
    //return comm_lr_b_test();
    //return ftp_test();
    //mnist = knn_mnist_formal();

    if(PARTY == 1) server_id = 0;
    else server_id = 1;

    //test_OT();
    //return 0;


    mnist = server_evaluate();
    return 0;
}

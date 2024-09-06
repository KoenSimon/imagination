//
// Created by kore on 9/11/21.
//

#include "MyTime.h"

using namespace std;


int MyTime::TimeRecorded(std::string info, char const * file_path){
    ofstream history_outfile;

    // 基于当前系统的当前日期/时间
    time_t now = time(0);

    // 把 now 转换为字符串形式
    char* dt = ctime(&now);

    std::cout << info << std::endl;
    std::cout << "Time：" << dt;


    // 以写模式打开文件
    history_outfile.open(file_path, std::ios::app);

    // 向文件写入用户输入的数据
    history_outfile << "Message: " << info << std::endl;
    history_outfile << "Time： " << dt;
    //history_outfile << "UTC 日期和时间："<< dt;






    clock_t dur;

    dur = clock();
    dur = dur - this_second;
    this_second = clock();

    //cout << "Server Time：" << dt << endl;

    cout<< "This phase took "<< double(dur)/CLOCKS_PER_SEC<<" s."<<endl << endl;  //输出时间（单位：ｓ）
    history_outfile << "This phase took " <<double(dur)/CLOCKS_PER_SEC<<" s." << dt;



    // 关闭打开的文件
    history_outfile.close();





    return 0;
}

int MyTime::TimeRecorded(std::string info, char const * file_path, int type){
    ofstream history_outfile;

    // 基于当前系统的当前日期/时间
    time_t now = time(0);

    // 把 now 转换为字符串形式
    char* dt = ctime(&now);

    std::cout << info << std::endl;
    std::cout << "Time：" << dt;


    // 以写模式打开文件
    history_outfile.open(file_path, std::ios::app);

    // 向文件写入用户输入的数据
    history_outfile << "Message: " << info << std::endl;
    history_outfile << "Time： " << dt;
    //history_outfile << "UTC 日期和时间："<< dt;





/*
    clock_t dur;

    dur = clock();
    dur = dur - this_second;
    this_second = clock();

    //cout << "Server Time：" << dt << endl;
    string phase_info;
    switch(type){
        case 0:
            clients_time += dur;
            phase_info = "Runtime of clients is ";
            break;
        case 10:
            server_time_0 += dur;
            phase_info = "Runtime of server 0 is ";
            break;
        case 11:
            server_time_1 += dur;
            phase_info = "Runtime of server 1 is ";
            break;
        case 12:
            server_time_0 += dur/2;
            server_time_1 += dur/2;
            break;
        case 2:
            controller_time += dur;
            phase_info = "Runtime of AC is ";
            break;
    }


    cout<< phase_info << double(dur)/CLOCKS_PER_SEC<<" s."<<endl << endl;  //输出时间（单位：ｓ）
    history_outfile << phase_info <<double(dur)/CLOCKS_PER_SEC<<" s." << dt;

*/

    // 关闭打开的文件
    history_outfile.close();

    return 0;
}

int MyTime::Time_begin(int role){
    //time_t this_time, time_distance;
    this_second = clock();
/*
    then = now;
    now = this_time;*/

    return 0;
}

int MyTime::Time_end(int role){
    clock_t another_second, time_distance;
    another_second = clock();

    //time_t this_time, time_distance;
    //this_time = time(0);

    time_distance = another_second - this_second;
/*
    switch(role){
        case 0:
            clients_time += time_distance;
            break;
        case 10:
            server_time_0 += time_distance;
            break;
        case 11:
            server_time_1 += time_distance;
            break;
        case 12:
            server_time_0 += time_distance/2;
            server_time_1 += time_distance/2;
            break;
        case 2:
            controller_time += time_distance;
            break;
    }

    //then = this->now;
    //now = this_time;
    this_second = another_second;*/

    return 0;
}

int MyTime::Time_report() {
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    //double server_time = this->server_time_0;
    //cout << server_time_0 << endl;
    //cout<<"Clients time = "<< to_string(static_cast<double>(this->clients_time/CLOCKS_PER_SEC)) <<" s"<<endl;  //输出时间（单位：ｓ）
    //cout<<"Server 0 time = "<< to_string(static_cast<double>(this->server_time_0)/CLOCKS_PER_SEC) <<" s"<<endl;  //输出时间（单位：ｓ）
    //cout<<"Server 1 time = "<< to_string(static_cast<double>(this->server_time_1)/CLOCKS_PER_SEC) <<" s"<<endl;  //输出时间（单位：ｓ）
    //cout<<"AC time = "<< to_string(static_cast<double>(this->controller_time)/CLOCKS_PER_SEC) <<" s"<<endl;  //输出时间（单位：ｓ）
    return 0;
}

//
// Created by kore on 9/11/21.
//

#ifndef PPKNN_MYTIME_H
#define PPKNN_MYTIME_H

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>

class MyTime {
public:

    clock_t this_second;//, server_time_0, server_time_1, clients_time, controller_time;

    //char * time_s = new char[1024];

    MyTime(){
        //time_t now;
        this->this_second = clock();
        /*
        this->clients_time = clock();
        this->server_time_0 = clock();
        this->server_time_1 = clock();
        this->controller_time = clock();
        time_s = ctime(&(now));
        */
    }



    int TimeRecorded(std::string info, char const * file_path);

    int TimeRecorded(std::string info, char const * file_path, int type);

    int Time_begin(int role);

    int Time_end(int role);

    int Time_report();

};


#endif //PPKNN_MYTIME_H

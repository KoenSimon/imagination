//
// Created by kore on 4/8/22.
//

#ifndef KNN_SERVER_SPEC_SORT_H
#define KNN_SERVER_SPEC_SORT_H


#include<iostream>
#include<vector>
#include<algorithm>
#include <queue>
#include <functional>
#include "crypto.h"
using namespace std;

class spec_sort{

public:
    bool sort_tree[10][10];
    spec_sort(bool in_sort_tree[10][10]){
        for(int i=0; i< 10; i++){
            for(int j=0; j< 10; j++){
                sort_tree[i][j] = in_sort_tree[i][j];
            }
        }
    }
    bool operator()(const int& a,const int& b){
        if(a < b){
            return !sort_tree[a][b];
        }else if(a == b){
            return false;
        }else{
            return !sort_tree[b][a];
        }

        if(sort_tree[a][b] > 0)
            return true;
        else
            return false;
    }
};


#endif //KNN_SERVER_SPEC_SORT_H

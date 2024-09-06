//
// Created by kore on 2/25/22.
//

#ifndef KNN_SERVER_POLY_H
#define KNN_SERVER_POLY_H
using namespace std;

class Poly {
public:
    Poly();		// 构造函数
    ~Poly();    // 析构函数
    void initParams(int n,int *indices,double *coef);//初始化
    void readPoly(char * p); //读取多项式
    void printPoly();//输出多项式
    void Simplification();//化简多项式
    Poly * togetherPoly(const Poly  &q);//合并多项式
    Poly * addition(const Poly  &q);//加法
    Poly * Subtraction(const Poly  &q);//减法
    Poly * multiplication(const Poly  &q);//乘法
    Poly * algorithm( const Poly &q, char op);//运算规则
    void clear();
private:
    int m_num;
    double *m_coefficient;
    int *m_index;
    int pri;
};


#endif //KNN_SERVER_POLY_H

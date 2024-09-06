//
// Created by kore on 2/25/22.
//
#include <fstream>
#include <iostream>
#include "Poly.h"

Poly::Poly()
{
    m_num = 0;
    m_coefficient = NULL;
    m_index = NULL;
    pri=0;
}
Poly::~Poly()
{
    //printf("hello\n");
    clear();
}
void   Poly::initParams(int n, int *indices, double *coef)
{
    m_num = n;
    m_index = indices;
    m_coefficient = coef;
}
void   Poly::readPoly(char * p)
{

    std::ifstream f(p);

    if (!f.is_open())
    {
        printf("Error to Read \n");
        exit(0);
    }

    f >> m_num;
    double *pa1 = new double[m_num];
    int *pa2 = new int[m_num];
    for (int i = 0; i < m_num; i++)
    {
        f >> pa1[i];
        f >> pa2[i];
    }
    f.close();
    m_coefficient = pa1;
    m_index = pa2;

}
void   Poly::printPoly()
{
    printf("\n");
    for (int j = 0; j < m_num; j++)
    {
        if ((m_coefficient[j]) != 0)
        {
            printf("(");
            printf("%.4lf", m_coefficient[j]);
            printf(")");
            printf("*x^");
            if (j < ((m_num)-1))
            {
                printf("%d+", m_index[j]);
            }
            else
            {
                printf("%d", m_index[j]);
            }
        }
    }
    pri++;
    printf("\n");
}
void   Poly::Simplification()
{
    double chang1;
    int chang2;
    //排序
    for (int i = 0; i < m_num; i++)
    {
        for (int j = i; j<m_num; j++)
        {
            if (m_index[i]>m_index[j])
            {
                chang2 = m_index[j]; m_index[j] = m_index[i]; m_index[i] = chang2;
                chang1 = m_coefficient[j]; m_coefficient[j] = m_coefficient[i]; m_coefficient[i] = chang1;
            }
        }
    }
    int m = 0;                                  //合并指数相同的项
    for (int i = 0; i < m_num; i++)
    {
        m_coefficient[m] = m_coefficient[i];
        m_index[m] = m_index[i];
        for (int j = i + 1; j < m_num; j++)
        {
            if (m_index[j] == m_index[i])
            {
                m_coefficient[m] = m_coefficient[m] + m_coefficient[j];
                i = j;
            }
        }
        m++;
    }
    m_num = m;

}
Poly * Poly::togetherPoly(const Poly  &q)
{

    double *pa1 = new double[m_num + q.m_num];
    int *pa2 = new int[m_num + q.m_num];

    for (int j = 0; j < (m_num); j++)
    {
        pa1[j] = m_coefficient[j];

        pa2[j] = m_index[j];
    }

    for (int j = 0; j < (q.m_num); j++)
    {
        pa1[j+m_num] = q.m_coefficient[j];

        pa2[j + m_num] = q.m_index[j];
    }

    Poly *z = new Poly;

    z->initParams(m_num + q.m_num, pa2, pa1);

    return z;
}
Poly * Poly::addition(const Poly  &q)
{

    //Step1: Put p q together
    Poly *res1 = Poly::togetherPoly(q);

    //Step2: Simplification

    res1->Simplification();

    return res1;
}
Poly * Poly::Subtraction(const Poly  &q)
{
    //Step1: Put p q together

    for (int i = 0; i < q.m_num; i++)
    {
        q.m_coefficient[i] = q.m_coefficient[i] * (-1);
    }
    Poly * res1 = Poly::togetherPoly(q);

    //Step2: Simplification

    res1->Simplification();

    return res1;
}
Poly * Poly::multiplication(const Poly  &q)
{

    double *pa1 = new double[(m_num)*(q.m_num)];
    int *pa2 = new int[(m_num)*(q.m_num)];
    int m = 0;
    for (int i = 0; i < m_num; i++)
    {
        for (int j = 0; j < q.m_num; j++)
        {
            pa2[m] = (m_index[i]) + (q.m_index[j]);
            pa1[m] = (m_coefficient[i])*(q.m_coefficient[j]);
            m++;
        }
    }

    Poly *z = new Poly;

    z->initParams(m_num *q.m_num, pa2, pa1);


    z->Simplification();

    return z;
}

Poly * Poly::algorithm(const Poly  &q, char op)
{

    if ('+' == op)
    {
        Poly *res = Poly::addition(q);
        printf("\n");
        return res;
    }
    else if ('-' == op)
    {

        Poly *res = Poly::Subtraction(q);
        printf("\n");
        return res;
    }
    else if ('*' == op)
    {
        Poly *res = Poly::multiplication(q);
        printf("\n");
        return res;
    }
    return NULL;
}
void   Poly::clear()
{
    if (m_coefficient != NULL)
    {
        delete[] m_coefficient;
    }
    if (m_index != NULL)
    {
        delete[] m_index;
    }
    m_num = 0;
    m_coefficient = NULL;
    m_index = NULL;
}

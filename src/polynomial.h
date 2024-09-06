//
// Created by kore on 3/23/22.
//

#ifndef KNN_SERVER_POLYNOMIAL_H
#define KNN_SERVER_POLYNOMIAL_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <gmp.h>
#include <gmpxx.h>

class polynomial
{
private:
    std::vector<mpz_class> coefficient;
    std::vector<int> radix;
public:
    polynomial(){}
    polynomial(int K, std::vector<int> r, std::vector<mpz_class> co){
        for (int i = 0; i<K; i++)
        {
            mpz_class tmpcoe;
            int tmprad;
            tmprad = r[i];
            tmpcoe = co[i];
            this->radix.push_back(tmprad);
            this->coefficient.push_back(tmpcoe);
        }
    }
    polynomial(const polynomial &poly){ *this = poly; }
    ~polynomial(){}
    //int output(mpz_t Kn, std::vector<mpz_t> * rad, std::vector<mpz_t> * coe);
    int output(int Kn, int rad[], mpz_class coe[]);
    int evaluate(mpz_t result, mpz_t input, mpz_t q);
    polynomial sort(mpz_class params[]);
    friend polynomial operator+(const polynomial &poly1, const polynomial &poly2);
    friend polynomial operator-(const polynomial &poly1, const polynomial &poly2);
    polynomial& operator+=(const polynomial &poly){ *this = *this + poly; return *this; }
    polynomial& operator-=(const polynomial &poly){ *this = *this - poly; return *this; }
    friend polynomial operator*(const polynomial &poly1, const polynomial &poly2);
    friend std::istream& operator>>(std::istream &is, polynomial &poly);
    friend std::ostream& operator<<(std::ostream &os, const polynomial &poly);
};
/*
int polynomial::output(int Kn, int rad[], int coe[]) {
    //mpz_t K;
    mpz_init(Kn);
    if (this->radix.size() != 0)
        mpz_set_ui(Kn, this->radix.size());

    for (int i = 0; i<this->radix.size(); i++)
    {
        mpz_t co;
        mpz_t ra;
        mpz_init(co);
        mpz_init(ra);
        mpz_set_ui(co, this->coefficient[i]);
        mpz_set_ui(ra, this->radix[i]);
        rad->push_back(ra);
        coe->push_back(co);
    }
    return 0;
}
*/

int polynomial::evaluate(mpz_t result, mpz_t input, mpz_t q){
    mpz_class x, y, tmp;
    mpz_t x2;
    mpz_init(x2);

    mpz_t modu;
    mpz_init(modu);
    mpz_sub_ui(modu, q, 1);
    mpz_set(modu, q);

    int inputs, exps;
    mpz_set_ui(result, 0);

    for (int i = 0; i<this->radix.size(); i++)
    {
        mpz_powm_ui(x2, input, this->radix[i], modu);
        //inputs = mpz_get_ui(input);
        //exps = pow(inputs, this->radix[i]);
        //mpz_set_ui();

        tmp = mpz_class(x2);//x2
        tmp *= this->coefficient[i];
        mpz_add(result, result, tmp.get_mpz_t());
    }
    //mpz_mod(result, result, q);

    mpz_clear(x2);
    return 0;
}

int polynomial::output(int Kn, int rad[], mpz_class coe[]) {
    if (this->radix.size() != 0)
        Kn = this->radix.size();

    for (int i = 0; i<this->radix.size(); i++)
    {
        coe[i] = this->coefficient[i];
        rad[i] = this->radix[i];
    }
    return 0;
}

polynomial polynomial::sort(mpz_class params[]) {
    int tmp_coe;
    polynomial tmp_poly;
    for (int i = 0; i<this->radix.size(); i++)
    {
        tmp_coe = this->radix[i];
        for (int j = i+1; j<this->radix.size(); j++){
            if(tmp_coe == this->radix[j]){
                this->coefficient[i] +=  this->coefficient[j];
                this->coefficient[j] = 0;
            }
        }
        if(this->coefficient[i] != 0){
            tmp_poly.radix.push_back(this->radix[i]);
            tmp_poly.coefficient.push_back(this->coefficient[i]);
            params[this->radix[i]] = this->coefficient[i];
        }
    }
    return tmp_poly;
}

polynomial operator+(const polynomial &poly1, const polynomial &poly2)
{
    polynomial polysum;

    int i = 0, j = 0;
    while (i<poly1.radix.size() || j<poly2.radix.size())
    {
        //only poly1
        if (i<poly1.radix.size() && j == poly2.radix.size())
        {
            polysum.radix.push_back(poly1.radix[i]);
            polysum.coefficient.push_back(poly1.coefficient[i]);
            i++;
        }
            //only poly2
        else if (j<poly2.radix.size() && i == poly1.radix.size())
        {
            polysum.radix.push_back(poly2.radix[j]);
            polysum.coefficient.push_back(poly2.coefficient[j]);
            j++;
        }
            //the radix of poly1 greater than poly2
        else if (poly1.radix[i] > poly2.radix[j])
        {
            polysum.radix.push_back(poly1.radix[i]);
            polysum.coefficient.push_back(poly1.coefficient[i]);
            i++;
        }
            //the radix of poly1 smaller than poly2
        else if (poly1.radix[i] < poly2.radix[j])
        {
            polysum.radix.push_back(poly2.radix[j]);
            polysum.coefficient.push_back(poly2.coefficient[j]);
            j++;
        }
            //the radix of poly1 equal to poly2
        else
        {
            if (poly1.coefficient[i] + poly2.coefficient[j] != 0)
            {
                polysum.radix.push_back(poly1.radix[i]);
                polysum.coefficient.push_back(poly1.coefficient[i] + poly2.coefficient[j]);
            }
            i++;
            j++;
        }
    }

    return polysum;
}

polynomial operator-(const polynomial &poly1, const polynomial &poly2)
{
    polynomial negativePoly;

    negativePoly = poly2;
    for (int i = 0; i < negativePoly.coefficient.size(); i++)
        negativePoly.coefficient[i] = -negativePoly.coefficient[i];

    return poly1 + negativePoly;
}

polynomial operator*(const polynomial &poly1, const polynomial &poly2)
{
    polynomial mul;

    int i = 0;
    while (i < poly2.coefficient.size())
    {
        polynomial part = poly1;
        mpz_class coefficient = poly2.coefficient[i];
        int radix = poly2.radix[i];

        for (int j = 0; j < part.coefficient.size(); j++)
        {
            part.radix[j] += radix;
            part.coefficient[j] *= coefficient;
        }
        mul += part;
        i++;
    }

    return mul;
}

std::istream& operator>>(std::istream &is, polynomial &poly)
{
    int k;
    is >> k;
    for (int i = 0; i<k; i++)
    {
        double tmpcoe;
        int tmprad;
        is >> tmprad >> tmpcoe;
        poly.radix.push_back(tmprad);
        poly.coefficient.push_back(tmpcoe);
    }

    return is;
}

std::ostream& operator<<(std::ostream &os, const polynomial &poly)
{
    os << poly.radix.size();
    if (poly.radix.size() != 0)
        std::cout << " ";
    for (int i = 0; i<poly.radix.size(); i++)
    {
        os << poly.radix[i] << " " << std::fixed << std::setprecision(1) << poly.coefficient[i];
        if (i != poly.radix.size() - 1)
            os << " ";
    }

    return os;
}


#endif //KNN_SERVER_POLYNOMIAL_H

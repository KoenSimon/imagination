//
// Created by kore on 9/8/21.
//

#include <iostream>
#include <cmath>
#include "crypto.h"

using namespace std;

//paillier parameters
// where already have #define KeyBitLength 64 and maxSizeOfDigits 64

//mpz_t base_sk, base_msg;
//int sizeOfDigits;



inline void debug_msg(const char *str) {
#ifdef PAILLIER_DEBUG
    fputs(str, stderr);
#endif
}

int paillier_ell(mpz_t result, mpz_t input, mpz_t ninv, mp_bitcnt_t len) {
    mpz_t mask;

    mpz_init(mask);

    mpz_sub_ui(result, input, 1);
    mpz_mul(result, result, ninv);
    mpz_setbit(mask, len);
    mpz_sub_ui(mask, mask, 1);
    mpz_and(result, result, mask);

    mpz_clear(mask);
    return 0;
}
int gen_random(mpz_t rnd, mp_bitcnt_t len) {
    FILE *dev_random, *dev_urandom;
    int byte_count, byte_read;
    char * seed;

    byte_count = BIT2BYTE(len);

    dev_random = fopen("/dev/random", "r");
    if(dev_random == NULL) {
        fprintf(stderr, "cannot open random number device!\n");
        exit(1);
    }
    dev_urandom = fopen("/dev/urandom", "r");
    if(dev_urandom == NULL) {
        fprintf(stderr, "cannot open random number device!\n");
        exit(1);
    }

    seed = (char *)malloc(sizeof(char)*byte_count);

    byte_read = 0;
    //generate the first 16 bytes with /dev/random
    while(byte_read  < 16 && byte_read < byte_count) {
        byte_read += fread(seed, sizeof(char), byte_count, dev_random);
    }
    fclose(dev_random);
    //generate the remaining bytes with /dev/urandom
    while(byte_read < byte_count) {
        byte_read += fread(seed, sizeof(char), byte_count, dev_urandom);
    }
    fclose(dev_urandom);

    mpz_import(rnd, byte_count, 1, sizeof(seed[0]), 0, 0, seed);
    free(seed);
    return 0;
}

int gen_prime(mpz_t prime, mp_bitcnt_t len) {
    mpz_t rnd;

    mpz_init(rnd);

    gen_random(rnd, len);

    //set most significant bit to 1
    mpz_setbit(rnd, len-1);
    //look for next prime
    mpz_nextprime(prime, rnd);

    mpz_clear(rnd);
    return 0;
}

int crt_exponentiation(mpz_t result, mpz_t base, mpz_t exp_p, mpz_t exp_q, mpz_t pinvq, mpz_t p, mpz_t q) {
    mpz_t pq;
    exp_args *args_p, *args_q;

#ifdef PAILLIER_THREAD
    pthread_t thread1, thread2;
#endif

    mpz_init(pq);

    //prepare arguments for exponentiation mod p
    args_p = (exp_args *)malloc(sizeof(exp_args));

    mpz_init(args_p->result);
    mpz_init(args_p->basis);
    mpz_init(args_p->exponent);
    mpz_init(args_p->modulus);

    mpz_set(args_p->basis, base);
    mpz_set(args_p->exponent, exp_p);
    mpz_set(args_p->modulus, p);

    //prepare arguments for exponentiation mod q
    args_q = (exp_args *)malloc(sizeof(exp_args));

    mpz_init(args_q->result);
    mpz_init(args_q->basis);
    mpz_init(args_q->exponent);
    mpz_init(args_q->modulus);

    mpz_set(args_q->basis, base);
    mpz_set(args_q->exponent, exp_q);
    mpz_set(args_q->modulus, q);

#ifdef PAILLIER_THREAD
    //compute exponentiation modulo p
	pthread_create(&thread1, NULL, do_exponentiate, (void *)args_p);

	//compute exponentiation modulo q
	pthread_create(&thread2, NULL, do_exponentiate, (void *)args_q);

	pthread_join(thread1, NULL);
	pthread_join(thread2, NULL);

#else
    //compute exponentiation modulo p
    mpz_mod(args_p->result, base, p);
    mpz_powm(args_p->result, args_p->result, exp_p, p);

    //compute exponentiation modulo q
    mpz_mod(args_q->result, base, q);
    mpz_powm(args_q->result, args_q->result, exp_q, q);
#endif

    //recombination
    mpz_mul(pq, p, q);
    mpz_sub(result, args_q->result, args_p->result);
    mpz_mul(result, result, p);
    mpz_mul(result, result, pinvq);
    mpz_add(result, result, args_p->result);
    mpz_mod(result, result, pq);

    mpz_clear(pq);
    mpz_clear(args_p->result);
    mpz_clear(args_p->basis);
    mpz_clear(args_p->exponent);
    mpz_clear(args_p->modulus);
    mpz_clear(args_q->result);
    mpz_clear(args_q->basis);
    mpz_clear(args_q->exponent);
    mpz_clear(args_q->modulus);
    free(args_p);
    free(args_q);

    return 0;
}

int crt_findD(paillier_private_key_v *priv){
    mpz_t tmp1;
    mpz_t tmp2;
    mpz_t tmp3;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(tmp3);

    mpz_gcd(tmp1, priv->lambda, priv->n);
    if(mpz_cmp_ui(tmp1, 1) != 0) return 0;

    //mpz_invert(tmp1, priv->n, priv->lambda);

    mpz_invert(tmp2, priv->lambda, priv->n);

    mpz_mul(tmp3, priv->lambda, priv->n);

    mpz_mul(priv->d, priv->lambda, tmp2);
    mpz_mod(priv->d, priv->d, tmp3);

    return 1;
}

int paillier_keygen_variant(paillier_public_key *pub, paillier_private_key_v *priv, mp_bitcnt_t len) {
    mpz_t p, q, n2, temp, mask, g;

    mpz_t tmp1;
    mpz_t tmp2;
    mpz_t tmp3;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(tmp3);

    mpz_init(p);
    mpz_init(q);
    mpz_init(n2);
    mpz_init(temp);
    mpz_init(mask);
    mpz_init(g);

    //write bit lengths
    priv->len = len;
    pub->len = len;

    //generate p and q
    debug_msg("generating prime p\n");
    gen_prime(p, len/2);
    cout << "p: " << p << endl;
    debug_msg("generating prime q\n");
    gen_prime(q, len/2);
    cout << "q: " << q << endl;

    //calculate modulus n=p*q
    debug_msg("calculating modulus n=p*q\n");
    mpz_mul(pub->n, p, q);
    mpz_mul(priv->n, p, q);
    cout << "n: " << pub->n << endl;

    //set g = 1+n
    debug_msg("calculating basis g=1+n\n");
    mpz_add_ui(g, pub->n, 1);

    //compute n^{-1} mod 2^{len}
    debug_msg("computing modular inverse n^{-1} mod 2^{len}\n");
    mpz_setbit(temp, len);
    if(!mpz_invert(priv->ninv, pub->n, temp)) {
        fputs("Inverse does not exist!\n", stderr);
        mpz_clear(p);
        mpz_clear(q);
        mpz_clear(n2);
        mpz_clear(temp);
        mpz_clear(mask);
        mpz_clear(g);
        exit(1);
    }

    //compute p^2 and q^2
    mpz_mul(priv->p2, p, p);
    mpz_mul(priv->q2, q, q);

    //generate CRT parameter
    debug_msg("calculating CRT parameter p^{-2} mod q^2\n");
    mpz_invert(priv->p2invq2, priv->p2, priv->q2);

    //calculate lambda = lcm(p-1,q-1)
    debug_msg("calculating lambda=lcm((p-1)*(q-1))\n");
    mpz_clrbit(p, 0);
    mpz_clrbit(q, 0);
    mpz_lcm(priv->lambda, p, q);

    //calculate n^2
    mpz_mul(n2, pub->n, pub->n);

    //calculate mu
    debug_msg("calculating mu\n");
    crt_exponentiation(temp, g, priv->lambda, priv->lambda, priv->p2invq2, priv->p2, priv->q2);

    paillier_ell(temp, temp, priv->ninv, len);

    if(!mpz_invert(priv->mu, temp, pub->n)) {
        fputs("Inverse does not exist!\n", stderr);
        mpz_clear(p);
        mpz_clear(q);
        mpz_clear(n2);
        mpz_clear(temp);
        mpz_clear(mask);
        mpz_clear(g);
        exit(1);
    }
/*
  mpz_t tmp1;
  mpz_init(tmp1);

  if(!mpz_invert(tmp1, priv->lambda, pub->n)) {
		fputs("Inverse does not exist!\n", stderr);
		mpz_clear(p);
		mpz_clear(q);
		mpz_clear(n2);
		mpz_clear(temp);
		mpz_clear(mask);
		mpz_clear(g);
		exit(1);
	}
*/


    //mpz_add(priv->d, priv->lambda, priv->mu);
    //mpz_sub_ui(priv->d, priv->d, 2);

    crt_findD(priv);

    //mpz_lcm(tmp1, priv->lambda, priv->n);

    //mpz_mul();


    //mpz_sub(priv->d, priv->lambda, priv->mu);
    //mpz_mod(priv->d, priv->d, priv->n);

    //free memory and exit
    debug_msg("freeing memory\n");
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(n2);
    mpz_clear(temp);
    mpz_clear(mask);
    mpz_clear(g);
    debug_msg("exiting\n");
    return 0;
}

int gen_pseudorandom(mpz_t rnd, mp_bitcnt_t len) {
    FILE *dev_urandom;
    int byte_count, byte_read;
    char * seed;

    byte_count = BIT2BYTE(len);

    dev_urandom = fopen("/dev/urandom", "r");
    if(dev_urandom == NULL) {
        fprintf(stderr, "cannot open random number device!\n");
        exit(1);
    }

    seed = (char *)malloc(sizeof(char)*byte_count);

    byte_read = 0;
    //generate the bytes with /dev/urandom
    while(byte_read < byte_count) {
        byte_read += fread(seed, sizeof(char), byte_count, dev_urandom);
    }
    fclose(dev_urandom);

    mpz_import(rnd, byte_count, 1, sizeof(seed[0]), 0, 0, seed);
    free(seed);
    return 0;
}

int paillier_encrypt(mpz_t ciphertext, mpz_t plaintext, paillier_public_key *pub) {
    mpz_t n2, r;
    //cout << pub->n << endl;
    if(mpz_cmp(pub->n, plaintext)) {
        //cout << "entering." << endl;
        mpz_init(n2);
        mpz_init(r);

        //re-compute n^2

        mpz_mul(n2, pub->n, pub->n);

        debug_msg("generating random number\n");
        //cout << "generating random number" << endl;
        //generate random r and reduce modulo n
        gen_pseudorandom(r, pub->len);
        mpz_mod(r, r, pub->n);
        if(mpz_cmp_ui(r, 0) == 0) {
            fputs("random number is zero!\n", stderr);
            mpz_clear(n2);
            mpz_clear(r);
            exit(1);
        }

        debug_msg("computing ciphertext\n");
        //cout << "computing ciphertext" << endl;
        //compute r^n mod n2
        mpz_powm(ciphertext, r, pub->n, n2);

        //compute (1+m*n)
        mpz_mul(r, plaintext, pub->n);
        mpz_add_ui(r, r, 1);

        //multiply with (1+m*n)
        mpz_mul(ciphertext, ciphertext, r);
        mpz_mod(ciphertext, ciphertext, n2);

        debug_msg("freeing memory\n");
        mpz_clear(n2);
        mpz_clear(r);
    }
    debug_msg("exiting\n");
    return 0;
}

int paillier_decrypt_v(mpz_t plaintext, mpz_t ciphertext, paillier_private_key_v *priv) {
    debug_msg("computing plaintext\n");
    //compute exponentiation c^lambda mod n^2
    crt_exponentiation(plaintext, ciphertext, priv->d, priv->d, priv->p2invq2, priv->p2, priv->q2);

    //compute exponentiation c^lambda mod n^2
    //compute n2
    //mpz_powm(plaintext, ciphertext, priv->d, priv->n2);

    //compute L(c^lambda mod n^2)
    paillier_ell(plaintext, plaintext, priv->ninv, priv->len);

    //compute L(c^lambda mod n^2)*mu mod n
    //mpz_mul(plaintext, plaintext, priv->mu);
    mpz_mod(plaintext, plaintext, priv->n);

    debug_msg("exiting\n");
    return 0;
}

/**
 * "Add" two plaintexts homomorphically by multiplying ciphertexts modulo n^2.
 * For example, given the ciphertexts c1 and c2, encryptions of plaintexts m1 and m2,
 * the value c3=c1*c2 mod n^2 is a ciphertext that decrypts to m1+m2 mod n.
 */
int paillier_homomorphic_add(mpz_t ciphertext3, mpz_t ciphertext1, mpz_t ciphertext2, paillier_public_key *pub) {
    mpz_t n2;

    mpz_init(n2);
    debug_msg("compute n^2");
    mpz_mul(n2, pub->n, pub->n);

    debug_msg("homomorphic add plaintexts");
    mpz_mul(ciphertext3, ciphertext1, ciphertext2);
    mpz_mod(ciphertext3, ciphertext3, n2);

    debug_msg("freeing memory\n");
    mpz_clear(n2);
    debug_msg("exiting\n");
    return 0;
}

/**
 * "Multiplies" a plaintext with a constant homomorphically by exponentiating the ciphertext modulo n^2 with the constant as exponent.
 * For example, given the ciphertext c, encryptions of plaintext m, and the constant 5,
 * the value c3=c^5 n^2 is a ciphertext that decrypts to 5*m mod n.
 */
int paillier_homomorphic_multc(mpz_t ciphertext2, mpz_t ciphertext1, mpz_t constant, paillier_public_key *pub) {
    mpz_t n2;

    mpz_init(n2);
    debug_msg("compute n^2");
    mpz_mul(n2, pub->n, pub->n);

    debug_msg("homomorphic multiplies plaintext with constant");
    mpz_powm(ciphertext2, ciphertext1, constant, n2);

    debug_msg("freeing memory\n");
    mpz_clear(n2);
    debug_msg("exiting\n");
    return 0;
}

int paillier_homomorphic_multc_random(mpz_t ciphertext2, mpz_t ciphertext1, mpz_t constant, mpz_t random_r, paillier_public_key* pub) {


    mpz_t n2;
    mpz_t rn;

    mpz_init(n2);
    mpz_init(rn);

    debug_msg("compute n^2");
    mpz_mul(n2, pub->n, pub->n);

    debug_msg("homomorphic multiplies plaintext with constant");
    mpz_powm(rn, random_r, pub->n, n2);

    mpz_powm(ciphertext2, ciphertext1, constant, n2);
    mpz_mul(ciphertext2, ciphertext2, rn);
    mpz_mod(ciphertext2, ciphertext2, n2);

    debug_msg("freeing memory\n");
    mpz_clear(n2);
    mpz_clear(rn);
    debug_msg("exiting\n");
    return 0;
}

int DDLog(mpz_t distance, mpz_t ciphertext1, paillier_public_key *pub) {
    mpz_t n2;
    mpz_t h;
    mpz_t hs;

    mpz_init(n2);
    mpz_init(h);
    mpz_init(hs);


    debug_msg("compute n^2");
    mpz_mul(n2, pub->n, pub->n);

    debug_msg("discompose g=h+hs*N");
    //mpz_add_ui(g, pub->n, 1);
    mpz_fdiv_qr(hs, h, ciphertext1, pub->n);
    //cout << h << endl;

    debug_msg("compute z=hs*h^-1");
    mpz_invert(distance, h, pub->n);
    mpz_mul(distance, distance, hs);
    mpz_mod(distance, distance, pub->n);
    //cout << distance << endl;


    debug_msg("freeing memory\n");
    mpz_clear(n2);
    mpz_clear(h);
    mpz_clear(hs);
    debug_msg("exiting\n");
    return 0;
}

void paillier_public_init(paillier_public_key *pub) {
    mpz_init(pub->n);
}


void paillier_private_init_v(paillier_private_key_v *priv) {
    mpz_init(priv->lambda);
    mpz_init(priv->mu);
    mpz_init(priv->p2);
    mpz_init(priv->q2);
    mpz_init(priv->p2invq2);
    mpz_init(priv->ninv);
    mpz_init(priv->n);
    mpz_init(priv->d);
}

//void decompose_into_digits(mpz_t number, vector<mpz_t> * digits, mpz_t base_sk, paillier_public_key *pub){
void decompose_into_digits(mpz_t number, mpz_t digits[], paillier_public_key *pub, mpz_t base_sk, int * sizeOfDigits){
    mpz_t remainder, quotient;
    mpz_init(remainder);
    mpz_init(quotient);
    int j=0;
    int qs[maxSizeOfDigits];

    // to be optimize
    //decompose number
    //gmp_printf("d is: %Zd\n", number);
    //cout << "the size of d is " << sizeof(number) << endl;
    gmp_printf("base_sk is: %Zd\n", base_sk);

    //cout << "d is: " << number << endl;
    mpz_fdiv_qr(quotient, remainder, number, base_sk);
    //cout << "decompose into: " << endl;
    while(mpz_cmp_ui(quotient, 0) != 0){
        qs[j] = mpz_get_ui(remainder);
        if(j<maxSizeOfDigits-1) j++;//furthudr change
        else{
            cout << "Too many digits for key! Try to increase Bsk." << endl;
            break;
        }
        //cout << "next step: " << quotient << ", " << remainder << endl;
        mpz_fdiv_qr(quotient, remainder, quotient, base_sk);
    }
    //cout << "next step: " << quotient << ", " << remainder << endl;
    qs[j] =mpz_get_ui(remainder);
    (*sizeOfDigits) = j+1;

    for(int i=j; i>=0; i--){
        mpz_init_set_ui(digits[i], qs[i]);
        //cout << "put " << digits[i] << " into: " << i << endl;
    }
/*
  for(int i=j+1; i<KeyBitLength; i++){
     mpz_init_set_ui(digits[i], 0);
     //cout << "put " << digits[i] << " into: " << i << endl;
  }*/
    cout << "totally digits: " << sizeOfDigits << endl;


}

//modified funcs
void setup_PKI_out(hss_pk_p *hpk, hss_ek_p *hek0, hss_ek_p *hek1, int * size, mpz_t base_sk, mpz_t base_msg){
    setup_PKI(hpk, hek0, hek1, size, base_sk, base_msg);
    //*size = sizeOfDigits;
}

void setup_PKI(hss_pk_p *hpk, hss_ek_p *hek0, hss_ek_p *hek1, int * sizeOfDigits, mpz_t base_sk, mpz_t base_msg){

    int result;
    paillier_public_key pub;
    paillier_private_key_v priv;

    //initialize keys
    cout << "System keygen: " << endl;
    paillier_public_init(&pub);
    paillier_private_init_v(&priv);
    //paillier_public_init(hpk->pk);
    //generate public key and private key of PAILLIER_DEBUG
    result = paillier_keygen_variant(&pub, &priv, PaillierKeySize/2);

    //--Subtractively secret share the digits of dF
    //cout << "-Subtractively secret share the digits of d." << endl;
    mpz_t tmp;
    //mpz_t base_sk, base_msg;

    mpz_init_set_ui(base_sk, BASE_SK);
    mpz_init(base_msg);
    mpz_init(tmp);

    mpz_t share0, share1;
    mpz_init(share0);
    mpz_init(share1);

    double pow2k, kappa=32;
    pow2k = pow(2, kappa);

    gmp_randstate_t state;
    gmp_randinit_mt(state);

    mpz_t remainder;
    mpz_init(remainder);

    //cout << "-Calculate b_msg." << endl;
    mpz_mul_ui(tmp, base_sk, pow2k);
    mpz_fdiv_qr(base_msg, remainder, pub.n, tmp);
    //if(mpz_cmp_ui(remainder,0) != 0){
    //  cout << "N/Bsk*2k is not integer." << endl;
    //}


    ///calculate base_msg from base_sk
    //mpz_mul_ui(base_msg, base_sk, pow2k);


    //cout << "-decompose begin." << endl;
    //cout << "where d is :" << priv.d << endl;
    mpz_t digits_of_d[maxSizeOfDigits];//furthudr change
    decompose_into_digits(priv.d, digits_of_d, &pub, base_sk, sizeOfDigits);
    //cout << "-decompose end." << endl;

    //vector<mpz_t> shares_of_digits_of_d_0;
    //vector<mpz_t> shares_of_digits_of_d_1;
    mpz_t shares_of_digits_of_d_0[maxSizeOfDigits];
    mpz_t shares_of_digits_of_d_1[maxSizeOfDigits];

    //cout << "prepare share of di." << endl;
    //cout << sizeOfDigits  << endl;
    int ist;

    /*
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    mpz_urandomm(randtmp, rstate, (hpk->pk)->n);
    mpz_set(mof11->sharex, randtmp);
    mpz_add_ui(randtmp, randtmp, 1);
    mpz_set(mof12->sharex, randtmp);
    */

    for(ist=0; ist<(*sizeOfDigits); ist++){
        //cout << i ;
        mpz_set_ui(share1, 123);
        //mpz_urandomm(share1, state, base_msg);
        //cout << i ;
        mpz_sub(share0 ,share1, digits_of_d[ist]);
        //shares_of_digits_of_d_0.push_back(share0);
        //shares_of_digits_of_d_1.push_back(share1);

        mpz_init_set(shares_of_digits_of_d_0[ist],share0);
        mpz_init_set(shares_of_digits_of_d_1[ist],share1);
    }
    //cout << digits_of_d[ist] <<  "digits of d 0 " << shares_of_digits_of_d_0[0] << "digits of d 1 " << shares_of_digits_of_d_1[0] << endl;

    for(; ist<KeyBitLength; ist++){
        //cout << i ;
        mpz_set_ui(share1, 321);
        //mpz_urandomm(share1, state, base_msg);
        //cout << i ;
        mpz_sub_ui(share0 ,share1, 0);
        //shares_of_digits_of_d_0.push_back(share0);
        //shares_of_digits_of_d_1.push_back(share1);

        mpz_init_set(shares_of_digits_of_d_0[ist],share0);
        mpz_init_set(shares_of_digits_of_d_1[ist],share1);
    }

    //mpz_urandomm(tmp, state, base_msg);
    ///...

    //--Sample a key for PRF
    //cout << "-Sample a key for PRF." << endl;
    mpz_t prf_key_num;
    mpz_init(prf_key_num);
    mpz_set_ui(prf_key_num, 202105);//seed, for further modification

    /*
    //the way to use seed as a prf_key
    gmp_randstate_t prf_state;
    gmp_randinit_mt(prf_state);
    gmp_randseed(prf_state, prf_key_num);*/

    //--pack ek0 ek1
    //prf_key_num
    //vector<mpz_t> shares_of_digits_of_d_0;
    //vector<mpz_t> shares_of_digits_of_d_1;
    mpz_init_set(hek0->kprf, prf_key_num); //changed from no init
    mpz_init_set(hek1->kprf, prf_key_num);

    for(ist=0; ist<KeyBitLength; ist++){
        mpz_init_set(hek0->ds[ist], shares_of_digits_of_d_0[ist]);
        mpz_init_set(hek1->ds[ist], shares_of_digits_of_d_1[ist]);
    }

    //Encrypt the digits of d
    //cout << "-Encrypt the digits of d." << endl;
    //int result;
    //vector<mpz_t> encryption_of_digits_of_d;
    mpz_t encryption_of_digits_of_d[maxSizeOfDigits];

    //for(int i; i<digits_of_d.size(); i++){
    for(ist=0; ist<(*sizeOfDigits); ist++){
        mpz_t tmped;
        mpz_init(tmped);
        result = paillier_encrypt(tmped, digits_of_d[ist], &pub);
        //encryption_of_digits_of_d.push_back(tmp);
        mpz_init_set(encryption_of_digits_of_d[ist],tmped);
    }

    for(; ist<KeyBitLength; ist++){
        mpz_t tmped, tmped2;
        mpz_init(tmped);
        mpz_init_set_ui(tmped2, 0);
        result = paillier_encrypt(tmped, tmped2, &pub);
        //encryption_of_digits_of_d.push_back(tmp);
        mpz_init_set(encryption_of_digits_of_d[ist],tmped);
    }


    //pack pk0 pk1
    //pk encryption_of_digits_of_d
    //cout << "-pack pk0 pk1." << endl;
    (hpk->pk)->len = PaillierKeySize/2;
    mpz_init_set((hpk->pk)->n, pub.n);
    //cout << "the n is:" << pub.n << endl;
    for(int i=0; i<KeyBitLength; i++){
        mpz_init_set(hpk->Ds[i], encryption_of_digits_of_d[i]);
        //cout << "Ds:" << hpk->Ds[i] << endl;
    }

}

void input_paillier(hss_input_p *Ix, hss_pk_p *hpk, mpz_t number, int * sizeOfDigits){
//(vector<mpz_t> shares_of_digits, mpz_t share, mpz_t number){


    //paillier_public_key pub = hpk->pk;

    int result;

    //gmp_randstate_t random_state;
    //gmp_randinit_mt(random_state);
    mpz_t u_r;
    mpz_init(u_r);
    mpz_set_ui(u_r, 16);
    //mpz_urandomm(u_r, random_state, pub->n);//furthur

    mpz_init(Ix->encx);

    //cout << "check point: " << (hpk->pk)->n << endl;
    result = paillier_encrypt(Ix->encx, number, hpk->pk);
    //cout << "mul: " << Ix->encx << endl;
    for(int i=0; i<(*sizeOfDigits); i++){
        mpz_t share;
        mpz_init(share);
        //result = paillier_homomorphic_multc(share, encryption_of_digits_of_d[i], number, &pub);
        //re-randomizing
        //result = paillier_homomorphic_multc_random(share, hpk->Ds[i], u_r, number, pub);
        result = paillier_homomorphic_multc(share, hpk->Ds[i], number, hpk->pk);
        //shares_of_digits.push(share);
        mpz_init_set(Ix->encdx[i], share);
        //cout << "mulc: " << Ix->encdx[i] << endl;
    }
    mpz_clear(u_r);
    //mpz_clear(share);
}

void input_paillier_value(hss_input_p *Ix, hss_pk_p *hpk, mpz_t number, int * sizeOfDigits){
//(vector<mpz_t> shares_of_digits, mpz_t share, mpz_t number){
    //paillier_public_key pub = hpk->pk;

    int result;

    //mpz_init(Ix->encx);

    //cout << "check point: " << (hpk->pk)->n << endl;
    result = paillier_encrypt(Ix->encx, number, hpk->pk);
    //cout << "mul: " << Ix->encx << endl;
    for(int i=0; i<(*sizeOfDigits); i++){
        //mpz_t share;
        //mpz_init(share);
        //result = paillier_homomorphic_multc(share, encryption_of_digits_of_d[i], number, &pub);
        //re-randomizing
        //result = paillier_homomorphic_multc_random(share, hpk->Ds[i], u_r, number, pub);
        result = paillier_homomorphic_multc(Ix->encdx[i], hpk->Ds[i], number, hpk->pk);
        //shares_of_digits.push(share);
        //mpz_set(Ix->encdx[i], share);
        //cout << "mulc: " << Ix->encdx[i] << endl;
    }
    //mpz_clear(u_r);
    //mpz_clear(share);
}

void input_paillier_clear(hss_input_p *Ix, int * sizeOfDigits){
    mpz_clear(Ix->encx);
    for(int i=0; i<(*sizeOfDigits); i++){
        mpz_clear(Ix->encdx[i]);
    }
}

void mult_shares_paillier(hss_men_p * Mx2, hss_men_p * Mx, hss_input_p *Ix, hss_ek_p * ek, paillier_public_key *pub, int serverId, int * sizeOfDigits, mpz_t base_sk){

    mpz_t tmp1,tmp2,tmp3,tmp4;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init_set_ui(tmp3, 8);
    mpz_init(tmp4);
    //mpz_init(result);

    //cout << "--reconstruct <yd>:" << endl;
    mpz_set_ui(tmp1, 0);
    for(int i=(*sizeOfDigits)-1; i>0; i--){
        //cout << Mx->sharedx[KeyBitLength-i-1] << " ";
        mpz_add(tmp1, Mx->sharedx[i], tmp1);
        mpz_mul(tmp1, tmp1, base_sk);
        //cout << "recover from: " << Mx->sharedx[sizeOfDigits-i] << "and have: " << tmp1 << endl;

    }

    mpz_add(tmp4, Mx->sharedx[0], tmp1);

    //share_divc(tmp4, tmp4, 256);
    //cout << "recover from: " << Mx->sharedx[0] << "and have: " << tmp4 << endl;
    //cout << tmp4 << endl;


    //cout << "--yield new result <z>." << endl;
    mpz_t n2;
    mpz_init(n2);
    mpz_mul(n2, pub->n, pub->n);
    //cout << tmp4 << endl;
    mpz_powm(tmp2, Ix->encx, tmp4, n2);
    //cout << tmp2 << endl;
    //cout << pub->n << endl;
    DDLog(tmp1, tmp2, pub);
    //cout << "Z: " << tmp1 << endl;


    //cout << "--generate random num." << endl;
    //gmp_randstate_t prf_state;
    //gmp_randinit_mt(prf_state);
    //mpz_add_ui(tmp3, ek->kprf, 1);
    //gmp_randseed(prf_state, tmp3);

    //cout << "--test." << endl;
    //cout << pub->n << endl;
    //mpz_urandomm(tmp2, prf_state, pub->n);
    mpz_set_ui(tmp2, 33333);

    mpz_add(tmp1, tmp1, tmp3);//suppose to be tmp2, but here cheat
    //mpz_init(Mx2->sharex);
    //cout << "--rerandomize." << endl;
    mpz_mod(Mx2->sharex, tmp1, pub->n);
    //cout << "z:" << Mx2->sharex << endl;
    //cout << "--yield new result <dz>." << endl;

    for(int i=0; i<(*sizeOfDigits); i++){
        //mpz_init(Mx2->sharedx[i]);
        //cout << "in: " << Ix->encdx[i] << endl;
        mpz_powm(tmp2, Ix->encdx[i], tmp4, n2);//furthur
        DDLog(tmp1, tmp2, pub);
        //cout << tmp1 << endl;
        mpz_add(tmp1, tmp1, tmp3);//suppose to be tmp2, but here cheat
        mpz_mod(Mx2->sharedx[i], tmp1, pub->n);  //only add here
        //cout << "zd:" << Mx2->sharedx[i] << endl;

    }


    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
    mpz_clear(tmp4);
    mpz_clear(n2);

}

void mult_shares_paillier_simpler(mpz_t result, hss_men_p * Mx, mpz_t Ixx, hss_ek_p * ek, paillier_public_key *pub, int serverId, int * sizeOfDigits, mpz_t base_sk){

    mpz_t tmp1,tmp2,tmp3,tmp4;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init_set_ui(tmp3, 8);
    mpz_init(tmp4);
    //mpz_init(result);

    //cout << "--reconstruct <yd>:" << endl;
    mpz_set_ui(tmp1, 0);
    for(int i=(*sizeOfDigits)-1; i>0; i--){
        //cout << Mx->sharedx[KeyBitLength-i-1] << " ";
        mpz_add(tmp1, Mx->sharedx[i], tmp1);
        mpz_mul(tmp1, tmp1, base_sk);
        //cout << "recover from: " << Mx->sharedx[sizeOfDigits-i] << "and have: " << tmp1 << endl;
    }

    mpz_add(tmp4, Mx->sharedx[0], tmp1);

    //cout << "recovered yd: " << tmp4 << endl;


    mpz_t n2;
    mpz_init(n2);
    mpz_mul(n2, pub->n, pub->n);
    //cout << tmp4 << endl;
    mpz_powm(tmp2, Ixx, tmp4, n2);
    //cout << tmp2 << endl;
    //cout << pub->n << endl;
    DDLog(tmp1, tmp2, pub);

    mpz_set_ui(tmp2, 33333);

    mpz_add(tmp1, tmp1, tmp3);//suppose to be tmp2, but here cheat
    //mpz_init(Mx2->sharex);
    //cout << "--rerandomize." << endl;
    mpz_mod(result, tmp1, pub->n);
    //cout << "z:" << Mx2->sharex << endl;
    //cout << "--yield new result <dz>." << endl;


    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
    mpz_clear(tmp4);
    mpz_clear(n2);

}

void mult_shares_paillier_prepare(hss_men_p * Mx, mpz_t Mxx, hss_ek_p * ek, paillier_public_key *pub, int serverId, int * sizeOfDigits, mpz_t base_sk) {
    mpz_t tmp1,tmp2,tmp3,tmp4;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init_set_ui(tmp3, 8);
    mpz_init(tmp4);
    //mpz_init(result);

    //cout << "--reconstruct <yd>:" << endl;
    mpz_set_ui(tmp1, 0);
    for(int i=(*sizeOfDigits)-1; i>0; i--){
        //cout << Mx->sharedx[KeyBitLength-i-1] << " ";
        mpz_add(tmp1, Mx->sharedx[i], tmp1);
        mpz_mul(tmp1, tmp1, base_sk);
        //cout << "recover from: " << Mx->sharedx[sizeOfDigits-i] << "and have: " << tmp1 << endl;
    }

    mpz_add(Mxx, Mx->sharedx[0], tmp1);

    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
    mpz_clear(tmp4);
}

void mult_shares_paillier_efficient(mpz_t result, mpz_t Mxx, mpz_t Ixx, hss_ek_p * ek, paillier_public_key *pub, int serverId, int * sizeOfDigits, mpz_t base_sk) {
    mpz_t tmp1,tmp2,tmp3,tmp4;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init_set_ui(tmp3, 8);
    mpz_init(tmp4);
    //mpz_init(result);


    mpz_t n2;
    mpz_init(n2);
    mpz_mul(n2, pub->n, pub->n);
    //cout << tmp4 << endl;
    mpz_powm(tmp2, Ixx, Mxx, n2);
    //cout << tmp2 << endl;
    //cout << pub->n << endl;
    DDLog(tmp1, tmp2, pub);

    mpz_set_ui(tmp2, 33333);

    mpz_add(tmp1, tmp1, tmp3);//suppose to be tmp2, but here cheat
    //mpz_init(Mx2->sharex);
    //cout << "--rerandomize." << endl;
    mpz_mod(result, tmp1, pub->n);
    //cout << "z:" << Mx2->sharex << endl;
    //cout << "--yield new result <dz>." << endl;


    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
    mpz_clear(tmp4);
    mpz_clear(n2);
}

void fromPaillierToAddShare(mpz_t shareP, mpz_t encP, hss_men_p * Mof1, hss_ek_p * ek, paillier_public_key *pub, int serverId, int * sizeOfDigits, mpz_t base_sk){

    mpz_t tmp1,tmp2,tmp3,tmp4;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init_set_ui(tmp3, 8);
    mpz_init(tmp4);
    //mpz_init(result);

    //cout << "--reconstruct <yd>:" << endl;
    mpz_set_ui(tmp1, 0);
    for(int i=(*sizeOfDigits)-1; i>0; i--){
        //cout << Mx->sharedx[KeyBitLength-i-1] << " ";
        mpz_add(tmp1, Mof1->sharedx[i], tmp1);
        mpz_mul(tmp1, tmp1, base_sk);
        //cout << "recover from: " << Mx->sharedx[sizeOfDigits-i] << "and have: " << tmp1 << endl;

    }

    mpz_add(tmp4, Mof1->sharedx[0], tmp1);

    //share_divc(tmp4, tmp4, 256);
    //cout << "recover from: " << Mx->sharedx[0] << "and have: " << tmp4 << endl;
    //cout << tmp4 << endl;


    //cout << "--yield new result <z>." << endl;
    mpz_t n2;
    mpz_init(n2);
    mpz_mul(n2, pub->n, pub->n);
    mpz_powm(tmp2, encP, tmp4, n2);
    //cout << tmp2 << endl;
    //cout << pub->n << endl;
    DDLog(tmp1, tmp2, pub);
    //cout << "Z: " << tmp1 << endl;


    //cout << "--generate random num." << endl;
    //gmp_randstate_t prf_state;
    //gmp_randinit_mt(prf_state);
    //mpz_add_ui(tmp3, ek->kprf, 1);
    //gmp_randseed(prf_state, tmp3);

    //cout << "--test." << endl;
    //cout << pub->n << endl;
    //mpz_urandomm(tmp2, prf_state, pub->n);
    mpz_set_ui(tmp2, 33333);

    mpz_add(tmp1, tmp1, tmp3);//suppose to be tmp2, but here cheat
    //mpz_init(Mx2->sharex);
    //cout << "--rerandomize." << endl;
    mpz_mod(shareP, tmp1, pub->n);
    //cout << "z:" << Mx2->sharex << endl;
    //cout << "--yield new result <dz>." << endl;


    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
    mpz_clear(tmp4);
    mpz_clear(n2);

}


void convert_shares_paillier(hss_men_p * Mx, hss_men_p * Mof1, hss_input_p *Ix, hss_ek_p * ek, paillier_public_key *pub, int serverId, int * sizeOfDigits, mpz_t base_sk){
    //mpz_t share0, share1;
    //mpz_init(share0);
    //mpz_init(share1);

/*
  //the way to use seed as a prf_key
  gmp_randstate_t prf_state;
  gmp_randinit_mt(prf_state);
  gmp_randseed(prf_state, prf_key_num);
  mpz_urandomm(share0, prf_state, pub->n);
  mpz_add_ui(share1, share0, 1);
*/

    //mpz_t tmp;
    //mpz_init(tmp);
    //mpz_add_ui(tmp, ek->kprf, serverId);//cheat

    mult_shares_paillier(Mx, Mof1, Ix, ek, pub, serverId, sizeOfDigits, base_sk);


    //mult_shares_paillier();

    //mpz_clear(share0);
    //mpz_clear(share1);

}


void add_shares_paillier(hss_men_p * Mz, hss_men_p * Mx, hss_men_p * My, hss_ek_p * ek, paillier_public_key *pub, int serverId){
    //mpz_init(Mz->sharex);
    mpz_add(Mz->sharex, Mx->sharex, My->sharex);
    mpz_mod(Mz->sharex, Mz->sharex, pub->n);
    for(int i=0; i<KeyBitLength; i++){
        //mpz_init(Mz->sharedx[i]);
        mpz_add(Mz->sharedx[i], Mx->sharedx[i], My->sharedx[i]);
        mpz_mod(Mz->sharedx[i], Mz->sharedx[i], pub->n);
    }
}

void sub_shares_paillier(hss_men_p * Mz, hss_men_p * Mx, hss_men_p * My, hss_ek_p * ek, paillier_public_key *pub, int serverId, int * sizeOfDigits){
    //mpz_init(Mz->sharex);
    mpz_sub(Mz->sharex, Mx->sharex, My->sharex);
    mpz_mod(Mz->sharex, Mz->sharex, pub->n);
    //cout << Mz->sharex << ";" << endl;
    for(int i=0; i<(*sizeOfDigits); i++){
        //mpz_init(Mz->sharedx[i]);
        mpz_sub(Mz->sharedx[i], Mx->sharedx[i], My->sharedx[i]);
        mpz_mod(Mz->sharedx[i], Mz->sharedx[i], pub->n);
    }
}

void add_inputs_paillier(hss_input_p *Iz, hss_input_p *Ix, hss_input_p *Iy, hss_ek_p * ek, paillier_public_key *pub, int serverId){
    mpz_init(Iz->encx);
    paillier_homomorphic_add(Iz->encx, Ix->encx, Iy->encx, pub);
    for(int i=0; i<KeyBitLength; i++){
        mpz_init(Iz->encdx[i]);
        paillier_homomorphic_add(Iz->encdx[i], Ix->encdx[i], Iy->encdx[i], pub);
        //mpz_mod(Mx2->sharedx[i], pub->n);
    }
}

void output_paillier(mpz_t out, hss_men_p * Mz, hss_ek_p * ek, paillier_public_key *pub, int serverId){
    mpz_init(out);
    mpz_set(out, Mz->sharex);
}

void men_share_init(hss_men_p* Mz){
    mpz_init_set_ui(Mz->sharex, 0);
    for (int i = 0; i < KeyBitLength; i++) {
        mpz_init_set_ui(Mz->sharedx[i], 0);
    }
}

void men_share_zeros(hss_men_p* Mz){
    mpz_set_ui(Mz->sharex, 0);
    for (int i = 0; i < KeyBitLength; i++) {
        mpz_set_ui(Mz->sharedx[i], 0);
    }
}

void men_share_clear(hss_men_p* Mz){
    mpz_clear(Mz->sharex);
    for (int i = 0; i < KeyBitLength; i++) {
        mpz_clear(Mz->sharedx[i]);
    }
}


void men_share_multc(hss_men_p* Mz, hss_men_p* Mx, int constNum, paillier_public_key *pub, int * sizeOfDigits){
    mpz_mul_ui(Mz->sharex, Mx->sharex, constNum);
    mpz_mod(Mz->sharex, Mz->sharex, pub->n);
    for (int i = 0; i < (*sizeOfDigits); i++) {
        mpz_mul_ui(Mz->sharedx[i], Mx->sharedx[i], constNum);
        mpz_mod(Mz->sharedx[i], Mz->sharedx[i], pub->n);
    }
}

void share_divc(mpz_t sharex, mpz_t share, int constNum){
    mpz_t quotient, remainder, tmp;
    mpz_init(quotient);
    mpz_init(remainder);
    mpz_init(tmp);

    mpz_fdiv_qr_ui(quotient,remainder,share,constNum);
    //cout << "share: " << share << endl;
    //cout << "q: " << quotient << endl;
    //cout << "r: " << remainder << endl;


    if(mpz_cmp_ui(share,0) >0){
        //mpz_fdiv_q_ui(tmp,base_sk,2);
        mpz_set_ui(tmp, int (256/2));
        //cout << "base_sk/2: " << tmp << endl;
        if(mpz_cmp(remainder,tmp) >0)mpz_add_ui(quotient,quotient,1);
    }else{
        mpz_neg(remainder,remainder);
        mpz_fdiv_q_ui(tmp,share,2);
        if(mpz_cmp(remainder,tmp) >0)mpz_sub_ui(quotient,quotient,1);
    }
    mpz_set(sharex, quotient);

    mpz_clear(quotient);
    mpz_clear(remainder);
    mpz_clear(tmp);
}

void men_share_divc(hss_men_p* Mz, int constNum){
    mpz_t tmp;
    mpz_init(tmp);
    //cout << "in zt: " << Mz->sharex << endl;
    share_divc(Mz->sharex, Mz->sharex, constNum);
    //cout << "out zt: "<< Mz->sharex << endl;
    for (int i = 0; i < KeyBitLength; i++) {
    //cout << "in zdit: " << Mz->sharedx[i] << endl;
    share_divc(Mz->sharedx[i], Mz->sharedx[i], constNum);
    //cout << "out zdit: " << Mz->sharedx[i] << endl;
    }

    mpz_clear(tmp);
}
/*
void men_share_divc(hss_men_p* Mz, int constNum){//origin
    mpz_t tmp;
    mpz_init(tmp);
    cout << "in zt: " << Mz->sharex << endl;
    share_divc(Mz->sharex, Mz->sharex, constNum);
    cout << "out zt: "<< Mz->sharex << endl;
    for (int i = 0; i < KeyBitLength; i++) {
      cout << "in zdit: " << Mz->sharedx[i] << endl;
      share_divc(Mz->sharedx[i], Mz->sharedx[i], constNum);
      cout << "out zdit: " << Mz->sharedx[i] << endl;
    }

    mpz_clear(tmp);
}*/
/*
void men_share_multc(hss_men_p* Mz, hss_men_p* Mx, int constNum, paillier_public_key *pub){
    mpz_t tmp;
    mpz_init_set_ui(tmp, constNum);
    //cout << "in zt: " << Mz->sharex << endl;
    paillier_homomorphic_multc(Mz->sharex, Mx->sharex, tmp, pub);
    mpz_mod(Mz->sharex, Mz->sharex, pub->n);
    //share_divc(Mz->sharex, Mz->sharex, constNum);
    //cout << "out zt: "<< Mz->sharex << endl;
    for (int i = 0; i < KeyBitLength; i++) {
      //cout << "in zdit: " << Mz->sharedx[i] << endl;
      paillier_homomorphic_multc(Mz->sharedx[i], Mx->sharedx[i], tmp, pub);
      mpz_mod(Mz->sharedx[i], Mz->sharedx[i], pub->n);
      //cout << "out zdit: " << Mz->sharedx[i] << endl;
    }

    mpz_clear(tmp);


}*/

void men_share_divc_cheat(hss_men_p* Mz, int constNum, paillier_private_key_v *priv){
    mpz_t tmp;
    mpz_init(tmp);

    share_divc(Mz->sharex, Mz->sharex, constNum);
    for (int i = 0; i < KeyBitLength; i++) {
        mpz_mul(Mz->sharedx[i], Mz->sharedx[i], priv->d);
    }
    mpz_clear(tmp);
}

void input_copy_pailliar(hss_input_p* Mz, hss_input_p* Mx, int * sizeOfDigits){
    mpz_set(Mz->encx, Mx->encx);
    //cout << "still ok" << endl;
    //mpz_mod(Mz->sharex, Mz->sharex, pub->n);
    for (int i = 0; i < (*sizeOfDigits); i++) {
        mpz_set(Mz->encdx[i], Mx->encdx[i]);
        //mpz_mod(Mz->sharedx[i], Mz->sharedx[i], pub->n);
    }
}

void men_copy_pailliar(hss_men_p* Mz, hss_men_p* Mx, int * sizeOfDigits){
    mpz_set(Mz->sharex, Mx->sharex);
    //mpz_mod(Mz->sharex, Mz->sharex, pub->n);
    for (int i = 0; i < (*sizeOfDigits); i++) {
        mpz_set(Mz->sharedx[i], Mx->sharedx[i]);
        //mpz_mod(Mz->sharedx[i], Mz->sharedx[i], pub->n);
    }
}


// end modification

int fromINTtoMEN(hss_men_p * mnum11, hss_men_p * mnum12, hss_men_p * mof11, hss_men_p * mof12, hss_pk_p* hpk, hss_ek_p * ek0, hss_ek_p * ek1, mpz_t num, MyTime* iktime,
                 int * sizeOfDigits, mpz_t base_sk){
    hss_input_p inum1;
    mpz_t randtmp;
    mpz_init(randtmp);

    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    mpz_urandomm(randtmp, rstate, (hpk->pk)->n);
    mpz_set(mof11->sharex, randtmp);
    mpz_add_ui(randtmp, randtmp, 1);
    mpz_set(mof12->sharex, randtmp);

    //prepare inputs
    //mpz_set(hpk.pk->n, nst);
    //cout << "check n:" << (hpk->pk)->n << endl;
    input_paillier(&inum1, hpk, num, sizeOfDigits);

    //cout << "check input values1: " << inum1.encx << endl;

    //convert into memory shares
    //cout << "convert into memory shares." << endl;

    iktime->Time_begin(12);

    men_share_init(mnum11);
    men_share_init(mnum12);
    convert_shares_paillier(mnum11, mof11, &inum1, ek0, hpk->pk, 0, sizeOfDigits, base_sk);
    convert_shares_paillier(mnum12, mof12, &inum1, ek1, hpk->pk, 1, sizeOfDigits, base_sk);

    iktime->Time_end(12);

    mpz_clear(randtmp);
    return 1;
}

int fromINTtoMEN_value(hss_men_p * mnum11, hss_men_p * mnum12, hss_men_p * mof11, hss_men_p * mof12, hss_pk_p* hpk, hss_ek_p * ek0, hss_ek_p * ek1, mpz_t num, MyTime* iktime,
                       int * sizeOfDigits, mpz_t base_sk){
    hss_input_p inum1;
    mpz_t randtmp1;
    mpz_init(randtmp1);
/*
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    mpz_urandomm(randtmp1, rstate, (hpk->pk)->n);
    mpz_set(mof11->sharex, randtmp1);
    mpz_add_ui(randtmp1, randtmp1, 1);
    mpz_set(mof12->sharex, randtmp1);
*/
    //prepare inputs
    //mpz_set(hpk.pk->n, nst);
    //cout << "check n:" << (hpk->pk)->n << endl;
    input_paillier(&inum1, hpk, num, sizeOfDigits);

    //cout << "check input values1: " << inum1.encx << endl;

    //convert into memory shares
    //cout << "convert into memory shares." << endl;

    //men_share_init(mnum11);
    //men_share_init(mnum12);
    iktime->Time_begin(12);
    convert_shares_paillier(mnum11, mof11, &inum1, ek0, hpk->pk, 0, sizeOfDigits, base_sk);
    convert_shares_paillier(mnum12, mof12, &inum1, ek1, hpk->pk, 1, sizeOfDigits, base_sk);
    iktime->Time_end(12);
    //cout << "Convert Over" << endl;

    //input_paillier_clear(&inum1);

    //mpz_clear(randtmp1);
    return 1;
}


int fromINTtoMEN_single(hss_men_p * mnum11, hss_men_p * mnum12, hss_men_p * mof11, hss_men_p * mof12, hss_pk_p* hpk, hss_ek_p * ek0, hss_ek_p * ek1, mpz_t num, MyTime* iktime,
                        int server_id, int * sizeOfDigits, mpz_t base_sk){
    hss_input_p inum1;
    mpz_t randtmp;
    mpz_init(randtmp);



    //prepare inputs
    //mpz_set(hpk.pk->n, nst);
    //cout << "check n:" << (hpk->pk)->n << endl;
    input_paillier(&inum1, hpk, num, sizeOfDigits);

    //cout << "check input values1: " << inum1.encx << endl;

    //convert into memory shares
    //cout << "convert into memory shares." << endl;

    iktime->Time_begin(12);

    men_share_init(mnum11);
    men_share_init(mnum12);
    convert_shares_paillier(mnum11, mof11, &inum1, ek0, hpk->pk, 0, sizeOfDigits, base_sk);
    convert_shares_paillier(mnum12, mof12, &inum1, ek1, hpk->pk, 1, sizeOfDigits, base_sk);

    iktime->Time_end(12);

    mpz_clear(randtmp);
    return 1;
}

int fromINTtoMEN_value_single(hss_men_p * mnum11, hss_men_p * mnum12, hss_men_p * mof11, hss_men_p * mof12, hss_pk_p* hpk, hss_ek_p * ek0, hss_ek_p * ek1, mpz_t num, MyTime* iktime,
                              int server_id, int * sizeOfDigits, mpz_t base_sk){
    hss_input_p inum1;
    mpz_t randtmp1;
    mpz_init(randtmp1);
/*
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    mpz_urandomm(randtmp1, rstate, (hpk->pk)->n);
    mpz_set(mof11->sharex, randtmp1);
    mpz_add_ui(randtmp1, randtmp1, 1);
    mpz_set(mof12->sharex, randtmp1);
*/
    //prepare inputs
    //mpz_set(hpk.pk->n, nst);
    //cout << "check n:" << (hpk->pk)->n << endl;
    input_paillier(&inum1, hpk, num, sizeOfDigits);

    //cout << "check input values1: " << inum1.encx << endl;

    //convert into memory shares
    //cout << "convert into memory shares." << endl;

    //men_share_init(mnum11);
    //men_share_init(mnum12);
    iktime->Time_begin(12);
    convert_shares_paillier(mnum11, mof11, &inum1, ek0, hpk->pk, 0, sizeOfDigits, base_sk);
    convert_shares_paillier(mnum12, mof12, &inum1, ek1, hpk->pk, 1, sizeOfDigits, base_sk);
    iktime->Time_end(12);
    //cout << "Convert Over" << endl;

    //input_paillier_clear(&inum1);

    //mpz_clear(randtmp1);
    return 1;
}

void input_share_multc(hss_input_p * Iz, hss_input_p * Ix, int constNum, paillier_public_key *pub, int * sizeOfDigits){
    mpz_t tmp;
    mpz_init_set_ui(tmp, constNum);

    //cout << "in z: " << Iz->encx << endl;
    paillier_homomorphic_multc(Iz->encx, Ix->encx, tmp, pub);

    //cout << "length: " << sizeOfDigits << endl;

    for (int i = 0; i < (*sizeOfDigits); i++) {
        //cout << "  " << i << endl;
        //cout << "in zdit: " << Iz->encdx[i] << endl;
        paillier_homomorphic_multc(Iz->encdx[i], Ix->encdx[i], tmp, pub);
        mpz_mod(Iz->encdx[i], Iz->encdx[i], pub->n);
        //cout << "out zdit: " << Iz->encdx[i] << endl;
    }
    //cout << "over " << endl;
    mpz_clear(tmp);
}

int server_output(mpz_t data[], int server_id, char * output_file_name){
    unsigned int result_flag, sum_out = 0;

    FILE *fp_out = NULL;
    fp_out = fopen(output_file_name, "w+");

    for(int i=0; i<10; i++) {

        result_flag = mpz_out_raw(fp_out, data[i]);
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

int server_input(mpz_t data[], int server_id, char * output_file_name){
    unsigned int result_flag, sum_in = 0;
    FILE *fp_in = NULL;

    fp_in = fopen("raw_mpz.dat", "r");

    int result_num;
    for(int i=0; i<10; i++) {
        //mpzDataIn >> bits2;
        mpz_init(data[i]);
        //mpzDataIn.read(bits2, sizeof(mpz_t));
        //result_num = mpz_set_str(colors[i], bits2, 2);
        //if(result_num == -1){
        //    cerr << "Invalid number!" << endl;
        //    return -1;
        //}
        result_flag = mpz_inp_raw(data[i], fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }
        cout << data[i] << endl;
    }

    cout << "Result: Read " << sum_in << " bytes." << endl;

    fclose(fp_in);

    return 0;
}

int server_result_output(node_paillier (*input_qX1)[train_number+100], int server_id, char * output_file_name){
    unsigned int result_flag, sum_out = 0;

    FILE *fp_out = NULL;
    fp_out = fopen(output_file_name, "w+");

    for(int i=0; i<train_number; i++) {

        result_flag = mpz_out_raw(fp_out, (*input_qX1)[i].dis);
        if(result_flag == 0){
            cerr << "Unable to write!" << endl;
            return -1;
        }else{
            sum_out += result_flag;
        }

        result_flag = mpz_out_raw(fp_out, (*input_qX1)[i].labels);
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

int server_result_input(node_paillier (*input_qX1)[train_number+100], int server_id, char * output_file_name){
    unsigned int result_flag, sum_in = 0;
    FILE *fp_in = NULL;

    fp_in = fopen(output_file_name, "r");

    int result_num;
    for(int i=0; i<train_number; i++) {
        //mpzDataIn >> bits2;
        //mpz_init(data[i]);
        //mpzDataIn.read(bits2, sizeof(mpz_t));
        //result_num = mpz_set_str(colors[i], bits2, 2);
        //if(result_num == -1){
        //    cerr << "Invalid number!" << endl;
        //    return -1;
        //}
        result_flag = mpz_inp_raw((*input_qX1)[i].dis, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }

        result_flag = mpz_inp_raw((*input_qX1)[i].labels, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }

        //cout << data[i] << endl;
    }

    cout << "Result: Read " << sum_in << " bytes." << endl;

    fclose(fp_in);

    return 0;
}



int client_write_input(std::vector<std::vector<hss_input_p> >* Xi1, mpz_t Yw1[], int server_id,
                       char * output_file1_name, int * sizeOfDigits) {
    unsigned int result_flag, sum_out = 0;

    FILE * fp_out = NULL;
    fp_out = fopen(output_file1_name, "w+");


    for (int i = 0; i < train_number; i++)  //for each train sample
    {
        /*cout << "Sample " << i << "," << endl;
        vector<hss_input_p> temp_vec1;
        vector<hss_input_p> temp_vec2;*/

        for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi

            result_flag = mpz_out_raw(fp_out, (*Xi1)[i][t].encx);
            if(result_flag == 0){
                cerr << "Unable to write!" << endl;
                return -1;
            }else{
                sum_out += result_flag;
            }

            for(int j=0; j<(*sizeOfDigits); j++){
                result_flag = mpz_out_raw(fp_out, (*Xi1)[i][t].encdx[j]);
                if(result_flag == 0){
                    cerr << "Unable to write!" << endl;
                    return -1;
                }else{
                sum_out += result_flag;
                }
            }

        }

        result_flag = mpz_out_raw(fp_out, (Yw1)[i]);
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

int server_read_input_DO(std::vector<std::vector<hss_input_p> >* Xi1, mpz_t Yw1[], int server_id,
                      char * output_file1_name, int * sizeOfDigits){
    unsigned int result_flag, sum_in = 0;

    FILE *fp_in = NULL;
    fp_in = fopen(output_file1_name, "r");


    mpz_t tmpNum;
    mpz_init(tmpNum);

    cout << "Check digits size:" << *sizeOfDigits << endl;

    for (int i = 0; i < train_number; i++)  //for each train sample
    {

        //cout << i << endl;
        for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi

            //cout << t << endl;

            result_flag = mpz_inp_raw(tmpNum, fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }

            mpz_set((*Xi1)[i][t].encx, tmpNum);


            //cout << *sizeOfDigits << endl;
            for(int j=0; j<(*sizeOfDigits); j++){
                result_flag = mpz_inp_raw(tmpNum, fp_in);
                if(result_flag == 0){
                    cerr << "Unable to read!" << endl;
                    cerr << j << endl;
                    return -1;
                }else{
                    sum_in += result_flag;
                }
                mpz_set((*Xi1)[i][t].encdx[j], tmpNum);
            }

        }

        result_flag = mpz_inp_raw(tmpNum, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }
        mpz_set(Yw1[i], tmpNum);


    }

    cout << "Result: Read " << sum_in << " bytes." << endl;
    fclose(fp_in);

    return 0;
}

int server_read_input_DO_array(std::array<std::array<hss_input_p, attribute_number> , train_number >* Xi1, std::array<mpz_t, train_number> * Yw1,
                         char * output_file1_name, int * sizeOfDigits){
    unsigned int result_flag, sum_in = 0;

    FILE *fp_in = NULL;
    fp_in = fopen(output_file1_name, "r");


    mpz_t tmpNum;
    mpz_init(tmpNum);

    cout << "Check digits size:" << *sizeOfDigits << endl;

    for (int i = 0; i < train_number; i++)  //for each train sample
    {

        //cout << i << endl;
        for (int t = 0; t < attribute_number-1; t++) { //read a train samples t to Xi

            //cout << t << endl;

            result_flag = mpz_inp_raw(tmpNum, fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }

            mpz_init_set((*Xi1)[i][t].encx, tmpNum);


            //cout << *sizeOfDigits << endl;
            for(int j=0; j<(*sizeOfDigits); j++){
                result_flag = mpz_inp_raw(tmpNum, fp_in);
                if(result_flag == 0){
                    cerr << "Unable to read!" << endl;
                    cerr << j << endl;
                    return -1;
                }else{
                    sum_in += result_flag;
                }
                mpz_init_set((*Xi1)[i][t].encdx[j], tmpNum);
            }

        }

        result_flag = mpz_inp_raw(tmpNum, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }
        mpz_init_set((*Yw1)[i], tmpNum);
    }

    cout << "Result: Read " << sum_in << " bytes." << endl;
    fclose(fp_in);

    return 0;
}

int server_read_input_client(std::vector<std::vector<hss_input_p> >* Xi1, int server_id,
    char * output_file1_name, int * sizeOfDigits){
    unsigned int result_flag, sum_in = 0;

    FILE *fp_in = NULL;
    fp_in = fopen(output_file1_name, "r");


    mpz_t tmpNum;
    mpz_init(tmpNum);

    //cout << keySize << "," << sizeOfDigits << endl;


    for (int i = 0; i < test_number; i++)  //for each train sample
    {

        //cout << i << endl;
        for (int t = 0; t < attribute_number; t++) { //read a train samples t to Xi



            result_flag = mpz_inp_raw(tmpNum, fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                cerr << "i:" << i << endl;
                cerr << "t:" << t << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }

            mpz_set((*Xi1)[i][t].encx, tmpNum);

            for(int j=0; j<(*sizeOfDigits); j++){
                result_flag = mpz_inp_raw(tmpNum, fp_in);
                if(result_flag == 0){
                    cerr << "Unable to read!" << endl;
                    cerr << j << endl;
                    return -1;
                }else{
                    sum_in += result_flag;
                }
                mpz_set((*Xi1)[i][t].encdx[j], tmpNum);
            }

        }

    }

    cout << "Result: Read " << sum_in << " bytes." << endl;
    fclose(fp_in);

    return 0;
}

int server_read_input(std::vector<std::vector<hss_input_p> >* Xi1, mpz_t Yw1[], int server_id,
        char * output_file1_name, int * sizeOfDigits){
    unsigned int result_flag, sum_in = 0;

    FILE *fp_in = NULL;
    fp_in = fopen(output_file1_name, "r");




    mpz_t tmpNum;
    mpz_init(tmpNum);


    for (int i = 0; i < train_number; i++)  //for each train sample
    {

        //cout << i << endl;
        for (int t = 0; t < attribute_number-1; t++) { //read a train samples t to Xi



            result_flag = mpz_inp_raw(tmpNum, fp_in);
            if(result_flag == 0){
                cerr << "Unable to read!" << endl;
                return -1;
            }else{
                sum_in += result_flag;
            }

            mpz_set((*Xi1)[i][t].encx, tmpNum);


            for(int j=0; j<(*sizeOfDigits); j++){
                result_flag = mpz_inp_raw(tmpNum, fp_in);
                if(result_flag == 0){
                    cerr << "Unable to read!" << endl;
                    return -1;
                }else{
                    sum_in += result_flag;
                }
                mpz_set((*Xi1)[i][t].encdx[j], tmpNum);
            }

        }

        result_flag = mpz_inp_raw(tmpNum, fp_in);
        if(result_flag == 0){
            cerr << "Unable to read!" << endl;
            return -1;
        }else{
            sum_in += result_flag;
        }
        mpz_set(Yw1[i], tmpNum);


    }

    cout << "Result: Read " << sum_in << " bytes." << endl;
    fclose(fp_in);

    return 0;
}

int subtractive_share(mpz_t x1, mpz_t x2, int x, mpz_t n){
    //cout << "" << endl;
    mpz_t tmp_ran;
    mpz_init(tmp_ran);
    mpz_set_ui(tmp_ran, x);
    
    //gmp_randstate_t my_state;
    //gmp_randinit_mt(my_state);

    //mpz_urandomm(x1, my_state, n);

    mpz_set_ui(x1, 33221);
    mpz_add(x2, x1, tmp_ran);

    //mpz_clear(tmp_ran);
    return 0;
}

int prand_perm(mpz_t result[], mpz_t array[], mpz_t key, int len){
    long int a, b, init_x,max_len;
    long int prp_index[len];

    //tmp set...
    a = 5;
    b = 3;
    init_x = 0;
    prp_index[0] = (a * init_x + b) % len;
    for(int i=1; i< len; i++){
        prp_index[i] = (a * prp_index[i-1] + b) % (len);
        //mpz_set(result[a * prp_index[i] + b], array[prp_index[i]]);
    }
    prp_index[0] = (a * prp_index[len-1] + b) % len;
    //cout << "Show prp generation..." << endl;
    for(int i=0; i< len; i++){
        //prp_index[i] = a * prp_index[i] + b;
        mpz_set(result[prp_index[i]], array[i]);
        //cout << prp_index[i] << ", ";
    }
    //cout << endl;
    //mpz_set(result[prp_index[0]], result[prp_index[len-1]]);

    return 0;
}

int PRShuffle(mpz_t data[], mpz_t label[], int len, gmp_randstate_t state) {
    mpz_t max_rand;
    mpz_t tmpresult;
    mpz_t tmpresult2;

    mpz_init(max_rand);
    mpz_init(tmpresult);
    mpz_init(tmpresult2);


    for (int i = 0, j; i < len; i++) {

        //j = (i == 0) ? 0 : (static_cast<int>((std::rand()) % (i + 1));
        if(i == 0){
            j = 0;
        }else{
            mpz_set_ui(max_rand, i+1);
            mpz_urandomm(tmpresult, state, max_rand);
            j = mpz_get_ui(tmpresult);
        }

        if (j != i){
            mpz_set(tmpresult2, data[i]);
            mpz_set(data[i], data[j]);
            mpz_set(data[j], tmpresult2);

            mpz_set(tmpresult2, label[i]);
            mpz_set(label[i], label[j]);
            mpz_set(label[j], tmpresult2);
        }
        //ind[i] = ind[j];
        //ind[j] = i;
    }

    mpz_clear(max_rand);
    mpz_clear(tmpresult);
    mpz_clear(tmpresult2);

    return 0;
}


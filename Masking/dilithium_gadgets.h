#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "random.h"
#include "convtable.h"
#include "convba_2014.h"
#include "utils.h"



#ifndef DILITHIUM_MODE
#define DILITHIUM_MODE 2
#endif


#define DIL_Q 8380417
#define DIL_Q_PRIME 134086672
#define ALPHA 4
//We assume alpha = 4 which means we can go up to 15 shares

#if DILITHIUM_MODE == 2
#define MU 18
#define K_APPROX 41
#define A_APPROX 262401
#define K_EXACT 45
#define A_EXACT 4198404
#define DELTA 44
#define D_BETA 78
#define D_GAMMA1 (1 << 17)
#define D_GAMMA2 95232
#define RHO 8


#elif DILITHIUM_MODE == 3
#define MU 20
#define K_APPROX 43
#define A_APPROX 1049601
#define K_EXACT 47
#define A_EXACT 16793615
#define DELTA 16
#define D_BETA 196
#define D_GAMMA1 (1 << 19)
#define D_GAMMA2 261888
#define RHO 9

#elif DILITHIUM_MODE == 5
#define MU 20
#define K_APPROX 43
#define A_APPROX 1049601
#define K_EXACT 47
#define A_EXACT 16793615
#define DELTA 16
#define D_BETA 120
#define D_GAMMA1 (1 << 19)
#define D_GAMMA2 261888
#define RHO 8
#endif


#define PORTABLE

#ifndef PORTABLE
#define uint128_t __uint128_t
#define int128_t __int128_t
#endif


void securediv(uint32_t* x, uint32_t* y, int log_alpha, int inv_alpha, int q, int n);

void generic_shift(uint32_t* x, uint32_t* y, int k, int q, int n);


void reject_sampling(uint32_t* x, uint32_t* z, int mode, int q, int n);
int single_coeff_reject(uint32_t* x, int mode, int q, int n);
int old_rejection_sampling(uint32_t* x, int mode, int n);
int old_rejection_sampling_shared_output(uint32_t* x, uint32_t* y, int mode, int n);
int full_zero_test(uint32_t* coeffs, uint32_t size, int n);


void gen_y(uint32_t* y, int n);
void gen_y_fast(uint32_t* y, int n);
void gen_y_det(uint64_t* x, uint32_t* y, int n);



void approximate_modulus_switching(uint64_t* x, uint32_t* y, int n);
void exact_modulus_switching(uint64_t* x, uint32_t* y, int n);
void MSB_estimation(uint32_t* x, uint32_t* y, int n);
void LMSwitch(uint32_t* x, uint32_t* y, int q, int rho, int n);



void decompose1(uint32_t* x, uint32_t* high, uint32_t* low, int n);
void decompose2(uint32_t* x, uint32_t* high, uint32_t* low, int n);
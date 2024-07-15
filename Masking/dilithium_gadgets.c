#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "random.h"
#include "convtable.h"
#include "convba_2014.h"
#include "utils.h"
#include "dilithium_gadgets.h"





void generic_1bit_shift(uint32_t* x, uint32_t* y, int q, int n){
  /* Shiftmod Shift of 1 bit from mod q to mod q/2 for any even q*/
  uint32_t b[n], a[n], z[n];

  for(int i=0; i < n; ++i) b[i] = x[i]&1;
  bool2ArithSPOGmodq(b, a, q, n);


  for(int i=0; i < n; ++i) z[i] = (x[i] + q - a[i])%q;

  for(int i=0; i < n-1; ++i){
    z[n-1] = (z[n-1]   + (z[i]&1))%q;
    z[i]   = (z[i] + q - (z[i]&1))%q;
  }
  for(int i=0; i < n; ++i) y[i] = z[i]>>1;
}

void generic_1bit_shift64(uint64_t* x, uint64_t* y, uint64_t q, int n){
  /* Shift of 1 bit from mod q to mod q/2 for any even q*/
  uint64_t b[n], a[n], z[n];

  for(int i=0; i < n; ++i) b[i] = x[i]&1;
  bool2ArithSPOGmodq64(b, a, q, n);


  for(int i=0; i < n; ++i) z[i] = (x[i] + q - a[i])%q;

  for(int i=0; i < n-1; ++i){
    z[n-1] = (z[n-1]   + (z[i]&1))%q;
    z[i]   = (z[i] + q - (z[i]&1))%q;
  }
  for(int i=0; i < n; ++i) y[i] = z[i]>>1;
}


void generic_shift(uint32_t* x, uint32_t* y, int k, int q, int n){
  /* Shift of k bits from mod 2^k * q to mod q for any q*/
  for(int i=0; i < (k>>1); ++i){
    generic_1bit_shift(x, y, (1<<(k-2*i))*q, n);
    generic_1bit_shift(y, x, (1<<(k-(2*i+1)))*q, n);
  }
  if (k&1){
    generic_1bit_shift(x, y, 2*q, n);
  } else {
    for(int i=0; i < n; ++i) y[i] = x[i];
  }

}

void generic_shift64(uint64_t* x, uint64_t* y, int k, uint64_t q, int n){
  /* Shift of k bits from mod 2^k * q to mod q for any q*/
  for(int i=0; i < (k>>1); ++i){
    generic_1bit_shift64(x, y, (1<<(k-2*i))*q, n);
    generic_1bit_shift64(y, x, (1<<(k-(2*i+1)))*q, n);
  }
  if (k&1){
    generic_1bit_shift64(x, y, 2*q, n);
  } else {
    for(int i=0; i < n; ++i) y[i] = x[i];
  }

}


int old_rejection_sampling(uint32_t* x, int mode, int n){
  /* Assumes input masked mod+ q, returns 1 if  value is valid
     Method converting A2B and looking at sign bit, used in [ABC+22]
  */


  const int BITSIZE = 24;
  const int beta = D_BETA;
  uint32_t a = 0, res = 0;
  if (mode == 0)  a = D_GAMMA1 - beta;
  else if (mode == 1) a = ((DIL_Q-1)/(2*DELTA)) - beta;

  uint32_t upper[n], b[n], y[n];
  for(int i=1; i < n; ++i) upper[i] = 0;
  upper[0] = (~((a<<1)-1)+1)%(1<<BITSIZE);


  x[0] = (x[0] + a - 1)%DIL_Q;


  ConvertABModp(x, b, DIL_Q, BITSIZE, n);

  SecAdd(b, upper, y, BITSIZE, n);

  for(int i=0; i < n; ++i) res ^= (y[i] >> (BITSIZE-1));

  return res;


}



int old_rejection_sampling_shared_output(uint32_t* x, uint32_t* y, int mode, int n){
  /* Assumes input masked mod+ q, returns a Boolean masking of 1 if  value is valid
     Method converting A2B and looking at sign bit, used in [ABC+22]
  */


  const int BITSIZE = 24;
  const int beta = D_BETA;
  uint32_t a = 0, res = 0;
  if (mode == 0)  a = D_GAMMA1 - beta;
  else if (mode == 1) a = ((DIL_Q-1)/(2*DELTA)) - beta;

  uint32_t upper[n], b[n];
  for(int i=1; i < n; ++i) upper[i] = 0;
  upper[0] = (~((a<<1)-1)+1)%(1<<BITSIZE);


  x[0] = (x[0] + a - 1)%DIL_Q;


  ConvertABModp(x, b, DIL_Q, BITSIZE, n);

  SecAdd(b, upper, y, BITSIZE, n);

  for(int i=0; i < n; ++i)  (y[i] >>= (BITSIZE-1));

  return res;


}


void gen_y(uint32_t* y, int n){
  /* Generates a random y masked mod q in the interval [-2^{mu-1}, 2^{mu-1}) [CGTZ23] */
  int k=K_EXACT;
  uint64_t x[n];
  uint64_t arith_x[n];
  for(int i=0; i < n; ++i) x[i] = rand32()%(1<<MU);
  impconvBA64(arith_x, x, n);
  for(int i=0; i < n; ++i) arith_x[i] %= (1LLU<<k);
  exact_modulus_switching(arith_x, y, n);
  y[0] = (y[0] + DIL_Q - (1<<(MU-1)))%DIL_Q;
}




void gen_y_det(uint64_t* x, uint32_t* y, int n){
  /* Generates a y masked mod q in the interval [-2^{mu-1}, 2^{mu-1}) 
     from a Boolean masked value in {0,1}^mu
    */
  int k=K_EXACT;
  uint64_t arith_x[n];

  impconvBA64(arith_x, x, n);
  for(int i=0; i < n; ++i) arith_x[i] %= (1LLU<<k);

  exact_modulus_switching(arith_x, y, n);
  y[0] = (y[0] + DIL_Q - (1<<(MU-1)))%DIL_Q;

}


uint64_t portable_128_bit_mul_shift(uint64_t a, uint64_t b, int k){
  // Computes (a*b)>>k with 64-bit a and b. Result must be less than 2^64.
  // https://stackoverflow.com/questions/25095741/how-can-i-multiply-64-bit-operands-and-get-128-bit-result-portably

    uint64_t lo_lo = (a & 0xFFFFFFFF) * (b & 0xFFFFFFFF);
    uint64_t hi_lo = (a >> 32)        * (b & 0xFFFFFFFF);
    uint64_t lo_hi = (a & 0xFFFFFFFF) * (b >> 32);
    uint64_t hi_hi = (a >> 32)        * (b >> 32);

    /* Now add the products together. These will never overflow. */
    uint64_t cross = (lo_lo >> 32) + (hi_lo & 0xFFFFFFFF) + lo_hi;
    uint64_t upper = (hi_lo >> 32) + (cross >> 32)        + hi_hi;
    uint64_t lower = (cross << 32) | (lo_lo & 0xFFFFFFFF);

    return ((upper&((1<<k)-1)) << (64-k)) | (lower >> k);


}


void exact_modulus_switching(uint64_t* x, uint32_t* y, int n){
  //assume alpha = 4
  uint32_t k = K_EXACT;
  uint32_t a = A_EXACT;

  uint32_t z[n];
  
//k+23+ 23 = 45+46 = 91. 91- 41=50;
// size = alpha+log(a)+log(q);

  for(int i=0; i < n; ++i){

    #ifndef PORTABLE
    uint128_t temp  = (uint128_t)x[i]*a*DIL_Q;
    temp >>= (k-ALPHA);
    #else
    uint64_t temp = portable_128_bit_mul_shift(x[i], (uint64_t)a*DIL_Q , k-ALPHA);
    #endif
    z[i] = (uint32_t)(temp % DIL_Q_PRIME);

  }
  z[0] = (z[0] + n - 1)%DIL_Q_PRIME;
  generic_shift(z, y, ALPHA, DIL_Q, n);
}


void rejection_sampling(uint32_t* x, uint32_t* res, int mode, int n){
  const int beta = 196;
  int a = 0;
  if (mode == 0)      a = (1<<17) - beta;
  else if (mode == 1) a = ((DIL_Q-1)/32) - beta;
  const int p = 10687;
  const int l = 10; //p*2^l = 10943488 > q
  int rho = l + ALPHA; // 14 (p*2^rho is 28 bits)
  int modulus = p*(1<<rho);
  int u = ((int64_t)modulus*a/DIL_Q);


  uint32_t y[n], ny[n], positive_side[n], negative_side[n];
  y[0] = ((((int64_t)x[0]*modulus)/DIL_Q)+n-1)%modulus;
  for(int i=1; i < n; ++i) y[i] = (((int64_t)x[i]*modulus)/DIL_Q)%modulus;

  ny[0] = ((((int64_t)(DIL_Q-x[0])*modulus)/DIL_Q)+n-1)%modulus;
  for(int i=1; i < n; ++i) ny[i] = (((int64_t)(DIL_Q-x[i])*modulus)/DIL_Q)%modulus;


  y[0] = (y[0] + modulus - u)%modulus;


  generic_shift(y, positive_side, rho, p, n);



  ny[0] = (ny[0] + modulus - u)%modulus;
  generic_shift(ny, negative_side, rho, p, n);


  SecMultModp(positive_side, negative_side, res, p, n);

}




void decompose1(uint32_t* x, uint32_t* high, uint32_t* low, int n){
  const int RRHO = 25; //assume n <= 8.
  //const int DELTA = 16;
  

  #if ((DILITHIUM_MODE == 3) || (DILITHIUM_MODE == 5))
  uint32_t modulus = (DELTA << RHO);
  uint32_t y[n], z[n];

  //printArith(x, DIL_Q, n);

  for(int i=0; i < n; ++i) y[i] = ((((uint64_t)x[i]*DELTA)<<RRHO)/DIL_Q)%modulus;
  y[0] = (y[0] + n - 1 + (1<<(RHO-1)))%modulus;
  //printArith(y, modulus, n);
  generic_shift(y, z, RRHO, DELTA, n);
  //printArith(z, DELTA, n);
  refreshArithModp(z, DELTA, n);

  #else

  uint64_t modulus = (DELTA << RRHO);
  uint64_t y[n], z[n];

  for(int i=0; i < n; ++i) y[i] = ((((uint64_t)x[i]*DELTA)<<RRHO)/DIL_Q)%modulus;
  y[0] = (y[0] + n - 1 + (1<<(RRHO-1)))%modulus;
  generic_shift64(y, z, RRHO, DELTA, n);
  refreshArithModp64(z, DELTA, n);

  #endif


  *high = 0;
  for(int i=0; i < n; ++i) *high = (*high + z[i])%DELTA;
  
  for(int i=0; i < n; ++i) low[i] = x[i];
  low[0] = (low[0] + DIL_Q - *high*(DIL_Q-1)/DELTA)%DIL_Q;
  

}


void decompose2(uint32_t* x, uint32_t* high, uint32_t* low, int n){

  const int BITSIZE = 24;

  uint32_t s[n], b[n];
  #if (DILITHIUM_MODE == 2)
  uint32_t a[n], t[n], masked_high[n];
  #endif




  for(int i=0; i < n; ++i) s[i] = ((DIL_Q-DELTA)*(uint64_t)x[i])%DIL_Q;



  s[0] = (s[0] + (DIL_Q-1)/2)%DIL_Q;


  
  ConvertABModp(s, b, DIL_Q, BITSIZE, n);



  #if (DILITHIUM_MODE == 2)
  for(int j=0; j < n; ++j) t[j] = b[j]&1;
  bool2ArithSPOGmodq(t, masked_high, DELTA, n);
  //for(int i=0; i < n; ++i) masked_high[i] = a[i];
  
  for(int i=1; i < BITSIZE; ++i){
    for(int j=0; j < n; ++j) t[j] = (b[j] >> i)&1;
    bool2ArithSPOGmodq(t, a, DELTA, n);
    for(int j=0; j < n; ++j) masked_high[j] = (masked_high[j] + a[j] * (1<<i)%DELTA)%DELTA;
  }


  *high = 0;
  for(int i=0; i < n; ++i) *high = (*high + masked_high[i])%DELTA;
  #else
  *high = 0;
  for(int i=0; i < n; ++i) *high = (*high ^ (b[i]&0xF));
  #endif


  for(int i=0; i < n; ++i) low[i] = x[i];
  low[0] = (low[0] + DIL_Q - *high*(DIL_Q-1)/DELTA)%DIL_Q;
 
  
  
}




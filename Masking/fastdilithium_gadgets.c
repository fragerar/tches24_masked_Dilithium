#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "random.h"
#include "convtable.h"
#include "convba_2014.h"
#include "utils.h"
#include "dilithium_gadgets.h"


#ifdef COUNT
uint64_t count_rand = 0;
#endif




static int zero_test_mul(uint32_t* x, int q, int n){
  // returns 1 if x = 0 mod Q

  uint32_t u;
  for(int i=0; i < n; ++i){
    u = (rand32()%(q-1))+1;
    for(int j=0; j < n; ++j) x[j] = ((uint64_t)x[j]*u)%q;
    refreshArithModp(x, q, n);
  }


  uint32_t B = x[0];
  for(int i=1; i < n; ++i) B = (B+x[i])%q;
  
  if (B == 0) return 1;
  else        return 0;


}


int single_coeff_reject(uint32_t* x, int mode, int q, int n){
  // returns 1 if coefficient in range
  uint32_t y[n];
  reject_sampling(x, y, mode, q, n);
  return !zero_test_mul(y, q, n);

}

void MSB_estimation(uint32_t* x, uint32_t* y, int n){
  /* Algorithm 1 */
  int T[15] = {1,2,2,3,3,3,3,4,4,4,4,4,4,4,4};
  int alpha = T[n-1];
  int k = MU+alpha;

  uint32_t t[n], v[n];
  for(int i=0; i < n; ++i) t[i] = x[i] >> MU;
  t[0] = (t[0] + (1<<alpha) - 1)%(DIL_Q << alpha);
  generic_shift(t, v, alpha, DIL_Q, n);
  refreshArithModp(v, DIL_Q, n);
  for(int i=0; i < n ; ++i) y[i] = (x[i] + DIL_Q - ((uint64_t)v[i] << k)%DIL_Q)%DIL_Q;
  
}


void LMSwitch(uint32_t* x, uint32_t* y, int q, int rho, int n){
  // maps x mod q to x mod q*2^rho
  int T[15] = {1,2,3,3,4,4,4,4,5,5,5,5,5,5,5};
  int alpha = T[n-1];

  uint32_t delta[n];

  for(int i=0; i < n; ++i) y[i] = (x[i]<<(alpha))/q;
  y[0] = (y[0] + n-1+(1<<(alpha-2)))%(1<<(rho+alpha));



  generic_shift(y, delta, alpha, 1<<rho, n);
  refreshArithModp(delta, 1<<rho, n);

  for(int i=0; i < n; ++i) y[i] = ((uint64_t)x[i] + ((q<<rho) - q*delta[i]))%(q<<rho);

  
}


void reject_sampling(uint32_t* x, uint32_t* z, int mode, int q, int n){

  /* Algorithm 2*/
  uint32_t a = 0;

  uint32_t u[n], y[n], y_p[n], t1[n], t2[n];

  if (mode == 0){ //reject z
    a = D_GAMMA1 - D_BETA;
  } else {
    a = D_GAMMA2 - D_BETA;
  }

  

  LMSwitch(x, u, q, RHO, n);
  for(int i=1; i < n; ++i) y[i] = u[i]; 
  y[0] = (u[0] - a)%(q<<RHO);
  generic_shift(y, t1, RHO, q, n);
  for(int i=0; i < n; ++i) y_p[i] = ((q<<RHO) - u[i]);
  y_p[0] = (y_p[0] + (q<<RHO) - a)%(q<<RHO);
  refreshArithModp(y_p, q<<RHO, n);
  generic_shift(y_p, t2, RHO, q, n);
  SecMultModp(t1, t2, z, q, n);
  

}


void gen_y_fast(uint32_t* y, int n){
  /* Generates a random y masked mod q in the interval [-2^{mu-1}, 2^{mu-1}) */
  int k=MU+ALPHA;
  uint32_t x[n];
  uint32_t arith_x[n];

  for(int i=0; i < n; ++i) x[i] = rand32()%(1<<MU);
  impconvBA(arith_x, x, n);
  for(int i=0; i < n; ++i) arith_x[i] %= (1<<k);
  MSB_estimation(arith_x, y, n);
  y[0] = (y[0] + DIL_Q - (1<<(MU-1)))%DIL_Q;

}



int early_abort_RS(uint32_t* masked_poly, int k, int l, int q, int n){
  // z is polyvec_l and r0 is polyvec_k. Masked_poly has thus 256*(l+k)*n coefficients. In case of dilithium3, l=5 k=6.
  // Algorithm 4
  int i=0;
  int b=1;
  while (b && (i < 256*k*n)){
    b = single_coeff_reject(masked_poly+i, 0, q, n);
    i += n;
  }

  while (b && (i < 256*(l+k)*n)){
    b = single_coeff_reject(masked_poly+i, 1, q, n);
    i += n;
  }
  return b;

}

int early_abort_old_RS(uint32_t* masked_poly, int k, int l, int n){
  int i=0;
  int b=1;
  while (b && (i < 256*k*n)){
    b = old_rejection_sampling(masked_poly+i, 0, n);
    i += n;
  }

  while (b && (i < 256*(l+k)*n)){
    b = old_rejection_sampling(masked_poly+i, 1, n);
    i += n;
  }
  return b;
}

int no_abort_RS(uint32_t* masked_poly, int k, int l, int q, int n){
  // Algorithm 5
  uint32_t res[n], temp[n];
  int i=0;

  res[0] = 1;
  for(int i=1; i < n; ++i) res[i] = 0;

  while(i < 256*k*n){
    reject_sampling(masked_poly+i, temp, 0, q, n);
    secMulAssignment(res, temp, q, n);
    i += n;
  }

  while(i < 256*(l+k)*n){
    reject_sampling(masked_poly+i, temp, 1, q, n);
    secMulAssignment(res, temp, q, n);
    i += n;
  }
  return !zero_test_mul(res, q, n);
}


int no_abort_old_RS(uint32_t* masked_poly, int k, int l, int n){
  // Algorithm 5
  uint32_t res[n], temp[n];
  int i=0;

  res[0] = 1;
  for(int i=1; i < n; ++i) res[i] = 0;

  while(i < 256*k*n){
    old_rejection_sampling_shared_output(masked_poly+i, temp, 0, n);
    secMulAssignment(res, temp, DIL_Q, n);
    i += n;
  }

  while(i < 256*(l+k)*n){
    old_rejection_sampling_shared_output(masked_poly+i, temp, 1, n);
    secMulAssignment(res, temp, DIL_Q, n);
    i += n;
  }
  return !zero_test_mul(res, DIL_Q, n);
}


#ifdef TEST_FASTDIL

void bench_MSB_conv(){

  uint64_t start, stop;
  int max_n = 8;
  uint32_t x[max_n], y[max_n];
  uint64_t x64[max_n];
  uint32_t val;
  
  val = rand()%(1<<MU);

  int ITER=100000;

  printf("Avg speed exact_mod_switch: ");
  for(int j=2; j < max_n; ++j){
    share64((uint64_t)val, x64, j);
    refreshArith64(x64, K_EXACT, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) exact_modulus_switching(x64, y, j);
    stop = cpucycles();
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");


  printf("Avg speed MSB_conversion: ");
  for(int j=2; j < max_n; ++j){
    share(val, x, j);
    refreshArith(x, 22, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) MSB_estimation(x, y, j);
    stop = cpucycles();
    //printf(" %f (%i) |", (double)(stop-start)/(ITER), j);
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");
  printf("\n");


}



void bench_gen_y(){
  
  uint64_t start, stop;
  int max_n = 8;
  uint32_t y[max_n];
  



  int ITER=100000;


  printf("Avg speed gen_y exact_mod_switch: ");
  for(int j=2; j < max_n; ++j){
    start = cpucycles();
    for(int i=0; i < ITER; i++) gen_y(y, j);
    stop = cpucycles();
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");


  printf("Avg speed gen_y LMSwitch: ");
  for(int j=2; j < max_n; ++j){
    start = cpucycles();
    for(int i=0; i < ITER; i++) gen_y_fast(y, j);
    stop = cpucycles();
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");



}



void bench_rejection_sampling(){


  uint64_t start, stop;
  int max_n = 11;
  uint32_t x[max_n];
  uint32_t val;
  

  val = rand()%(DIL_Q);


  int ITER=10000;


  printf("Avg speed Reject[ABC+22]: ");
  for(int j=2; j < max_n; ++j){
    share(val, x, j);
    refreshArith(x, DIL_Q, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) old_rejection_sampling(x, 0, j);
    stop = cpucycles();
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");
  printf("Avg speed rejection_sampling: ");
  for(int j=2; j < max_n; ++j){
    share(val, x, j);
    refreshArith(x, DIL_Q, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) single_coeff_reject(x, 0, DIL_Q, j);
    stop = cpucycles();
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");
  printf("\n");





}

static int centered(uint32_t* x, int q, int n){
  // Returns unmasked value in centered representation mod q 
  int res = 0;
  for(int i=0; i < n; ++i) res = (res + x[i])%q;
  
  if (res >= q/2) res -= q;
  return res;

}




int read_sign_vectors(const char* file_name, uint32_t* out, int size, int n)
{
  // Size is the number of values to read, we assume that it is less or equal to the total number of values in the file.
  FILE* file = fopen (file_name, "r");
  int coef = 0;
  int i = 0;
  uint32_t x[n];
  int c=0;
  
  c += fscanf (file, "%d", &coef);    
  while ((!feof (file)) && (i < size))
    {  
      share((uint32_t)coef, x, n);
      refreshArithModp(x, DIL_Q, n);
      for(int j=0; j<n; ++j) out[i*n+j] = x[j];
      c += fscanf (file, "%d", &coef);  

      i++;
    }
  fclose (file);
  return c;
}



void bench_no_abort(){
  uint64_t start, stop;
  int ITER=30;
  int max_n = 10;
  int k=4, l=4;
  int size=256*(k+l);
  uint32_t* x;
  x =  (uint32_t*) malloc(sizeof(uint32_t)*(size*max_n*ITER));

  for(int i=0; i < size*max_n*ITER; ++i) x[i] = 0;
  printf("Init done\n");
  

  printf("Avg speed no abort reject: ");
  for(int j=2; j < max_n; ++j){

    read_sign_vectors("./ZR.txt", x, size*ITER, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) no_abort_RS(x+size*j*i, k, l, DIL_Q, j);
    stop = cpucycles();
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");


  printf("Avg speed no abort old reject [ABC22+]: ");
  for(int j=2; j < max_n; ++j){

    read_sign_vectors("./ZR.txt", x, size*ITER, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) no_abort_old_RS(x+size*j*i, k, l, j);
    stop = cpucycles();
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");
}




void bench_abort(){

  uint64_t start, stop;
  int ITER=30;
  int max_n = 10;
  int k=4, l=4;
  int size=256*(k+l);
  uint32_t* x;
  x =  (uint32_t*) malloc(sizeof(uint32_t)*(size*max_n*ITER));

  for(int i=0; i < size*max_n*ITER; ++i) x[i] = 0;
  printf("Init done\n");
  
  printf("Avg speed abort reject: ");
  for(int j=2; j < max_n; ++j){
    read_sign_vectors("./ZR.txt", x, size*ITER, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) early_abort_RS(x+size*j*i, k, l, DIL_Q, j);
    stop = cpucycles();
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");


  printf("Avg speed abort old reject [ABC22+]: ");
  for(int j=2; j < max_n; ++j){
    read_sign_vectors("./ZR.txt", x, size*ITER, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) early_abort_old_RS(x+size*j*i, k, l, j);
    stop = cpucycles();
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");

  //early_abort_RS(uint32_t* masked_poly, int k, int l, int q, int n)

}


static int test_MSB_conv(int n){
  uint32_t x[n], y[n];
  uint32_t val, res;

  int ITER=100000;


  printf("Test MSB_conversion: ");
  for(int i=0; i < ITER; i++){
    val = rand()%(1<<MU);
    share(val, x, n);
    refreshArith(x, 22, n);
    MSB_estimation(x, y, n);
    res = addopmodp(y, DIL_Q, n);
    if (res != val){
      printf("Fail! Iteration %i: (val, res) = (%i, %i)\n",i, val, res);
      return 0;
    }
  }
  printf("Success!\n");

  return 1;

}


static int test_mod_switch(int n){
  printf("Test mod switch in range [-q/4; q/4] at order %i: ", n-1);
  uint32_t x[n], y[n];



  for(int i=0; i < DIL_Q/4; ++i){
    int val = i;
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);  
    LMSwitch(x, y, DIL_Q, RHO, n);

    int t = centered(y, DIL_Q<<RHO, n);
    if (t != val){
      printf("i: %i => %i\n", i, t);
      return 0;
    }

  } 

  for(int i=DIL_Q*3/4; i < DIL_Q; ++i){
    int val = i;
    int centered_val = val - DIL_Q;
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);  
    LMSwitch(x, y, DIL_Q, RHO, n);

    int t = centered(y, DIL_Q<<RHO, n);
    if (t != centered_val){
      printf("i: %i => %i\n", i, t);
      return 0;
    }

  } 
  printf("Success\n");
  return 1;
}


static int test_RS_z(int n){
  printf("Test rejection sampling for z at order %i: ", n-1);
  const int beta = D_BETA;
  const int b = D_GAMMA1 + beta;
  const int bound = D_GAMMA1 - beta;
  //const int p = 10687;
  uint32_t x[n]; 
  int t;
  //printf("Bound: %d\n", bound);
  for(int val=0; val < b+1; val++){
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);
    t = single_coeff_reject(x, 0, DIL_Q, n);

    if ((t == 0) && (val < bound)){
      printf("Fail. Positive Valid value rejected. Value: %i\n", val);
      return 0;
    }
    if ((t != 0) && (val >= bound)){
      printf("Fail. Positive Invalid value accepted. Value: %i\n", val);
      return 0;
    }
  } 


  for(int val = DIL_Q-1; val > DIL_Q - b - 1; --val){
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);
    t = single_coeff_reject(x, 0, DIL_Q, n);
    if ((t == 0) && (val > DIL_Q - bound)){
      printf("Fail. Negative Valid value rejected. Value: %i\n", val);
      return 0;
    }
    if ((t != 0) && (val <= DIL_Q - bound)){
      printf("Fail. Negative Invalid value accepted. Value: %i\n", val);
      return 0;
    }
  }
  printf("Success\n");

  return 1;

}


static int test_RS_r(int n){
  printf("Test rejection sampling for r at order %i: ", n-1);
  const int beta = D_BETA;
  const int b = (DIL_Q-1)/(2*DELTA) + beta;
  const int bound = (DIL_Q-1)/(2*DELTA) - beta;
  uint32_t x[n];
  int t;
  //printf("Bound: %d\n", bound);
  for(int val=0; val < b+1; val++){
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);

    t = single_coeff_reject(x, 1, DIL_Q, n);
    if ((t == 0) && (val < bound)){
      printf("Fail. Valid value rejected. Value: %i\n", val);
      return 0;
    }
    if ((t != 0) && (val >= bound)){
      printf("Fail. Invalid value accepted. Value: %i\n", val);
      return 0;
    }
  } 

 

  for(int val = DIL_Q-1; val > DIL_Q - b - 1; --val){
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);
    t = single_coeff_reject(x, 1, DIL_Q, n);
    if ((t == 0) && (val > DIL_Q - bound)){
      printf("Fail. Valid value rejected. Value: %i\n", val);
      return 0;
    }
    if ((t != 0) && (val <= DIL_Q - bound)){
      printf("Fail. Invalid value accepted. Value: %i\n", val);
      return 0;
    }
  }

  printf("Success\n");

  return 1;

}


static int test_gen_y(int n){
  printf("Test gen y at order %i: ", n-1);
  uint32_t y[n];

  for(int i=0; i < 1000000; ++i){
    gen_y(y, n);
    if ((addopmodp(y, DIL_Q, n) >= (1<<(MU-1))) && (addopmodp(y, DIL_Q, n) < (DIL_Q - (1<<(MU-1))))){
      printf("Fail\n");
      return 0;
    }
  }
  printf("Success\n");
  return 1;

}



void test_all(){
  for(int i=2; i < 10; ++i){
    printf("Tests at order: %i\n", i-1);
    printf("\n");
    test_MSB_conv(i);
    test_mod_switch(i);
    test_RS_z(i);
    test_RS_r(i);
    test_gen_y(i);
    printf("\n\n");
  }
}

void bench_all(){
  bench_MSB_conv();
  bench_rejection_sampling();
  bench_no_abort();
  bench_abort();
}

int main(){
  srand(0);
  printf("Hello world\n");
  printf("Dilithium mode: %d\n", DILITHIUM_MODE);

  //test_all();
  bench_all();

  return 0;
}
#endif
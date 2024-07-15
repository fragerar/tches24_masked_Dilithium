#include "masking_interface.h"
#include "../Masking/dilithium_gadgets.h"
#include "params.h"
#include "polyvec.h"
#include "fips202.h"





int32_t canon_to_centered(uint32_t x){
  int32_t res = (int32_t) x;
  if (res >= Q/2) res -= Q;
  return res;
}

int32_t center(int32_t x){
  /* Maps elements of Z_q in [-Q, ..., Q] to representatives in [-Q/2, ..., Q/2[ */
  x += Q;
  x %= Q;
  if (x > Q/2) x -= Q;
  return x;
}

uint32_t center_to_canon(int32_t x){
  /* Maps elements of Z_q in [-Q/2, ..., Q/2[ to representatives in  [0, ..., Q[*/
  if (x < 0) x += Q; 
  return x;
}



void masked_sample_y(masked_polyvecl* masked_y){
  /* Randomized version of Dilithium. */
  polyvecl* y;
  uint32_t masked_coeff[N_SHARES];

  for(int i=0; i < L; ++i){
    for(int j=0; j < N; ++j){
      #ifdef FAST
      gen_y_fast(masked_coeff, N_SHARES);
      #else
      gen_y(masked_coeff, N_SHARES);
      #endif
      for(int k=0; k < N_SHARES; ++k){
          y = &(masked_y->shares[k]);
          y->vec[i].coeffs[j] = canon_to_centered(masked_coeff[k]);
      }
    }
  }
}


void masked_sample_y_MLDSA(masked_polyvecl* masked_y, uint8_t* seed, uint16_t mu){
  /* Correspond to expand_mask in MLDSA */
  int out_size = 32*mu;
  uint8_t out[out_size*N_SHARES];
  int seed_size=64; //see FIPS 204
  int in_size = seed_size + 2; 
  uint8_t in[(in_size)*N_SHARES];
  uint32_t masked_coeff_modq[N_SHARES];
  uint64_t masked_coeff[N_SHARES];
  
  uint16_t y_raw_coeffs[(out_size/2)*N_SHARES];
  for(int i=0; i < N_SHARES; ++i){
    for(int j=0; j < seed_size; ++j){
      in[j+in_size*i] = seed[j+seed_size*i];
    }
    in[seed_size   + in_size*i] = 0;
    in[seed_size+1 + in_size*i] = 0;
  }
  for(uint16_t r=0; r < L; ++r){
    uint16_t n = r+mu;
    in[seed_size+1] = (uint8_t)(n>>8);
    in[seed_size  ] = (uint8_t)(n&0xFFFF);

    shake256_masked(out, out_size, in, in_size); //32MU bytes = 256MU bits 
 

    for(int i=0; i < out_size*N_SHARES; i+=2 ) y_raw_coeffs[i/2] = (((uint16_t)(out[i])) << 8) | ((uint16_t)(out[i+1]));
 

    for(int i=0; i < N; ++i){
 
      int offset = (16 - ((i+1)*mu)%16)%16;
      for(int j=0; j < N_SHARES; ++j){
 
        uint16_t a = y_raw_coeffs[(i*mu/16)+0 + j*out_size/2]; // outsize/2 because bits are packed by 16 and not 8.
        uint16_t b = y_raw_coeffs[(i*mu/16)+1 + j*out_size/2];
        uint32_t share = ((((uint32_t)a)<<(16-offset)) | ((uint32_t)(b>>offset)))%(1LU<<mu);
        
        masked_coeff[j] = (uint64_t)share;

      }
      gen_y_det(masked_coeff, masked_coeff_modq, N_SHARES);

      for(int j=0; j < N_SHARES; ++j) (masked_y->shares[j]).vec[r].coeffs[i] = canon_to_centered(masked_coeff_modq[j]);
    }

  }

}



int masked_rejection_sampling_z(masked_polyvecl* mz){
  /* Return 1 if reject */
  uint32_t temp[N_SHARES];

  for(int i=0; i < L; ++i){
    for(int j=0; j < N; ++j){
      for(int k=0; k < N_SHARES; ++k){
        temp[k] = center_to_canon(mz->shares[k].vec[i].coeffs[j]);
      }
    #ifdef FAST
    if ((single_coeff_reject(temp, 0, DIL_Q, N_SHARES)) == 0) return 1;
    #else
    if (old_rejection_sampling(temp, 0, N_SHARES) == 0) return 1;
    #endif
    }
  }
  return 0;
}

int masked_rejection_sampling_r(masked_polyveck* mr){
  /* Should return 1 if reject */
  uint32_t temp[N_SHARES];

  for(int i=0; i < K; ++i){
    for(int j=0; j < N; ++j){
      for(int k=0; k < N_SHARES; ++k){
        temp[k] = center_to_canon(mr->shares[k].vec[i].coeffs[j]);
      }
    #ifdef FAST
    if ((single_coeff_reject(temp, 1, DIL_Q, N_SHARES)) == 0) return 1;
    #else
    if (old_rejection_sampling(temp, 1, N_SHARES) == 0) return 1;
    #endif
    }
  }
  return 0;
}





void masked_decompose(polyveck* r1, masked_polyveck* mr0, masked_polyveck* mr){


  uint32_t temp_r[N_SHARES], temp_r0[N_SHARES];


  for(int i=0; i < K; ++i){
    for(int j=0; j < N; ++j){
      for(int k=0; k < N_SHARES; ++k) temp_r[k] = mr->shares[k].vec[i].coeffs[j];
      decompose2(temp_r, &(r1->vec[i].coeffs[j]), temp_r0, N_SHARES);
      for(int k=0; k < N_SHARES; ++k) mr0->shares[k].vec[i].coeffs[j] = temp_r0[k];
    }
  } 


}



void mask_bitstring(uint8_t* bs_out, uint8_t* bs_in, int len){
  uint8_t temp;
  for(int i=0; i < len; ++i){
    bs_out[i] = bs_in[i]; 
    for(int j=1; j < N_SHARES; ++j){
      temp = rand32()&0xFF;
      bs_out[i + len*j]  = temp;
      bs_out[i        ] ^= temp;

    }
  }

}


void unmask_polyvecl(masked_polyvecl* mpv, polyvecl* pv){
  int32_t temp;
  for(int i=0; i < L; ++i){
    for(int j=0; j < N; ++j){
      temp = 0;
      for(int k=0; k < N_SHARES; ++k) temp = (temp + (mpv->shares[k]).vec[i].coeffs[j])%Q;
      pv->vec[i].coeffs[j] = canon_to_centered((temp+2*Q)%Q);
    }
  }
}


void unmask_polyveck(masked_polyveck* mpv, polyveck* pv){
  int32_t temp;
  for(int i=0; i < K; ++i){
    for(int j=0; j < N; ++j){
      temp = 0;
      for(int k=0; k < N_SHARES; ++k) temp = (temp + (mpv->shares[k]).vec[i].coeffs[j])%Q;
      pv->vec[i].coeffs[j] = canon_to_centered((temp+2*Q)%Q);
    }
  }
}



void mask_polyvecl(masked_polyvecl* mpv, polyvecl* pv){
 /* Takes an unmasked polyvecl and mask it (arith mod q), should not be used outside of testing*/
  int32_t temp;
  for(int i=0; i < L; ++i){
    for(int j=0; j < N; ++j){
      ((mpv->shares[0]).vec[i]).coeffs[j] = (pv->vec[i]).coeffs[j];
    }
  }

  for(int k=1; k < N_SHARES; ++k){
    for(int i=0; i < L; ++i){
      for(int j=0; j < N; ++j){
        temp = ((int32_t)rand32())%Q;
        ((mpv->shares[k]).vec[i]).coeffs[j] = center(temp);
        ((mpv->shares[0]).vec[i]).coeffs[j] = center((((mpv->shares[0]).vec[i]).coeffs[j] - temp)%Q);
      }
    }
  }
}


void mask_polyveck(masked_polyveck* mpv, polyveck* pv){
 /* Takes an unmasked polyveck and mask it (arith mod q), should not be used outside of testing*/
  int32_t temp;
  for(int i=0; i < K; ++i){
    for(int j=0; j < N; ++j){
      ((mpv->shares[0]).vec[i]).coeffs[j] = (pv->vec[i]).coeffs[j];
    }
  }

  for(int k=1; k < N_SHARES; ++k){
    for(int i=0; i < K; ++i){
      for(int j=0; j < N; ++j){
        temp = ((int32_t)rand32())%Q;
        ((mpv->shares[k]).vec[i]).coeffs[j] = center(temp);
        ((mpv->shares[0]).vec[i]).coeffs[j] = center((((mpv->shares[0]).vec[i]).coeffs[j] - temp)%Q);
      }
    }
  }
}

void print_masked_polyvecl(masked_polyvecl* mpv){
  polyvecl t;
  unmask_polyvecl(mpv, &t);
  for(int i=0; i < L; ++i){
    printf("Pv[%i]: ", i);
    for(int j=0; j < 8; ++j) printf("%i ", t.vec[i].coeffs[j]);
    printf("\n");
  }
  printf("\n");
}

void print_polyvecl(polyvecl* t){
  for(int i=0; i < L; ++i){
    printf("Pv[%i]: ", i);
    for(int j=0; j < 8; ++j) printf("%i ", t->vec[i].coeffs[j]);
    printf("\n");
  }
  printf("\n");
}

#ifndef MASKING_INTERFACE_H
#define MASKING_INTERFACE_H

#ifndef MASKING_ORDER
#define MASKING_ORDER 2
#endif
#define N_SHARES (MASKING_ORDER + 1)

#include <stdint.h>
#include "polyvec.h"
#include "fips202.h"

#define FAST
//#define MLDSA_Y



typedef struct masked_polyvecl {polyvecl shares[N_SHARES];} masked_polyvecl;
typedef struct masked_polyveck {polyveck shares[N_SHARES];} masked_polyveck;

int32_t canon_to_centered(uint32_t x);
int32_t center(int32_t x);
void masked_sample_y(masked_polyvecl* masked_y);
void masked_decompose(polyveck* r1, masked_polyveck* r0, masked_polyveck* r);
int masked_rejection_sampling(masked_polyvecl* mz, masked_polyveck* mr);
int masked_rejection_sampling_z(masked_polyvecl* mz);
int masked_rejection_sampling_r(masked_polyveck* mr);


void unmask_polyvecl(masked_polyvecl* mpv, polyvecl* pv);
void unmask_polyveck(masked_polyveck* mpv, polyveck* pv);

void mask_bitstring(uint8_t* bs_out, uint8_t* bs_in, int len);
void mask_polyvecl(masked_polyvecl* mpv, polyvecl* pv);
void mask_polyveck(masked_polyveck* mpv, polyveck* pv);

void print_masked_polyvecl(masked_polyvecl* mpv);
void print_polyvecl(polyvecl* mpv);



// Paper Fast

int fast_no_abort_masked_rejection_sampling(masked_polyvecl* mz, masked_polyveck* mr);
int fast_early_abort_masked_rejection_sampling(masked_polyvecl* mz, masked_polyveck* mr);

void masked_sample_y_MLDSA(masked_polyvecl* masked_y, uint8_t* seed, uint16_t mu);


#endif

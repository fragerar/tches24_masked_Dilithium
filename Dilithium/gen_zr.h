#ifndef GEN_ZR
#define GEN_ZR

#include <stdio.h>

int crypto_sign_zr_print(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk);

#endif

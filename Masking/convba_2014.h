#ifndef CONVBA2014_H
#define CONVBA2014_H

void ConvertAB(uint32_t *A,uint32_t *z,int k,int n);
void ConvertABModp(uint32_t *A,uint32_t *z,uint32_t p,int k,int n);
void ConvertBA(uint32_t *x,uint32_t *A,int k,int n);
void ConvertBA_SPOG(uint32_t *x,uint32_t *y,int n);
void shift(uint32_t *x,uint32_t *y,int kin,int ell,int n);
void SecMul(uint32_t *a,uint32_t *b,uint32_t *c,int n);


void impconvBA(uint32_t *D_,uint32_t *x,int n); //[BCZ18]
void impconvBA64(uint64_t *D_,uint64_t *x,int n);
#endif
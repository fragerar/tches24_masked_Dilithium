#ifndef CONVTABLE_H
#define CONVTABLE_H


void arith2BoolOpti3NI(uint32_t *x,uint32_t *y,int k,int n);
void bool2ArithSPOGmodqMulti(uint32_t *x,uint32_t *y,int l,int q,int n);
void bool2ArithSPOGmodqMulti64(uint64_t *x, uint64_t *y, int l, uint64_t q, int n);
void bool2ArithSPOGmodq64(uint64_t *x, uint64_t *y, uint64_t q, int n);
void bool2ArithSPOGmodq(uint32_t *x,uint32_t *y,int q,int n);


#endif


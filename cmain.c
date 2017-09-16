#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include "ymm_mul.h"

void combine_29bit(unsigned int* data, unsigned int* result);
void split_29bit(unsigned int* data, unsigned int* result, int digits);
void multiply(unsigned int* a, unsigned int* b, unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int digits, int digitstimestwo);
void BigMultiply(unsigned int* A, unsigned int* B, unsigned int* T);


int main(){

    unsigned int* data_a = calloc((N+4), sizeof(int));
    unsigned int* data_b = calloc((N+4), sizeof(int));
    unsigned int* A = calloc((N+4), sizeof(int));
    unsigned int* B = calloc((N+4), sizeof(int));
    unsigned int* T = calloc((2*N+4), sizeof(int));
    unsigned int* result_t = calloc((2*N+4), sizeof(int));
	
    unsigned int* a = calloc((M+8),sizeof(int));
    unsigned int* b = calloc((M+8),sizeof(int));
    unsigned int* t = calloc((2*M+8),sizeof(long));
    unsigned int* u = calloc((2*M+8),sizeof(long));
    unsigned int* v = calloc((2*M+8),sizeof(long));
    unsigned int* w = calloc((2*M+8),sizeof(long));
    

    if(!(result_t && data_a && data_b && a && b && t && u && v && w && A && B && T)){
		puts("malloc error");
		exit(1);
    }

    int i,j;

    // initialize
    for(j=0; j<N-1; j++){
		data_a[j] = A[j] = 0x11111111;
		data_b[j] = B[j] = 0x11111111;
    }
    
    
    
    split_29bit(data_a, a, MDIGITS);
    for(j=0; j<M+1; j++) printf("%d : %x\n", j, a[j]);
    //split_29bit(data_b, b, MDIGITS);
    //for(j=0; j<M+1; j++) printf("%d : %x\n", j, b[j]);

    clock_t s, e;
    double total = 0.0;
    double time = 0.0;
    
    // optimized multiply
    
    for(i=0; i<30; i++){
	
	    for(j=0; j<2*(2*M+1); j++) t[j] = u[j] = v[j] = w[j]  = 0;

    	s = clock();
	    split_29bit(data_a, a, MDIGITS);
    	split_29bit(data_b, b, MDIGITS);

	    multiply(a, b, t, u, v, w, DIGITS, DIGITSTIMESTWO);

    	combine_29bit(t, result_t);
	    e = clock();
        //printf("%lf\n", (double)(e-s)/CLOCKS_PER_SEC);
	    total += (double)(e-s);

    }
    
    time = total / 30.0;
    printf("\nOptimized : Average time = %lf [ms]\n", (time/CLOCKS_PER_SEC)*1000);
    total = 0.0;
    
	for(j=0; j<2*(2*M+1); j++) t[j] = u[j] = v[j] = w[j]  = 0;

    split_29bit(data_a, a, MDIGITS);
    split_29bit(data_b, b, MDIGITS);

    multiply(a, b, t, u, v, w, DIGITS, DIGITSTIMESTWO);
    
    combine_29bit(t, result_t);
    
    // normal multiply 
    for(i=0; i<30; i++){

        for(j=0; j<(2*N+1); j++) T[j] = 0;

    	s = clock();
	    BigMultiply(A, B, T);
    	e = clock();
	    //printf("%lf\n", (double)(e-s)/CLOCKS_PER_SEC);
	    total += (double)(e-s);

    }
    time = total / 30.0;
    printf("\nNormal : Average time = %lf [ms]\n", (time/CLOCKS_PER_SEC)*1000);
    
    for(j=0; j<(2*N+1); j++) T[j] = 0;

    BigMultiply(A, B, T);
    // Error check
    int flag = 0;
    for(i=2*N; i>=0; i--){
	if(result_t[i] != T[i]){
	    printf("%d\n%x %x\n", i, result_t[i], T[i]);
	    flag = 1;
	}
    }
    if(!flag) puts("-------- No Error --------");

    //printf("digit = %d\ndigittimestwo = %d\n", DIGITS, DIGITSTIMESTWO);
/*
    //for(i=0; i<2*N+2; i++) printf("%d : %08x\n", i, result_t[i]);
    for(i=2*N; i>=0; i--)  printf("%08x ", result_t[i]); puts("");
    for(i=2*N; i>=0; i--)  printf("%08x ", T[i]); puts("");
*/  
    
    free(a);
    free(b);
    free(t);
    free(u);
    free(v);
    free(w);
    free(A);
    free(B);
    free(T);
    free(data_a);
    free(data_b);
    free(result_t);
    
    return 0;
}

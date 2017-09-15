#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#define BITSIZE 128 * 2
#define N BITSIZE/32 //dataの配列数
#define M BITSIZE/29 //resultの配列数を決める指標
//#define R BITSIZE - 58 * (BITSIZE / 58) //余りのビット数
#define MDIGITS 2*(BITSIZE/116) // 116 = 29*4. used in split_29bit

#define DIGITS M+1
#define DIGITSTIMESTWO ((M+1)*2)+1

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


void cloop(unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int i){

	__asm__ volatile(
	"mov %0, %%ebx;" //i
	"movq %1, %%rcx;" //t
	"movq %2, %%rdx;" //u
	"movq %3, %%r8;" //v
	"movq %4, %%r9;" //w
	::"r"(i), "r"(t), "r"(u), "r"(v), "r"(w)
	);

	__asm__ volatile(
	"movq (%rcx, %rbx, 8), %mm0;"
	"movq -8(%rdx, %rbx, 8), %mm1;"
	"movq -16(%r8, %rbx, 8), %mm2;"
	"movq -24(%r9, %rbx, 8), %mm3;"
	"paddq %mm1, %mm0;"
	"paddq %mm2, %mm0;"
	"paddq %mm3, %mm0;"
	"paddq %mm6, %mm0;"
	"movq %mm0, %mm6;"
	"pand %mm7, %mm6;"
	"movd %mm6, (%rcx, %rbx, 4);"
	"psrlq $29, %mm0;"
	

	"movq 8(%rcx, %rbx, 8), %mm6;"
	"movq (%rdx, %rbx, 8), %mm1;"
	"movq -8(%r8, %rbx, 8), %mm2;"
	"movq -16(%r9, %rbx, 8), %mm3;"
	"paddq %mm1, %mm6;"
	"paddq %mm2, %mm6;"
	"paddq %mm3, %mm6;"
	"paddq %mm0, %mm6;"
	"movq %mm6, %mm0;"
	"pand %mm7, %mm0;"
	"movd %mm0, 4(%rcx, %rbx, 4);"
	"psrlq $29, %mm6;"	

	);

}


void calc_carry(unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int digitstimestwo){


	/*
	  in this routine,
	  mm7 : 0x1fffffff
	  mm6 : carry
	  mm0 : t
	  mm1 : u
	  mm2 : v
	  mm3 : w

	  memo
	  "movq -24(%rbp), %rcx;" //t
	  "movq -32(%rbp), %rdx;" //u
	  "movq -40(%rbp), %r8;" //v
	  "movq -48(%rbp), %r9;" //w
	 */

	__asm__ volatile(
	"movq %0, %%rcx;" //t
	"movq %1, %%rdx;" //u
	"movq %2, %%r8;" //v
	"movq %3, %%r9;" //w
	::"r"(t), "r"(u), "r"(v), "r"(w)
	);

	 __asm__ volatile(
	"pxor %mm6, %mm6;"
	"mov $0x1fffffff, %eax;"
	"movd %eax, %mm7;"
	//[0]
	"movq (%rcx), %mm0;"
	"movq %mm0, %mm6;"
	"pand %mm7, %mm0;"
	"movd %mm0, (%rcx);"
	"psrlq $29, %mm6;"
	
	"mov $1, %rbx;"
	//[1]
	"movq (%rcx, %rbx, 8), %mm0;"
	"movq -8(%rdx, %rbx, 8), %mm1;"
	"paddq %mm1, %mm0;"
	"paddq %mm6, %mm0;"
	"movq %mm0, %mm6;"
	"pand %mm7, %mm6;"
	"movd %mm6, (%rcx, %rbx, 4);"
	"psrlq $29, %mm0;"

	"inc %ebx;"
	//[2]
	"movq (%rcx, %rbx, 8), %mm6;"
	"movq -8(%rdx, %rbx, 8), %mm1;"
	"movq -16(%r8, %rbx, 8), %mm2;"
	"paddq %mm1, %mm6;"
	"paddq %mm2, %mm6;"
	"paddq %mm0, %mm6;"
	"movq %mm6, %mm0;"
	"pand %mm7, %mm0;"
	"movd %mm0, (%rcx, %rbx, 4);"
	"psrlq $29, %mm6;"
	
	//"inc %ebx;"
	);

	int i;
	for(i=3; i<digitstimestwo; i+=2){
		cloop(t, u, v, w, i);
	}
	
	// %rbxを引き継いでいるようだが、怖いなら%rbx = i - 3をやってもいいかも
	__asm__ volatile(
	"sub $1, %rbx;"
	"movd %mm6, (%rcx, %rbx, 4);"
	);
}


void multiply_outer(unsigned int* a, unsigned int* b, unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int i){


	__asm__ volatile(
	//"mov %0, %%ebx;" //i
	"movq %0, %%rax;" //a
	"movq %1, %%rbx;" //b
	"movq %2, %%rcx;" //t
	"movq %3, %%rdx;" //u
	"movq %4, %%r8;" //v
	"movq %5, %%r9;" //w
	::"a"(a), "b"(b), "c"(t), "d"(u), "r"(v), "r"(w)
	);

	__asm__ volatile(

    "movq $0, %rax;"
"outer:"
	"push %rbx;"
	);

	__asm__ volatile(
	"movq -16(%rbp), %rbx;"
	);

	__asm__ volatile(
	// load b[i] ... b[i+7]
	"vmovdqu (%rbx, %rax, 4), %ymm1;"
	"pop %rbx;"
	// load t[8-15], u[8-15], u[0-7], t[0-7], respectively
	"vmovdqu 32(%r8, %rax, 8), %ymm4;" //latter v
	"vmovdqu 32(%r9, %rax, 8), %ymm5;" //latter w
	"vmovdqu 32(%rdx, %rax, 8), %ymm6;" //latter u
	"vmovdqu 32(%rcx, %rax, 8), %ymm7;" //latter t
	"vmovdqu (%r8, %rax, 8), %ymm8;" //former v
	"vmovdqu (%r9, %rax, 8), %ymm9;" //former w
	"vmovdqu (%rdx, %rax, 8), %ymm10;" //former u
	"vmovdqu (%rcx, %rax, 8), %ymm11;" //former t
	"movq $0, %rbx;"
	);

	__asm__ volatile(
"inner:"
	"movq %rax, %rdi;"
	"addq %rbx, %rdi;"
	"push %rax;"
	);

	__asm__ volatile("movq %0, %%rax;":: "a"(a));

	__asm__ volatile(
	// load a[j] ... a[j+7]
	"vmovdqu (%rax, %rbx, 4), %ymm0;"
	"pop %rax;"

	"vpxor %ymm15, %ymm15, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //00000000
	"vpermd %ymm0, %ymm12, %ymm2;" //07060504
	"vpmuludq %ymm3, %ymm2, %ymm3;" //A
	"vpaddq %ymm3, %ymm7, %ymm7;"
	
	"vpermd %ymm1, %ymm14, %ymm3;" //01010101
	"vpmuludq %ymm3, %ymm2, %ymm2;" //B
	"vpaddq %ymm2, %ymm6, %ymm6;"

    "vpsubq %ymm13, %ymm12, %ymm12;" // ymm12(07060504) - ymm13(04040404)
	"vpermd %ymm0, %ymm12, %ymm2;" //03020100
	"vpmuludq %ymm2, %ymm3, %ymm3;" //C
	"vpaddq %ymm3, %ymm10, %ymm10;"

	"vmovdqu %ymm10, (%rdx, %rdi, 8);"
	"vmovdqu 64(%rdx, %rdi, 8), %ymm10;"

	"vpermd %ymm1, %ymm15, %ymm3;" //00000000
	"vpmuludq %ymm2, %ymm3, %ymm3;" //D
	"vpaddq %ymm3, %ymm11, %ymm11;"

	"vmovdqu %ymm11, (%rcx, %rdi, 8);"
	"vmovdqu 64(%rcx, %rdi, 8), %ymm11;"

	"vpermd %ymm1, %ymm13, %ymm3;" //04040404
	"vpmuludq %ymm2, %ymm3, %ymm3;" //E
	"vpaddq %ymm3, %ymm7, %ymm7;"

	"vmovdqu %ymm7, 32(%rcx, %rdi, 8);"
	"vmovdqu 96(%rcx, %rdi, 8), %ymm7;"

    "vpaddq %ymm13, %ymm14, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //05050505
	"vpmuludq %ymm3, %ymm2, %ymm2;" //F
	"vpaddq %ymm2, %ymm6, %ymm6;"

	"vmovdqu %ymm6, 32(%rdx, %rdi, 8);"
	"vmovdqu 96(%rdx, %rdi, 8), %ymm6;"

    "vpaddq %ymm13, %ymm12, %ymm12;"
	"vpermd %ymm0, %ymm12, %ymm2;" //07060504
	"vpmuludq %ymm2, %ymm3, %ymm3;" //G
	"vpaddq %ymm3, %ymm10, %ymm10;"

	"vpermd %ymm1, %ymm13, %ymm3;" //04040404
	"vpmuludq %ymm2, %ymm3, %ymm3;" //H
	"vpaddq %ymm3, %ymm11, %ymm11;"

    "vpaddq %ymm14, %ymm14, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //02020202
	"vpmuludq %ymm3, %ymm2, %ymm3;" //I
	"vpaddq %ymm3, %ymm4, %ymm4;"

    "vpaddq %ymm14, %ymm15, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //03030303
	"vpmuludq %ymm3, %ymm2, %ymm2;" //J
	"vpaddq %ymm2, %ymm5, %ymm5;"

    "vpsubq %ymm13, %ymm12, %ymm12;"
	 //07-21 ymm1 -> ymm0 
	"vpermd %ymm0, %ymm12, %ymm2;" //03020100
	"vpmuludq %ymm2, %ymm3, %ymm3;" //K
	"vpaddq %ymm3, %ymm9, %ymm9;"

	"vmovdqu %ymm9, (%r9, %rdi, 8);"
	"vmovdqu 64(%r9, %rdi, 8), %ymm9;"

    "vpsubq %ymm14, %ymm15, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //02020202
	"vpmuludq %ymm2, %ymm3, %ymm3;" //L
	"vpaddq %ymm3, %ymm8, %ymm8;"

	"vmovdqu %ymm8, (%r8, %rdi, 8);"
	"vmovdqu 64(%r8, %rdi, 8), %ymm8;"

    "vpaddq %ymm13, %ymm15, %ymm15;"
    "vpermd %ymm1, %ymm15, %ymm3;" //06060606
	"vpmuludq %ymm2, %ymm3, %ymm3;" //O
	"vpaddq %ymm3, %ymm4, %ymm4;"

	"vmovdqu %ymm4, 32(%r8, %rdi, 8);"
	"vmovdqu 96(%r8, %rdi, 8), %ymm4;"
    
	"vpaddq %ymm15, %ymm14, %ymm15;"
    "vpermd %ymm1, %ymm15, %ymm3;" //07070707
	"vpmuludq %ymm3, %ymm2, %ymm2;" //P
	"vpaddq %ymm2, %ymm5, %ymm5;"

	"vmovdqu %ymm5, 32(%r9, %rdi, 8);"
	"vmovdqu 96(%r9, %rdi, 8), %ymm5;"

    "vpaddq %ymm13, %ymm12, %ymm12;"
	"vpermd %ymm0, %ymm12, %ymm2;" //07060504
	"vpmuludq %ymm2, %ymm3, %ymm3;" //Q
	"vpaddq %ymm3, %ymm9, %ymm9;"
    
	"vpsubq %ymm14, %ymm15, %ymm15;"
    "vpermd %ymm1, %ymm15, %ymm3;" //06060606
	"vpmuludq %ymm2, %ymm3, %ymm3;" //R
	"vpaddq %ymm3, %ymm8, %ymm8;"


	
	"add $8, %rbx;"

	// DIGITS = (BITSIZE / 29) + 1
	"cmp -52(%rbp), %ebx;"
	"jl inner;"

	// in outer
	"vmovdqu %ymm11, 64(%rcx, %rdi, 8);"
	"vmovdqu %ymm10, 64(%rdx, %rdi, 8);"	
	"vmovdqu %ymm8, 64(%r8, %rdi, 8);"
	"vmovdqu %ymm9, 64(%r9, %rdi, 8);"
	"add $8, %rax;"

	// DIGTTS = (BITSIZE / 29) + 1
	"cmp -52(%rbp), %eax;"
	"jl outer;"
	);
}




void multiply(unsigned int* a, unsigned int* b, unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w,  int digits, int digitstimestwo){

    __asm__ volatile("sub  $64, %rsp;");
    
    __asm__ volatile(
	"movq -8(%rbp), %rax;" //a
	"movq -16(%rbp), %rbx;" //b
	"movq -24(%rbp), %rcx;" //t
	"movq -32(%rbp), %rdx;" //u
	"movq -40(%rbp), %r8;" //v
	"movq -48(%rbp), %r9;" //w
	);
	
        // It is not desirable to use scalar operation,
	// since it can introduce bottleneck in using not only SSE but AVX-512.

	/*
	  Algorithm ( Concept )

	                7 6 5 4 3 2 1 0 ... (i)
	  *             7'6'5'4'3'2'1'0'... (ii)
	  ------------------------------
	                A A A A D D D D     t(rcx){ymm7, ymm11}
	              B B B B C C C C       u(rdx){ymm6, ymm10}
	            I I I I L L L L         v(r8){ymm4, ymm8}
	          J J J J K K K K           w(r9){ymm5, ymm9}
	        H H H H E E E E
	      G G G G F F F F
	    R R R R O O O O
	  Q Q Q Q P P P P
	  ------------------------------

	  Capital alphabet commented out in following code (e.g., //A)
	  corresponds to the alphabet shown in above figure.
	  Steplike number (e.g., 03020100) and repeated number (e.g., 01010101) 
	  indicates multiplicand and multiplier, respectively;
	  that is, the former is (i) and the latter is (ii).
	  

	 */
	

	__asm__ volatile(
	// make ymm12 = 03020100
	"mov $2, %esi;"
	"movd %esi, %xmm14;"
	"inc %esi;"
	"movd %esi, %xmm15;"
	"movlhps %xmm15, %xmm14;"
	"psrlq $1, %xmm15;"
	"pslldq $8, %xmm15;"
	"vinserti128 $1, %xmm14, %ymm15, %ymm12;"

	// make ymm13 = 04040404
    // and make ymm12 = 07060504 (ymm12 + ymm13)
    // and make ymm14 = 01010101
	"mov $4, %esi;"
	"movd %esi, %xmm14;"
	"vpbroadcastq %xmm14, %ymm13;" // ymm13 = 04040404
	"vpaddq %ymm13, %ymm12, %ymm12;" // ymm12 = 07060504
	
    "mov $1, %esi;"
	"movd %esi, %xmm14;"
	"vpbroadcastq %xmm14, %ymm14;" // ymm14 = 01010101
	);


	int i = 0;
	multiply_outer(a, b, t, u, v, w, i);
/*
	__asm__ volatile(

    "movq $0, %rax;"
"outer:"
	"push %rbx;"
	"movq -16(%rbp), %rbx;"
	// load b[i] ... b[i+7]
	"vmovdqu (%rbx, %rax, 4), %ymm1;"
	"pop %rbx;"
	// load t[8-15], u[8-15], u[0-7], t[0-7], respectively
	"vmovdqu 32(%r8, %rax, 8), %ymm4;" //latter v
	"vmovdqu 32(%r9, %rax, 8), %ymm5;" //latter w
	"vmovdqu 32(%rdx, %rax, 8), %ymm6;" //latter u
	"vmovdqu 32(%rcx, %rax, 8), %ymm7;" //latter t
	"vmovdqu (%r8, %rax, 8), %ymm8;" //former v
	"vmovdqu (%r9, %rax, 8), %ymm9;" //former w
	"vmovdqu (%rdx, %rax, 8), %ymm10;" //former u
	"vmovdqu (%rcx, %rax, 8), %ymm11;" //former t
	"movq $0, %rbx;"
"inner:"
	"movq %rax, %rdi;"
	"addq %rbx, %rdi;"
	"push %rax;"
	"movq -8(%rbp), %rax;"
	// load a[j] ... a[j+7]
	"vmovdqu (%rax, %rbx, 4), %ymm0;"
	"pop %rax;"

	"vpxor %ymm15, %ymm15, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //00000000
	"vpermd %ymm0, %ymm12, %ymm2;" //07060504
	"vpmuludq %ymm3, %ymm2, %ymm3;" //A
	"vpaddq %ymm3, %ymm7, %ymm7;"
	
	"vpermd %ymm1, %ymm14, %ymm3;" //01010101
	"vpmuludq %ymm3, %ymm2, %ymm2;" //B
	"vpaddq %ymm2, %ymm6, %ymm6;"

    "vpsubq %ymm13, %ymm12, %ymm12;" // ymm12(07060504) - ymm13(04040404)
	"vpermd %ymm0, %ymm12, %ymm2;" //03020100
	"vpmuludq %ymm2, %ymm3, %ymm3;" //C
	"vpaddq %ymm3, %ymm10, %ymm10;"

	"vmovdqu %ymm10, (%rdx, %rdi, 8);"
	"vmovdqu 64(%rdx, %rdi, 8), %ymm10;"

	"vpermd %ymm1, %ymm15, %ymm3;" //00000000
	"vpmuludq %ymm2, %ymm3, %ymm3;" //D
	"vpaddq %ymm3, %ymm11, %ymm11;"

	"vmovdqu %ymm11, (%rcx, %rdi, 8);"
	"vmovdqu 64(%rcx, %rdi, 8), %ymm11;"

	"vpermd %ymm1, %ymm13, %ymm3;" //04040404
	"vpmuludq %ymm2, %ymm3, %ymm3;" //E
	"vpaddq %ymm3, %ymm7, %ymm7;"

	"vmovdqu %ymm7, 32(%rcx, %rdi, 8);"
	"vmovdqu 96(%rcx, %rdi, 8), %ymm7;"

    "vpaddq %ymm13, %ymm14, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //05050505
	"vpmuludq %ymm3, %ymm2, %ymm2;" //F
	"vpaddq %ymm2, %ymm6, %ymm6;"

	"vmovdqu %ymm6, 32(%rdx, %rdi, 8);"
	"vmovdqu 96(%rdx, %rdi, 8), %ymm6;"

    "vpaddq %ymm13, %ymm12, %ymm12;"
	"vpermd %ymm0, %ymm12, %ymm2;" //07060504
	"vpmuludq %ymm2, %ymm3, %ymm3;" //G
	"vpaddq %ymm3, %ymm10, %ymm10;"

	"vpermd %ymm1, %ymm13, %ymm3;" //04040404
	"vpmuludq %ymm2, %ymm3, %ymm3;" //H
	"vpaddq %ymm3, %ymm11, %ymm11;"

    "vpaddq %ymm14, %ymm14, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //02020202
	"vpmuludq %ymm3, %ymm2, %ymm3;" //I
	"vpaddq %ymm3, %ymm4, %ymm4;"

    "vpaddq %ymm14, %ymm15, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //03030303
	"vpmuludq %ymm3, %ymm2, %ymm2;" //J
	"vpaddq %ymm2, %ymm5, %ymm5;"

    "vpsubq %ymm13, %ymm12, %ymm12;"
	 //07-21 ymm1 -> ymm0 
	"vpermd %ymm0, %ymm12, %ymm2;" //03020100
	"vpmuludq %ymm2, %ymm3, %ymm3;" //K
	"vpaddq %ymm3, %ymm9, %ymm9;"

	"vmovdqu %ymm9, (%r9, %rdi, 8);"
	"vmovdqu 64(%r9, %rdi, 8), %ymm9;"

    "vpsubq %ymm14, %ymm15, %ymm15;"
	"vpermd %ymm1, %ymm15, %ymm3;" //02020202
	"vpmuludq %ymm2, %ymm3, %ymm3;" //L
	"vpaddq %ymm3, %ymm8, %ymm8;"

	"vmovdqu %ymm8, (%r8, %rdi, 8);"
	"vmovdqu 64(%r8, %rdi, 8), %ymm8;"

    "vpaddq %ymm13, %ymm15, %ymm15;"
    "vpermd %ymm1, %ymm15, %ymm3;" //06060606
	"vpmuludq %ymm2, %ymm3, %ymm3;" //O
	"vpaddq %ymm3, %ymm4, %ymm4;"

	"vmovdqu %ymm4, 32(%r8, %rdi, 8);"
	"vmovdqu 96(%r8, %rdi, 8), %ymm4;"
    
	"vpaddq %ymm15, %ymm14, %ymm15;"
    "vpermd %ymm1, %ymm15, %ymm3;" //07070707
	"vpmuludq %ymm3, %ymm2, %ymm2;" //P
	"vpaddq %ymm2, %ymm5, %ymm5;"

	"vmovdqu %ymm5, 32(%r9, %rdi, 8);"
	"vmovdqu 96(%r9, %rdi, 8), %ymm5;"

    "vpaddq %ymm13, %ymm12, %ymm12;"
	"vpermd %ymm0, %ymm12, %ymm2;" //07060504
	"vpmuludq %ymm2, %ymm3, %ymm3;" //Q
	"vpaddq %ymm3, %ymm9, %ymm9;"
    
	"vpsubq %ymm14, %ymm15, %ymm15;"
    "vpermd %ymm1, %ymm15, %ymm3;" //06060606
	"vpmuludq %ymm2, %ymm3, %ymm3;" //R
	"vpaddq %ymm3, %ymm8, %ymm8;"


	
	"add $8, %rbx;"

	// DIGITS = (BITSIZE / 29) + 1
	"cmp -52(%rbp), %ebx;"
	"jl inner;"

	// in outer
	"vmovdqu %ymm11, 64(%rcx, %rdi, 8);"
	"vmovdqu %ymm10, 64(%rdx, %rdi, 8);"	
	"vmovdqu %ymm8, 64(%r8, %rdi, 8);"
	"vmovdqu %ymm9, 64(%r9, %rdi, 8);"
	"add $8, %rax;"

	// DIGTTS = (BITSIZE / 29) + 1
	"cmp -52(%rbp), %eax;"
	"jl outer;"
	);
*/

	//Finally, we do the 29bit carry calculation
	//and add the two result arrays

	calc_carry(t, u, v, w, digitstimestwo);

    
	__asm__ volatile("add $64, %rsp;");
}



void split_init(){
	
	asm(
		"movq $0x1fffffff, %%rax;"
		"movq %%rax, %%xmm1;"
		"movdqa %%xmm1, %%xmm2;"
		"psllq $29, %%xmm2;"
		"movq $0x7ffffffc, %%rax;"
		"movq %%rax, %%xmm3;"
		"pslldq $7, %%xmm3;"
		"movq $0xfffffff80, %%rax;"
		"movq %%rax, %%xmm4;"
		"pslldq $10, %%xmm4;"
		"movq $0x1fffffff0, %%rax;"
		"movq %%rax, %%xmm5;"
		"movq $0x3ffffffe, %%rax;"
		"movq %%rax, %%xmm6;"
		"pslldq $4, %%xmm6;"
		"movq $0x7ffffffc0, %%rax;"
		"movq %%rax, %%xmm7;"
		"pslldq $7, %%xmm7;"
		"movq $0xfffffff8, %%rax;"
		"movq %%rax, %%xmm8;"
		"pslldq $11, %%xmm8;"
	::: "%rax");

}


void split_main(unsigned int* data, unsigned int* result, int digits){

    __asm__ volatile(
	"movq %0, %%rbx;"
	"movq %1, %%rdx;"
	"xor %%rdi, %%rdi;"
	"xor %%rsi, %%rsi;"
	:
	:"r"(data),"r"(result)
	);

	 __asm__ volatile(
    "split:"
	"movdqu (%rbx, %rdi), %xmm0;"
	// [0]
	"movdqa %xmm0, %xmm10;"
	"pand %xmm1, %xmm10;"
	// [1]
	"movdqa %xmm0, %xmm9;"
	"pand %xmm2, %xmm9;"
	"psllq $3, %xmm9;"
	"pxor %xmm9, %xmm10;"
	// [2]
	"movdqa %xmm0, %xmm9;"
	"pand %xmm3, %xmm9;"
	"movq %xmm9, %xmm11;" //lower,  6-bit in 29-bit
	"psrlq $58, %xmm11;"
	"psrldq $8, %xmm9;" //upper, 23-bit in 29-bit
	"psllq $6, %xmm9;"
	"pxor %xmm11, %xmm9;"
	"pslldq $8, %xmm9;"
	"pxor %xmm9, %xmm10;"
	//[3]
	"movdqa %xmm0, %xmm9;"
	"pand %xmm4, %xmm9;"
	"psllq $9, %xmm9;"
	"pxor %xmm9, %xmm10;"

	"movdqa %xmm10, (%rdx, %rsi, 8);"
	"add $14, %rdi;"
	"add $2, %rsi;"
	"cmp -20(%rbp), %esi;"
	"jg _ret;"
	);

	__asm__ volatile(
	"movdqu (%rbx, %rdi), %xmm0;"
	// [0]
	"movdqa %xmm0, %xmm10;"
	"pand %xmm5, %xmm10;"
	"psrlq $4, %xmm10;"
        // [1]
	"movdqa %xmm0, %xmm9;"
	"pand %xmm6, %xmm9;"
	"psrlq $1, %xmm9;"
	"pxor %xmm9, %xmm10;"
	// [2]
	"movdqa %xmm0, %xmm9;"
	"pand %xmm7, %xmm9;"
	"movq %xmm9, %xmm11;" //lower,  2-bit in 29-bit
	"psrlq $62, %xmm11;"
	"psrldq $8, %xmm9;" //upper, 27-bit in 29-bit
	"psllq $2, %xmm9;"
	"pxor %xmm11, %xmm9;"
	"pslldq $8, %xmm9;"
	"pxor %xmm9, %xmm10;"
	//[3]
	"movdqa %xmm0, %xmm9;"
	"pand %xmm8, %xmm9;"
	"psllq $5, %xmm9;"
	"pxor %xmm9, %xmm10;"

	"movdqa %xmm10, (%rdx, %rsi, 8);"
	"add $15, %rdi;"
	"add $2, %rsi;"
	"cmp -20(%rbp), %esi;"
	"jle split;"

	);
	//}
    
	__asm__ volatile("_ret:");
}


void split_29bit(unsigned int* data, unsigned int* result, int digits){

	// initialize
	split_init();

	/*
	 in split,
	 xmm9 : tmp
	 xmm10: result
	 */

	split_main(data, result, digits);
    
}






void combine_29bit(unsigned int* data, unsigned int* result){

    int index, i, j;

    //下記for文で29通りの組み合わせがあるため。29bitという意味ではない
    const int LIMIT = (2*N) / 29;

    //あまりの組み合わせ処理用
    const int REM = 2*N - 29 * LIMIT;
    
    for(index = 0; index < LIMIT; index++){

	j = index * 29; // result
	i = index << 5; // data
	//i = index * 32; // data

	result[j] = (data[i] & 0x1fffffff) ^ ((data[i+1] & 0x7) << 29);
	result[j+1] = ((data[i+1] & 0x1ffffff8) >> 3) ^ ((data[i+2] & 0x3f) << 26);	
	result[j+2] = ((data[i+2] & 0x1fffffc0) >> 6) ^ ((data[i+3] & 0x1ff) << 23);
	result[j+3] = ((data[i+3] & 0x1ffffe00) >> 9) ^ ((data[i+4] & 0xfff) << 20);
	result[j+4] = ((data[i+4] & 0x1ffff000) >> 12) ^ ((data[i+5] & 0x7fff) << 17);
	result[j+5] = ((data[i+5] & 0x1fff8000) >> 15) ^ ((data[i+6] & 0x3ffff) << 14);
	result[j+6] = ((data[i+6] & 0x1ffc0000) >> 18) ^ ((data[i+7] & 0x1fffff) << 11);
	result[j+7] = ((data[i+7] & 0x1fe00000) >> 21) ^ ((data[i+8] & 0xffffff) << 8);
	result[j+8] = ((data[i+8] & 0x1f000000) >> 24) ^ ((data[i+9] & 0x7ffffff) << 5);

	result[j+9] = ((data[i+9] & 0x18000000) >> 27) ^ ((data[i+10] & 0x1fffffff) << 2) ^ ((data[i+11] & 0x1) << 31);

	result[j+10] = ((data[i+11] & 0x1ffffffe) >> 1) ^ ((data[i+12] & 0xf) << 28);
	result[j+11] = ((data[i+12] & 0x1ffffff0) >> 4) ^ ((data[i+13] & 0x7f) << 25);
	result[j+12] = ((data[i+13] & 0x1fffff80) >> 7) ^ ((data[i+14] & 0x3ff) << 22);
	result[j+13] = ((data[i+14] & 0x1ffffc00) >> 10) ^ ((data[i+15] & 0x1fff) << 19);
	result[j+14] = ((data[i+15] & 0x1fffe000) >> 13) ^ ((data[i+16] & 0xffff) << 16);
	result[j+15] = ((data[i+16] & 0x1fff0000) >> 16) ^ ((data[i+17] & 0x7ffff) << 13);
	result[j+16] = ((data[i+17] & 0x1ff80000) >> 19) ^ ((data[i+18] & 0x3fffff) << 10);
	result[j+17] = ((data[i+18] & 0x1fc00000) >> 22) ^ ((data[i+19] & 0x1ffffff) << 7);
	result[j+18] = ((data[i+19] & 0x1e000000) >> 25) ^ ((data[i+20] & 0xfffffff) << 4);

	result[j+19] = ((data[i+20] & 0x10000000) >> 28) ^ ((data[i+21] & 0x1fffffff) << 1) ^ ((data[i+22] & 0x3) << 30);

	result[j+20] = ((data[i+22] & 0x1ffffffc) >> 2) ^ ((data[i+23] & 0x1f) << 27);
	result[j+21] = ((data[i+23] & 0x1fffffe0) >> 5) ^ ((data[i+24] & 0xff) << 24);
	result[j+22] = ((data[i+24] & 0x1fffff00) >> 8) ^ ((data[i+25] & 0x7ff) << 21);
	result[j+23] = ((data[i+25] & 0x1ffff800) >> 11) ^ ((data[i+26] & 0x3fff) << 18);
	result[j+24] = ((data[i+26] & 0x1fffc000) >> 14) ^ ((data[i+27] & 0x1ffff) << 15);
	result[j+25] = ((data[i+27] & 0x1ffe0000) >> 17) ^ ((data[i+28] & 0xfffff) << 12);
	result[j+26] = ((data[i+28] & 0x1ff00000) >> 20) ^ ((data[i+29] & 0x7fffff) << 9);
	result[j+27] = ((data[i+29] & 0x1f800000) >> 23) ^ ((data[i+30] & 0x3ffffff) << 6);
	result[j+28] = ((data[i+30] & 0x1c000000) >> 26) ^ ((data[i+31] & 0x1fffffff) << 3);

	
    }

    j = index * 29;
    i = index << 5;
    int offset;
    
    if(REM == 0){ return; }
    else if(REM < 9){

	for(offset=0; offset<REM; offset++){
	    
	    result[j+offset] =
		((data[i+offset] & ((0x1fffffff >> offset*3) << offset*3)) >> offset*3) ^
		((data[i+offset+1] & (0x1fffffff >> (29 - (offset+1)*3))) << (29 - offset*3));
	}
    }
    else if(REM < 19){

	result[j] = (data[i] & 0x1fffffff) ^ ((data[i+1] & 0x7) << 29);
	result[j+1] = ((data[i+1] & 0x1ffffff8) >> 3) ^ ((data[i+2] & 0x3f) << 26);
	result[j+2] = ((data[i+2] & 0x1fffffc0) >> 6) ^ ((data[i+3] & 0x1ff) << 23);
	result[j+3] = ((data[i+3] & 0x1ffffe00) >> 9) ^ ((data[i+4] & 0xfff) << 20);
	result[j+4] = ((data[i+4] & 0x1ffff000) >> 12) ^ ((data[i+5] & 0x7fff) << 17);
	result[j+5] = ((data[i+5] & 0x1fff8000) >> 15) ^ ((data[i+6] & 0x3ffff) << 14);
	result[j+6] = ((data[i+6] & 0x1ffc0000) >> 18) ^ ((data[i+7] & 0x1fffff) << 11);
	result[j+7] = ((data[i+7] & 0x1fe00000) >> 21) ^ ((data[i+8] & 0xffffff) << 8);
	result[j+8] = ((data[i+8] & 0x1f000000) >> 24) ^ ((data[i+9] & 0x7ffffff) << 5);

	result[j+9] = ((data[i+9] & 0x18000000) >> 27) ^ ((data[i+10] & 0x1fffffff) << 2) ^ ((data[i+11] & 0x1) << 31);

	int shift = 0;
	
	for(offset=10; offset<REM; offset++, shift++){

	    result[j+offset] =
		((data[i+offset+1] & ((0x1fffffff >> (shift*3+1)) << (shift*3+1))) >> (shift*3+1)) ^
		((data[i+offset+2] & (0x1fffffff >> (28 - (shift+1)*3))) << (28 - shift*3));
	    
	}
    }
    else if(REM < 29){

	result[j] = (data[i] & 0x1fffffff) ^ ((data[i+1] & 0x7) << 29);
	result[j+1] = ((data[i+1] & 0x1ffffff8) >> 3) ^ ((data[i+2] & 0x3f) << 26);
	result[j+2] = ((data[i+2] & 0x1fffffc0) >> 6) ^ ((data[i+3] & 0x1ff) << 23);
	result[j+3] = ((data[i+3] & 0x1ffffe00) >> 9) ^ ((data[i+4] & 0xfff) << 20);
	result[j+4] = ((data[i+4] & 0x1ffff000) >> 12) ^ ((data[i+5] & 0x7fff) << 17);
	result[j+5] = ((data[i+5] & 0x1fff8000) >> 15) ^ ((data[i+6] & 0x3ffff) << 14);
	result[j+6] = ((data[i+6] & 0x1ffc0000) >> 18) ^ ((data[i+7] & 0x1fffff) << 11);
	result[j+7] = ((data[i+7] & 0x1fe00000) >> 21) ^ ((data[i+8] & 0xffffff) << 8);
	result[j+8] = ((data[i+8] & 0x1f000000) >> 24) ^ ((data[i+9] & 0x7ffffff) << 5);

	result[j+9] = ((data[i+9] & 0x18000000) >> 27) ^ ((data[i+10] & 0x1fffffff) << 2) ^ ((data[i+11] & 0x1) << 31);

	result[j+10] = ((data[i+11] & 0x1ffffffe) >> 1) ^ ((data[i+12] & 0xf) << 28);
	result[j+11] = ((data[i+12] & 0x1ffffff0) >> 4) ^ ((data[i+13] & 0x7f) << 25);
	result[j+12] = ((data[i+13] & 0x1fffff80) >> 7) ^ ((data[i+14] & 0x3ff) << 22);
	result[j+13] = ((data[i+14] & 0x1ffffc00) >> 10) ^ ((data[i+15] & 0x1fff) << 19);
	result[j+14] = ((data[i+15] & 0x1fffe000) >> 13) ^ ((data[i+16] & 0xffff) << 16);
	result[j+15] = ((data[i+16] & 0x1fff0000) >> 16) ^ ((data[i+17] & 0x7ffff) << 13);
	result[j+16] = ((data[i+17] & 0x1ff80000) >> 19) ^ ((data[i+18] & 0x3fffff) << 10);
	result[j+17] = ((data[i+18] & 0x1fc00000) >> 22) ^ ((data[i+19] & 0x1ffffff) << 7);
	result[j+18] = ((data[i+19] & 0x1e000000) >> 25) ^ ((data[i+20] & 0xfffffff) << 4);

	result[j+19] = ((data[i+20] & 0x10000000) >> 28) ^ ((data[i+21] & 0x1fffffff) << 1) ^ ((data[i+22] & 0x3) << 30);

	int shift = 0;

	for(offset=20; offset<REM; offset++, shift++){
	    
	    result[j+offset] =
		((data[i+offset+2] & ((0x1fffffff >> (shift*3+2)) << (shift*3+2))) >> (shift*3+2)) ^
		((data[i+offset+3] & (0x1fffffff >> (27 - (shift+1)*3))) << (27 - shift*3));

	}
    }
}


#define LoWord(result)  result & 0xffffffff
#define HiWord(result) (result & 0xffffffff00000000) >> 32


void BigMultiply(unsigned int* A, unsigned int* B, unsigned int* T){

    unsigned int i, j, carry;
    unsigned long result;

    for(i=0; i<N; i++){
    	carry = 0;
	    for(j=0; j<N; j++){
	        result = (long)A[j] * B[i] + T[i+j] + carry;
	    
    	    T[i+j] = LoWord(result);
	        carry = HiWord(result);
    	}
	    T[i+j] = carry;
    }
}



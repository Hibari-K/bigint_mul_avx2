#include<stdio.h>
#include<stdlib.h>

#include "ymm_mul.h"


void multiply_outer(unsigned int* a, unsigned int* b, unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int digits);
long multiply_inner(unsigned int* a, unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int digits, long index_outer);
void calc_carry(unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int digitstimestwo);
void cloop(unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int digitstimestwo);




void multiply(unsigned int* a, unsigned int* b, unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w,  int digits, int digitstimestwo){

   /*
	 It is not desirable to use scalar operation,
	 since it can introduce bottleneck in using not only SSE but AVX-512.

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
	
	// ready for multiply_outer
	__asm__ volatile(
	// make ymm12 = 03020100
	"mov $2, %%esi;"
	"movd %%esi, %%xmm14;"
	"inc %%esi;"
	"movd %%esi, %%xmm15;"
	"movlhps %%xmm15, %%xmm14;"
	"psrlq $1, %%xmm15;"
	"pslldq $8, %%xmm15;"
	"vinserti128 $1, %%xmm14, %%ymm15, %%ymm12;"

	// make ymm13 = 04040404
    // and make ymm12 = 07060504 (ymm12 + ymm13)
    // and make ymm14 = 01010101
	"mov $4, %%esi;"
	"movd %%esi, %%xmm14;"
	"vpbroadcastq %%xmm14, %%ymm13;" // ymm13 = 04040404
	"vpaddq %%ymm13, %%ymm12, %%ymm12;" // ymm12 = 07060504
	
    "mov $1, %%esi;"
	"movd %%esi, %%xmm14;"
	"vpbroadcastq %%xmm14, %%ymm14;" // ymm14 = 01010101
	:::"%rsi"
	);


	multiply_outer(a, b, t, u, v, w, digits);

	//Finally, we do the 29bit carry calculation
	//and add the two result arrays

	calc_carry(t, u, v, w, digitstimestwo);
    
}

void multiply_outer(unsigned int* a, unsigned int* b, unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int digits){


	long index_outer, index_inner, rdi;

	for(index_outer = 0; index_outer < digits; index_outer += 8){

		__asm__ volatile(
			// load b[i] ... b[i+7]
			"vmovdqu (%0, %5, 4), %%ymm1;"
	
			// load t[8-15], u[8-15], u[0-7], t[0-7], respectively
			"vmovdqu 32(%3, %5, 8), %%ymm4;" //latter v
			"vmovdqu 32(%4, %5, 8), %%ymm5;" //latter w
			"vmovdqu 32(%2, %5, 8), %%ymm6;" //latter u
			"vmovdqu 32(%1, %5, 8), %%ymm7;" //latter t
			"vmovdqu (%3, %5, 8), %%ymm8;" //former v
			"vmovdqu (%4, %5, 8), %%ymm9;" //former w
			"vmovdqu (%2, %5, 8), %%ymm10;" //former u
			"vmovdqu (%1, %5, 8), %%ymm11;" //former t
			::"r"(b), "r"(t), "r"(u), "r"(v), "r"(w), "r"(index_outer)
		);


		index_inner = multiply_inner(a, t, u, v, w, digits, index_outer);

		rdi = index_inner + index_outer;

		__asm__ volatile(
			"vmovdqu %%ymm11, 64(%0, %4, 8);"
			"vmovdqu %%ymm10, 64(%1, %4, 8);"	
			"vmovdqu %%ymm8, 64(%2, %4, 8);"
			"vmovdqu %%ymm9, 64(%3, %4, 8);"
			::"r"(t), "r"(u), "r"(v), "r"(w), "r"(rdi)
		);
	}
}



long multiply_inner(unsigned int* a, unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int digits, long index_outer){


	long index_inner;
	for(index_inner = 0; index_inner < digits; index_inner += 8){

	long index = index_inner + index_outer;

	__asm__ volatile(
	// load a[j] ... a[j+7]
	"vmovdqu (%0, %5, 4), %%ymm0;"

	// calculate multiply
	"vpxor %%ymm15, %%ymm15, %%ymm15;"
	"vpermd %%ymm1, %%ymm15, %%ymm3;" //00000000
	"vpermd %%ymm0, %%ymm12, %%ymm2;" //07060504
	"vpmuludq %%ymm3, %%ymm2, %%ymm3;" //A
	"vpaddq %%ymm3, %%ymm7, %%ymm7;"
	
	"vpermd %%ymm1, %%ymm14, %%ymm3;" //01010101
	"vpmuludq %%ymm3, %%ymm2, %%ymm2;" //B
	"vpaddq %%ymm2, %%ymm6, %%ymm6;"

    "vpsubq %%ymm13, %%ymm12, %%ymm12;" // ymm12(07060504) - ymm13(04040404)
	"vpermd %%ymm0, %%ymm12, %%ymm2;" //03020100
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //C
	"vpaddq %%ymm3, %%ymm10, %%ymm10;"

	"vmovdqu %%ymm10, (%2, %6, 8);"
	"vmovdqu 64(%2, %6, 8), %%ymm10;"

	"vpermd %%ymm1, %%ymm15, %%ymm3;" //00000000
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //D
	"vpaddq %%ymm3, %%ymm11, %%ymm11;"

	"vmovdqu %%ymm11, (%1, %6, 8);"
	"vmovdqu 64(%1, %6, 8), %%ymm11;"

	"vpermd %%ymm1, %%ymm13, %%ymm3;" //04040404
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //E
	"vpaddq %%ymm3, %%ymm7, %%ymm7;"

	"vmovdqu %%ymm7, 32(%1, %6, 8);"
	"vmovdqu 96(%1, %6, 8), %%ymm7;"

    "vpaddq %%ymm13, %%ymm14, %%ymm15;"
	"vpermd %%ymm1, %%ymm15, %%ymm3;" //05050505
	"vpmuludq %%ymm3, %%ymm2, %%ymm2;" //F
	"vpaddq %%ymm2, %%ymm6, %%ymm6;"

	"vmovdqu %%ymm6, 32(%2, %6, 8);"
	"vmovdqu 96(%2, %6, 8), %%ymm6;"

    "vpaddq %%ymm13, %%ymm12, %%ymm12;"
	"vpermd %%ymm0, %%ymm12, %%ymm2;" //07060504
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //G
	"vpaddq %%ymm3, %%ymm10, %%ymm10;"

	"vpermd %%ymm1, %%ymm13, %%ymm3;" //04040404
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //H
	"vpaddq %%ymm3, %%ymm11, %%ymm11;"

    "vpaddq %%ymm14, %%ymm14, %%ymm15;"
	"vpermd %%ymm1, %%ymm15, %%ymm3;" //02020202
	"vpmuludq %%ymm3, %%ymm2, %%ymm3;" //I
	"vpaddq %%ymm3, %%ymm4, %%ymm4;"

    "vpaddq %%ymm14, %%ymm15, %%ymm15;"
	"vpermd %%ymm1, %%ymm15, %%ymm3;" //03030303
	"vpmuludq %%ymm3, %%ymm2, %%ymm2;" //J
	"vpaddq %%ymm2, %%ymm5, %%ymm5;"

    "vpsubq %%ymm13, %%ymm12, %%ymm12;"
	"vpermd %%ymm0, %%ymm12, %%ymm2;" //03020100
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //K
	"vpaddq %%ymm3, %%ymm9, %%ymm9;"

	"vmovdqu %%ymm9, (%4, %6, 8);"
	"vmovdqu 64(%4, %6, 8), %%ymm9;"

    "vpsubq %%ymm14, %%ymm15, %%ymm15;"
	"vpermd %%ymm1, %%ymm15, %%ymm3;" //02020202
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //L
	"vpaddq %%ymm3, %%ymm8, %%ymm8;"

	"vmovdqu %%ymm8, (%3, %6, 8);"
	"vmovdqu 64(%3, %6, 8), %%ymm8;"

    "vpaddq %%ymm13, %%ymm15, %%ymm15;"
    "vpermd %%ymm1, %%ymm15, %%ymm3;" //06060606
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //O
	"vpaddq %%ymm3, %%ymm4, %%ymm4;"

	"vmovdqu %%ymm4, 32(%3, %6, 8);"
	"vmovdqu 96(%3, %6, 8), %%ymm4;"
    
	"vpaddq %%ymm15, %%ymm14, %%ymm15;"
    "vpermd %%ymm1, %%ymm15, %%ymm3;" //07070707
	"vpmuludq %%ymm3, %%ymm2, %%ymm2;" //P
	"vpaddq %%ymm2, %%ymm5, %%ymm5;"

	"vmovdqu %%ymm5, 32(%4, %6, 8);"
	"vmovdqu 96(%4, %6, 8), %%ymm5;"

    "vpaddq %%ymm13, %%ymm12, %%ymm12;"
	"vpermd %%ymm0, %%ymm12, %%ymm2;" //07060504
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //Q
	"vpaddq %%ymm3, %%ymm9, %%ymm9;"
    
	"vpsubq %%ymm14, %%ymm15, %%ymm15;"
    "vpermd %%ymm1, %%ymm15, %%ymm3;" //06060606
	"vpmuludq %%ymm2, %%ymm3, %%ymm3;" //R
	"vpaddq %%ymm3, %%ymm8, %%ymm8;"

	::"r"(a), "r"(t), "r"(u), "r"(v), "r"(w), "r"(index_inner), "r"(index)
	);
	}

	// subtract extra 8 which is added in "for" expression. 
	return index_inner - 8;
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
	 */

	 __asm__ volatile(
		"pxor %%mm6, %%mm6;"
		"mov $0x1fffffff, %%eax;"
		"movd %%eax, %%mm7;"
		//[0]
		"movq (%0), %%mm0;"
		"movq %%mm0, %%mm6;"
		"pand %%mm7, %%mm0;"
		"movd %%mm0, (%0);"
		"psrlq $29, %%mm6;"
		
		"mov $1, %%rbx;"
		//[1]
		"movq (%0, %%rbx, 8), %%mm0;"
		"movq -8(%1, %%rbx, 8), %%mm1;"
		"paddq %%mm1, %%mm0;"
		"paddq %%mm6, %%mm0;"
		"movq %%mm0, %%mm6;"
		"pand %%mm7, %%mm6;"
		"movd %%mm6, (%0, %%rbx, 4);"
		"psrlq $29, %%mm0;"
	
		"inc %%ebx;"
		//[2]
		"movq (%0, %%rbx, 8), %%mm6;"
		"movq -8(%1, %%rbx, 8), %%mm1;"
		"movq -16(%2, %%rbx, 8), %%mm2;"
		"paddq %%mm1, %%mm6;"
		"paddq %%mm2, %%mm6;"
		"paddq %%mm0, %%mm6;"
		"movq %%mm6, %%mm0;"
		"pand %%mm7, %%mm0;"
		"movd %%mm0, (%0, %%rbx, 4);"
		"psrlq $29, %%mm6;"

		::"r"(t), "r"(u), "r"(v)
		:"%rax", "%rbx"
	);

	
	cloop(t, u, v, w, digitstimestwo);
	

	long rbx = digitstimestwo - 1;

	__asm__ volatile(
		"movd %%mm6, (%0, %1, 4);"
		::"r"(t), "r"(rbx)
	);
}



void cloop(unsigned int* t, unsigned int* u, unsigned int* v, unsigned int* w, int digitstimestwo){

	long i;
	for(i=3; i<=digitstimestwo; i+=2){

		__asm__ volatile(
			"movq (%0, %4, 8), %%mm0;"
			"movq -8(%1, %4, 8), %%mm1;"
			"movq -16(%2, %4, 8), %%mm2;"
			"movq -24(%3, %4, 8), %%mm3;"
			"paddq %%mm1, %%mm0;"
			"paddq %%mm2, %%mm0;"
			"paddq %%mm3, %%mm0;"
			"paddq %%mm6, %%mm0;"
			"movq %%mm0, %%mm6;"
			"pand %%mm7, %%mm6;"
			"movd %%mm6, (%0, %4, 4);"
			"psrlq $29, %%mm0;"


			"movq 8(%0, %4, 8), %%mm6;"
			"movq (%1, %4, 8), %%mm1;"
			"movq -8(%2, %4, 8), %%mm2;"
			"movq -16(%3, %4, 8), %%mm3;"
			"paddq %%mm1, %%mm6;"
			"paddq %%mm2, %%mm6;"
			"paddq %%mm3, %%mm6;"
			"paddq %%mm0, %%mm6;"
			"movq %%mm6, %%mm0;"
			"pand %%mm7, %%mm0;"
			"movd %%mm0, 4(%0, %4, 4);"
			"psrlq $29, %%mm6;"

			::"r"(t), "r"(u), "r"(v), "r"(w), "r"(i)
		);
	}
}

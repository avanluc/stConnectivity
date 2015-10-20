#pragma once

/* DEBUG CONFIG */
#define				     BFS	0
#define				   DEBUG	1
#define				  ATOMIC	1
#define 		   BOTTOM_UP	1
#define				  N_TEST 	100

/* CUDA CONFIG */
#define					 Tid 	threadIdx.x
#define					 Bid 	blockIdx.x
#define			  BLOCK_SIZE	256
#define		   Thread_Per_SM	2048
#define				N_OF_SMs	12
#define	   	  	 SMem_Per_SM	49152
#define 			Int_Size	4
#define		  MAX_CONCURR_TH	(Thread_Per_SM * N_OF_SMs)
#define		 SMem_Per_Thread	(SMem_Per_SM / Thread_Per_SM)
#define					GTid	(Bid * BLOCK_SIZE) + Tid

const int IntSMem_Per_Thread  =  SMem_Per_Thread / Int_Size;
const int REG_QUEUE  = 	32;

#define                  	MIN_V(a, b)		((a) > (b) ? (b) : (a))
#define                  	MAX_V(a, b)		((a) < (b) ? (b) : (a))
#define    	   MAX_CONCURR_BL(BlockDim)		( MAX_CONCURR_TH / (BlockDim) )
#define		   SMem_Per_Block(BlockDim)		( MIN_V( SMem_Per_Thread * (BlockDim) , 49152) )
#define     IntSMem_Per_Block(BlockDim)		( MIN_V( IntSMem_Per_Thread * (BlockDim) , 49152 / Int_Size) )
#define     SMem_Per_BlockNum(BlockNum)		( MAX_V( SMem_Per_SM / (BlockNum) , SMem_Per_SM / (2 * Int_Size) ) )
#define  IntSMem_Per_BlockNum(BlockNum)		( SMem_Per_BlockNum(BlockNum) / Int_Size )

#define				TEMP_POS 	0
#define			END_TEMP_POS 	(TEMP_POS + (35 * Int_Size))
#define       	  F2Size_POS	END_TEMP_POS
//#define			   F1_OFFSET	(F2Size_POS + (2 * Int_Size))
//#define			   EXIT_FLAG	(F2Size_POS + Int_Size)
#define			   F1_OFFSET	(F2Size_POS + Int_Size)
#define		   FRONTIER_SIZE 	((SMem_Per_Block(BLOCK_SIZE) - F1_OFFSET))
#define			   F2_OFFSET	(F1_OFFSET + FRONTIER_SIZE)
#define	BLOCK_FRONTIER_LIMIT 	(FRONTIER_SIZE / Int_Size)

#define cudaAssert(condition, pos) \
  if (!(condition)){ printf("Assertion %s failed!\tpos = %d\n", #condition, pos); asm("trap;"); }

const int SOURCES[] = {10, 50, 100, 500, 1000, 2000, 4000, 6000, 8000, 10000};
const int LENGTH = sizeof(SOURCES) / sizeof(int);

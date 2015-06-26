#pragma once

/* DEBUG CONFIG */
#define				   DEBUG	1
#define				   ATOMIC	1
#define			SINGLE_BLOCK	0
#define				  N_TEST 	1000

/* CUDA CONFIG */
#define					 Tid 	threadIdx.x
#define			  BLOCK_SIZE	512
#define		   Thread_Per_SM	2048
#define				N_OF_SMs	12
#define	   	  	 SMem_Per_SM	49152
#define		  MAX_CONCURR_TH	(Thread_Per_SM * N_OF_SMs)
#define		 SMem_Per_Thread	(SMem_Per_SM / Thread_Per_SM)

const int IntSMem_Per_Thread  =  SMem_Per_Thread / 4;
const int REG_QUEUE  = 	32;


#define                  	MIN_V(a, b)		((a) > (b) ? (b) : (a))
#define                  	MAX_V(a, b)		((a) < (b) ? (b) : (a))
#define    	   MAX_CONCURR_BL(BlockDim)		( MAX_CONCURR_TH / (BlockDim) )
#define		   SMem_Per_Block(BlockDim)		( MIN_V( SMem_Per_Thread * (BlockDim) , 49152) )
#define     IntSMem_Per_Block(BlockDim)		( MIN_V( IntSMem_Per_Thread * (BlockDim) , 49152 / 4) )
#define     SMem_Per_BlockNum(BlockNum)		( MAX_V( SMem_Per_SM / (BlockNum) , SMem_Per_SM / 8 ) )
#define  IntSMem_Per_BlockNum(BlockNum)		( SMem_Per_BlockNum(BlockNum) / 4 )

#define				TEMP_POS 	0
#define			END_TEMP_POS 	(TEMP_POS + (34 * 4))
#define       	  F2Size_POS	END_TEMP_POS
#define			   F1_OFFSET	(F2Size_POS + 8)
#define		   FRONTIER_SIZE 	((SMem_Per_Block(BLOCK_SIZE) - F1_OFFSET))
#define			   F2_OFFSET	(F1_OFFSET + FRONTIER_SIZE)
#define	BLOCK_FRONTIER_LIMIT 	(FRONTIER_SIZE / 4)

#define cudaAssert(condition, pos) \
  if (!(condition)){ printf("Assertion %s failed!\tpos = %d\n ", #condition, pos); asm("trap;"); }

// #define				TEMP_POS 	0
// #define			END_TEMP_POS 	(TEMP_POS + (34 * 4))
// #define       	  F2Size_POS	END_TEMP_POS
// #define			   F1_OFFSET	(F2Size_POS + 8)
// #define		   FRONTIER_SIZE 	((SMem_Per_SM - F1_OFFSET))
// #define			   F2_OFFSET	(FRONTIER_SIZE + F1_OFFSET)
// #define	BLOCK_FRONTIER_LIMIT 	(FRONTIER_SIZE / 4)


// #define	   SM_BYTE_PER_BLOCK	49152
// #define       	  F2Size_POS	(SM_BYTE_PER_BLOCK - (sizeof(int) * 2))
// #define	BLOCK_FRONTIER_LIMIT 	(F2Size_POS / (2 * sizeof(int)))
// #define			   F1_OFFSET	0
// #define			   F2_OFFSET	(F2Size_POS / 2)

/*
#define HASHTABLE_BLOCK_POS  0
#define END_OF_HASHTABLE	(4096 * 8)	// 8: long long int size
#define       F2Size_POS	END_OF_HASHTABLE
#define         TEMP_POS	(F2Size_POS + 4)
#define     END_TEMP_POS	(TEMP_POS + 32 * 4)
#define    FRONTIER_SIZE	(((49152 - END_TEMP_POS) / 2) - 2)		//-2 align
#define     F1_BLOCK_POS	(END_TEMP_POS)
#define     F2_BLOCK_POS	(F1_BLOCK_POS + FRONTIER_SIZE)
*/

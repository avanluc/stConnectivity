#pragma once

#define				   DEBUG	1
#define				  N_TEST 	50

#define					 Tid 	threadIdx.x
#define			  BLOCK_SIZE	1024
#define		   Thread_Per_SM	2048
#define				N_OF_SMs	12
#define		  MAX_CONCURR_TH	(Thread_Per_SM * N_OF_SMs)
#define	MAX_CONCURR_BL(BlockDim)	( MAX_CONCURR_TH / (BlockDim) )

#define	   SM_BYTE_PER_BLOCK	49152
#define				TEMP_POS 	0
#define			END_TEMP_POS 	(TEMP_POS + (34 * 4))
#define       	  F2Size_POS	END_TEMP_POS
#define			   F1_OFFSET	(F2Size_POS + 8)
#define			   F2_OFFSET	((SM_BYTE_PER_BLOCK - F1_OFFSET) / 2)
#define	BLOCK_FRONTIER_LIMIT 	((F2_OFFSET - F1_OFFSET) / 4)


const int REG_QUEUE  = 	32;


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

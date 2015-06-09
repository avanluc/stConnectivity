#pragma once

#define				   DEBUG	1
#define				  N_TEST 	1

#define					 Tid 	threadIdx.x
#define   TOTAL_SM_PER_BLOCK	49152
#define	BLOCK_FRONTIER_LIMIT 	TOTAL_SM_PER_BLOCK / (2 * sizeof(int))
#define			SMEMORY_SIZE	1024
#define			  BLOCK_SIZE	1024
#define		   Thread_Per_SM	2048
#define				N_OF_SMs	12
#define		  MAX_CONCURR_TH	(Thread_Per_SM * N_OF_SMs)
#define	MAX_CONCURR_BL(BlockDim)	( MAX_CONCURR_TH / (BlockDim) )

#define	   SM_BYTE_PER_BLOCK	49152
#define			   F1_OFFSET	0
#define			   F2_OFFSET	SM_BYTE_PER_BLOCK/2

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

const int REG_QUEUE  = 	32;
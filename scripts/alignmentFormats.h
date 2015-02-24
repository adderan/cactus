#ifndef PAIRWISEALIGNMENT_H
#define PAIRWISEALIGNMENT_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sonLib.h"


//MAF alignment line for one sequence
typedef struct _mafAlignment 
{
	char *name;
	int start;
	int size;
	char strand;
	int length;
	char *seq;
} MafAlignment;

//MAF alignment of two sequences.
typedef struct _pairwiseMafBlock 
{
	MafAlignment *a1;
	MafAlignment *a2;
	int score;
} PairwiseMafBlock;

//cigar representation of one pairwise alignment
typedef struct _cigarAlignment 
{
	char *queryname, *targetname;
	int querystart, querystop, targetstart, targetstop;
	int score;
	char querystrand, targetstrand;
	char *cigarstring;
} CigarAlignment;
void printMafBlock(PairwiseMafBlock *block);
CigarAlignment *mafToCigar(PairwiseMafBlock *block);
void printCigar(CigarAlignment *cigar);
void cigarAlignmentDestruct(CigarAlignment *alignment);
PairwiseMafBlock *parseMaf(FILE *f);
void pairwiseMafBlockDestruct(PairwiseMafBlock *block);

#endif

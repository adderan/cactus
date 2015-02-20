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

MafAlignment *mafAlignmentConstruct(char *line) 
{
	int seqlen = strlen(line);
	char strand;
	char firstchar;
	char *name = st_calloc(1, seqlen + 1);
	char *seq = st_calloc(1, seqlen + 1);
	int start, size, length;
	sscanf(line, "%c %s %d %d %c %d %s", &firstchar, name, &start, &size, &strand, &length, seq);
	MafAlignment *alignment = st_malloc(sizeof(MafAlignment));
	alignment->name = name;
	alignment->start = start;
	alignment->size = size;
	alignment->strand = strand;
	//strcpy(alignment->strand, strand);
	alignment->length = length;
	alignment->seq = seq;
	return(alignment);
}

PairwiseMafBlock *pairwiseMafBlockConstruct(char *scoreline, char *line1, char *line2) 
{
	char *space = " ";
	char *eq = "=";

	//parse score line
	strtok(scoreline, space);
	char *scoreEquation = strtok(NULL, space);
	char *scoreword = strtok(scoreEquation, eq);
	assert(stString_eq(scoreword, "score"));
	char *scorestr = strtok(NULL, eq);
	int score = atoi(scorestr);

	//parse alignment lines
	MafAlignment *alignment1 = mafAlignmentConstruct(line1);
	MafAlignment *alignment2 = mafAlignmentConstruct(line2);
	PairwiseMafBlock *block = st_malloc(sizeof(PairwiseMafBlock));
	block->a1 = alignment1;
	block->a2 = alignment2;
	block->score = score;
	return(block);
}


void pairwiseMafBlockDestruct(PairwiseMafBlock *block) 
{
	free(block->a1->name);
	free(block->a2->name);
	free(block->a1->seq);
	free(block->a2->seq);
	free(block->a1);
	free(block->a2);
	free(block);
}
void printMafBlock(PairwiseMafBlock *block) 
{
	printf("a score=%d\n", block->score);
	printf("s %s %d %d %c %d %s\n", block->a1->name, block->a1->start, block->a1->size, block->a1->strand, block->a1->length, block->a1->seq);
	printf("s %s %d %d %c %d %s\n", block->a2->name, block->a2->start, block->a2->size, block->a2->strand, block->a2->length, block->a2->seq);
}

char getCigarMode(char c1, char c2) 
{
	if(c1 == '-') {
		return('I');
	}
	if(c2 == '-') {
		return('D');
	}
	else {
		return('M');
	}
}
		
char *makeCigarString(char *target, char *query) 
{
	int len = strlen(target);
	char mode = getCigarMode(target[0], query[0]);
	int counter = 1;
	int i;
	char *cigarstring = st_calloc(len * 4, 1);
	for(i = 1; i < len; i++) {
		char newmode = getCigarMode(target[i], query[i]);
		if(newmode != mode) {
			char *newstr = stString_print("%c %d ", mode, counter);
			strcat(cigarstring, newstr);
			free(newstr);
			counter = 0;
			mode = newmode;
		}
		counter++;
	}
	char *newstr = stString_print("%c %d", mode, counter);
	strcat(cigarstring, newstr);
	free(newstr);
	return(cigarstring);
}

CigarAlignment *mafToCigar(PairwiseMafBlock *block)
{
	CigarAlignment *cigarAlignment = st_malloc(sizeof(CigarAlignment));
	cigarAlignment->targetname = stString_copy(block->a1->name);
	cigarAlignment->queryname = stString_copy(block->a2->name);

	cigarAlignment->targetstart = block->a1->start;
	cigarAlignment->querystart = block->a2->start;
	cigarAlignment->targetstop = cigarAlignment->targetstart + block->a1->size;
	cigarAlignment->querystop = cigarAlignment->querystart + block->a2->size;

	//wrap around if strand is negative
	if(block->a1->strand ==  '-') {
		cigarAlignment->targetstart = block->a1->length - cigarAlignment->targetstart;
		cigarAlignment->targetstop = block->a1->length - cigarAlignment->targetstop;
	}
	if(block->a2->strand == '-') {
		cigarAlignment->querystart = block->a2->length - cigarAlignment->querystart;
		cigarAlignment->querystop = block->a2->length - cigarAlignment->querystop;
	}
	
	cigarAlignment->targetstrand = block->a1->strand;
	cigarAlignment->querystrand = block->a2->strand;
	cigarAlignment->score = block->score;
	cigarAlignment->cigarstring = makeCigarString(block->a1->seq, block->a2->seq);
	return(cigarAlignment);
}
void printCigar(CigarAlignment *cigar)
{
	printf("cigar: %s %d %d %c %s %d %d %c %d %s\n", cigar->queryname, cigar->querystart, cigar->querystop, cigar->querystrand, cigar->targetname, cigar->targetstart, cigar->targetstop, cigar->targetstrand, cigar->score, cigar->cigarstring);
}

void cigarAlignmentDestruct(CigarAlignment *alignment) 
{
	free(alignment->targetname);
	free(alignment->queryname);
	free(alignment->cigarstring);
	free(alignment);
}

int main()
{
	char *line = NULL;
	size_t size;
	int len;

	char *scoreline = NULL;
	char *line1 = NULL;
	char *line2 = NULL;
	while((len =  (getline(&line, &size, stdin))) != -1) {
		if(len == 0) continue;
		if(line[0] == '\n') continue;
		if(line[0] == '#') continue;

		if (!scoreline && strlen(line) != 0 && line[0] == 'a') {
			scoreline = st_calloc(1, len + 1);
			strcpy(scoreline, line);
		}
		else if (scoreline && !line1 && strlen(line) != 0 && line[0] == 's') {
			line1 = st_calloc(1, len + 1);	
			strcpy(line1, line);
		}
		else if (scoreline && line1 && !line2 && strlen(line) != 0 && line[0] == 's') {
			line2 = st_calloc(1, len + 1);
			strcpy(line2, line);
		}

		if (scoreline && line1 && line2) {
			PairwiseMafBlock *block = pairwiseMafBlockConstruct(scoreline, line1, line2);
			CigarAlignment *cigar = mafToCigar(block);
			printCigar(cigar);

			pairwiseMafBlockDestruct(block);
			cigarAlignmentDestruct(cigar);
			
			free(scoreline);
			free(line1);
			free(line2);
			scoreline = NULL;
			line1 = NULL;
			line2 = NULL;
		}
	
	}
	return(0);
}


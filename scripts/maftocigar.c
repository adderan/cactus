#include "alignmentFormats.h"
#include <stdio.h>
#include <stdlib.h>


int isTrivial(CigarAlignment *c)
{
	if(c->querystart == c->targetstart && c->querystop == c->targetstop) return 1;
	else return 0;
}
int main(int argc, char **argv)
{
	int notrivial = 0;
	PairwiseMafBlock *block;
	if(argc > 1 && strcmp(argv[1], "--notrivial") == 0) {
		notrivial = 1;
		fprintf(stderr, "using --notrivial.\n");
	}
	while((block = parseMaf(stdin)) != NULL) {
		CigarAlignment *ca = mafToCigar(block);
		if(notrivial && isTrivial(ca)) continue;
		printCigar(ca);
		cigarAlignmentDestruct(ca);
		pairwiseMafBlockDestruct(block);
	}
	return(0);
}

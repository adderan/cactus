#include "alignmentFormats.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	//PairwiseMafBlock *block;
	//block = parseMaf(stdin);
	//pairwiseMafBlockDestruct(block);
	PairwiseMafBlock *block;
	while((block = parseMaf(stdin)) != NULL) {
		CigarAlignment *ca = mafToCigar(block);
		printCigar(ca);
		cigarAlignmentDestruct(ca);
		pairwiseMafBlockDestruct(block);
	}
	return(0);
}

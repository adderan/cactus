#include "alignmentFormats.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	PairwiseMafBlock *block;
	while((block = parseMaf(stdin)) != NULL) {
		CigarAlignment *ca = mafToCigar(block);
		printCigar(ca);
		cigarAlignmentDestruct(ca);
		pairwiseMafBlockDestruct(block);
	}
	return(0);
}

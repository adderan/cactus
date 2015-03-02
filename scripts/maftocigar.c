#include "alignmentFormats.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>


int isTrivial(CigarAlignment *c)
{
	if(c->querystart == c->targetstart && c->querystop == c->targetstop) return 1;
	else return 0;
}
int main(int argc, char **argv)
{
	static int notrivial;
	PairwiseMafBlock *block;
	int c;
	while(1) {
		static struct option opts[] = 
			{{"notrivial", no_argument, &notrivial, 1}};
		int option_index = 0;
		c = getopt_long(argc, argv, "", opts, &option_index);
		if(c == -1) break;
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

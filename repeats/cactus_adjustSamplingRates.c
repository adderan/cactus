#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>

#include "commonC.h"
#include "hashTableC.h"
#include "bioioC.h"
#include "cactus.h"
#include "avl.h"
#include "pairwiseAlignment.h"
#include "blastAlignmentLib.h"


int main(int argc, char **argv) {
	struct option opts[] = { {"initialSamplingRates", required_argument, NULL, '1'},
                             {"alignments", required_argument, NULL, '2'},
							 {"seed", required_argument, NULL, '3'},
							 {"seq1", required_argument, NULL, '4'},
							 {"seq2", required_argument, NULL, '5'},
                             {0, 0, 0, 0} };
	char *initialSamplingRatesFilename;
	char *alignmentsFilename;
	char *seed;
	char *seq1Filename;
	char *seq2Filename;
    int flag;
    while((flag = getopt_long(argc, argv, "", opts, NULL)) != -1) {
        switch(flag) {
        case '1':
            initialSamplingRatesFilename = stString_copy(optarg);
            break;
        case '2':
            alignmentsFilename = stString_copy(optarg);
            break;
		case '3':
			seed = stString_copy(optarg);
			break;
        }
		
    }
	FILE *alignmentsFile = fopen(alignmentsFilename, "r");
	FILE *initialSamplingRatesFile = fopen(alignmentsFilename, "r");
	FILE *seq1File = fopen(seq1Filename, "r");
	FILE *seq2File = fopen(seq2Filename, "r");

    stList *seq1 = stList_construct();
	stList *seqLengths1 = stList_construct();
	stList *seq1Names = stList_construct();
	stList *seq2 = stList_construct();
	stList *seqLengths2 = stList_construct();
	stList *seq2Names = stList_construct();
	//fastaRead(seq1File, seq1, seqLengths1, seq1Names);
	//fastaRead(seq2File, seq2, seqLengths2, seq2Names);

	struct PairwiseAlignment *pairwiseAlignment;
	while((pairwiseAlignment = cigarRead(alignmentsFile)) != NULL) {
		fprintf(stderr, "%s\n", pairwiseAlignment->contig2);
	}
}

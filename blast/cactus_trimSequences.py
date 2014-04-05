#!/usr/bin/env python
from argparse import ArgumentParser
from collections import defaultdict
from operator import itemgetter

def windowFilter(windowSize, threshold, blockDict, seqLengths):
    ret = defaultdict(list)
    for seq, blocks in blockDict.items():
        curBlock = 0
        inRegion = False
        regionStart = 0
        for i in xrange(seqLengths[seq]):
            score = 0
            while curBlock < len(blocks) and blocks[curBlock][1] < i:
                curBlock += 1
            for block in blocks[curBlock:]:
                if block[0] > i + windowSize:
                    break
                size = min(block[1], i + windowSize) - max(i, block[0])
                score += block[2]*size
            score /= windowSize
            if score >= threshold and not inRegion:
                regionStart = i
                inRegion = True
            elif score < threshold and inRegion:
                ret[seq].append((regionStart, i + windowSize))
                inRegion = False
    return ret

def uniquifyBlocks(blocksDict, mergeDistance):
    """Take list of blocks and return sorted list of non-overlapping and
    blocks (merging blocks that are mergeDistance or less apart)."""
    ret = {}
    for chr, blocks in blocksDict.items():
        blocks = sorted(blocks, key=itemgetter(0))
        newBlocks = []
        prevBlock = None
        for block in blocks:
            if prevBlock is None:
                prevBlock = block
            else:
                if prevBlock[1] < block[0] - mergeDistance:
                    newBlocks.append(prevBlock)
                    prevBlock = block
                else:
                    prevBlock = (prevBlock[0], block[1])
        if prevBlock is not None:
            newBlocks.append(prevBlock)
        ret[chr] = newBlocks
    return ret

def getSeparateBedBlocks(bedFile):
    """Get dict of sequence -> (start, stop, score) regions from bed file,
    counting BED12 exons separately."""
    ret = defaultdict(list)
    for line in bedFile:
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        fields = line.split('\t')
        chr = fields[0]
        start = int(fields[1])
        stop = int(fields[2])
        score = int(fields[4])
        if len(fields) <= 9:
            ret[chr].append((start, stop, score))
        else:
            assert(len(fields) == 12)
            blockSizes = map(int, filter(lambda x: x != '',
                                         fields[10].split(',')))
            blockStarts = map(int, filter(lambda x: x != '',
                                          fields[11].split(',')))
            for blockStart, blockSize in zip(blockSizes, blockStarts):
                nonRelativeBlockStart = start + blockStart
                nonRelativeBlockEnd = nonRelBlockStart + blockSize
                ret[chr].append((nonRelativeBlockStart, nonRelativeBlockEnd,
                                 score))
    return ret

def getSeqLengths(fastaFile):
    ret = defaultdict(int)
    chr = ""
    for line in fastaFile:
        line = line.strip()
        if line[0] == '>':
            chr = line.split()[0][1:]
        else:
            ret[chr] += len(line)
    return ret

def complementBlocks(blocksDict, seqLengths):
    """Complement a sorted block-dict."""
    ret = defaultdict(list)
    for chr, blocks in blocksDict.items():
        len = seqLengths[chr]
        start = 0
        for block in blocks:
            if block[0] == start:
                start = block[1]
            else:
                ret[chr].append((start, block[0]))

def printTrimmedSeq(header, seq, blocks):
    for block in blocks:
        print ">%s:%d-%d" % (header, block[0], block[1])
        print seq[block[0]:block[1]]

def printTrimmedFasta(fastaFile, toTrim):
    pos = 0
    header = None
    seq = None
    for line in fastaFile:
        line = line.strip()
        if line[0] == '>':
            header = line[1:].split()[0]
            assert header in toTrim
            if seq is not None:
                printTrimmedSeq(header, seq, toTrim[header])
            seq = ""
            continue
        seq += line
    if seq is not None:
        printTrimmedSeq(header, seq, toTrim[header])

def main():
    argParser = ArgumentParser()
    argParser.add_argument("--flanking", help="Amount of flanking sequence to "
                           "leave on either end of the trimmed regions",
                           default=0)
    argParser.add_argument("--minSize", help="Minimum size of trimmed regions"
                           " before overlap (smaller regions will be dropped)",
                           default=0)
    argParser.add_argument("--complement", action="store_true", 
                           help="BED file specifies regions to "
                           "keep instead of regions to trim")
    argParser.add_argument("fasta", help="Input sequence")
    argParser.add_argument("bed", help="Regions to trim (BED format)")
    argParser.add_argument("--windowSize", type=int, default=10,
                           help="Window size for averaging "
                           "coverage over")
    argParser.add_argument("--threshold", type=int, default=1,
                           help="Coverage threshold")
    opts = argParser.parse_args()

    bedFile = open(opts.bed)
    fastaFile = open(opts.fasta)
    seqLengths = getSeqLengths(fastaFile)
    toTrim = windowFilter(opts.windowSize, opts.threshold,
                          getSeparateBedBlocks(bedFile), seqLengths)
    if opts.complement:
        toTrim = complementBlocks(toTrim, seqLengths)
    toTrim = uniquifyBlocks(toTrim, opts.flanking)
    
    # filter based on size
    toTrim = dict((k, filter(lambda x: (x[1] - x[0]) >= opts.minSize, v))
                  for k, v in toTrim.items())

    # extend blocks to include flanking regions
    toTrim = dict((k, map(lambda x: (max(x[0] - opts.flanking, 0),
                                     min(x[1] + opts.flanking, seqLengths[k])),
                          v))
                   for k, v in toTrim.items())

    fastaFile.seek(0)
    printTrimmedFasta(fastaFile, toTrim)

if __name__ == '__main__':
    main()
#!/usr/bin/env python
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""Script for running an all against all (including self) set of alignments on a set of input
sequences. Uses the toil framework to parallelise the blasts.
"""
import os
import sys
from optparse import OptionParser
from sonLib.bioio import TempFileTree
from sonLib.bioio import logger
from sonLib.bioio import system, popenCatch
from sonLib.bioio import getLogLevelString
from sonLib.bioio import makeSubDir
from sonLib.bioio import catFiles
from sonLib.bioio import getTempFile, getTempDirectory
from sonLib.bioio import nameValue
from toil.job import Job
from cactus.shared.common import WritePermanentFile

class BlastOptions:
    def __init__(self, chunkSize=10000000, overlapSize=10000, 
                 lastzArguments="", compressFiles=True, realign=False, realignArguments="",
                 minimumSequenceLength=1, memory=None,
                 # Trim options for trimming ingroup seqs:
                 trimFlanking=10, trimMinSize=20,
                 trimWindowSize=10, trimThreshold=1,
                 # Trim options for trimming outgroup seqs (options
                 # other than flanking sequence will need a check to
                 # remove alignments that don't qualify)
                 # HACK: outgroup flanking is only set so high by
                 # default because it's needed for the tests (which
                 # don't use realign.)
                 trimOutgroupFlanking=2000):
        """Class defining options for blast
        """
        self.chunkSize = chunkSize
        self.overlapSize = overlapSize
        
        if realign:
            self.blastString = "cactus_lastz --format=cigar %s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] | cactus_realign %s SEQ_FILE_1 SEQ_FILE_2 > CIGARS_FILE"  % (lastzArguments, realignArguments) 
        else:
            self.blastString = "cactus_lastz --format=cigar %s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE"  % lastzArguments 
        if realign:
            self.selfBlastString = "cactus_lastz --format=cigar %s SEQ_FILE[multiple][nameparse=darkspace] SEQ_FILE[nameparse=darkspace] --notrivial | cactus_realign %s SEQ_FILE > CIGARS_FILE" % (lastzArguments, realignArguments)
        else:
            self.selfBlastString = "cactus_lastz --format=cigar %s SEQ_FILE[multiple][nameparse=darkspace] SEQ_FILE[nameparse=darkspace] --notrivial > CIGARS_FILE" % lastzArguments
        self.compressFiles = compressFiles
        self.minimumSequenceLength = 10
        self.memory = memory
        self.trimFlanking = trimFlanking
        self.trimMinSize = trimMinSize
        self.trimThreshold = trimThreshold
        self.trimWindowSize = trimWindowSize
        self.trimOutgroupFlanking = trimOutgroupFlanking

class BlastFlower(Job):
    """Take a reconstruction problem and generate the sequences in chunks to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, cactusDisk, flowerName, finalResultsFileID, blastOptions):
        Job.__init__(self)
        self.cactusDisk = cactusDisk
        self.flowerName = flowerName
        self.finalResultsFileID = finalResultsFileID
        self.blastOptions = blastOptions
        blastOptions.roundsOfCoordinateConversion = 2
        
    def run(self, fileStore):
        chunksDir = getTempDirectory(rootDir=fileStore.getLocalTempDir())
        chunks = [ chunk for chunk in popenCatch("cactus_blast_chunkFlowerSequences %s '%s' %s %i %i %i %s" % \
                                                          (getLogLevelString(), self.cactusDisk, self.flowerName, 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          self.blastOptions.minimumSequenceLength,
                                                          chunksDir)).split("\n") if chunk != "" ]
        logger.info("Broken up the flowers into individual 'chunk' files")
        chunkIDs = [fileStore.writeGlobalFile(chunk) for chunk in chunks]
        self.addChild(MakeBlastsAllAgainstAll(self.blastOptions, chunkIDs, self.finalResultsFileID))
class BlastSequencesAllAgainstAllWrapper(Job):
    """Runs BlastSequencesAllAgainstAll on permanent files (not in the fileStore)
    """
    def __init__(self, sequenceFiles1, finalResultsFile, blastOptions):
        Job.__init__(self)
        self.sequenceFiles1 = sequenceFiles1
        self.finalResultsFile = finalResultsFile
        self.blastOptions = blastOptions
        blastOptions.roundsOfCoordinateConversion = 1
    def run(self, fileStore):
        sequenceFileIDs1 = [fileStore.writeGlobalFile(path) for path in self.sequenceFiles1]
        finalResultsFileID = fileStore.getEmptyFileStoreID()
        self.addChild(BlastSequencesAllAgainstAll(sequenceFileIDs1, finalResultsFileID, self.blastOptions))
        self.addFollowOn(WritePermanentFile(finalResultsFileID, self.finalResultsFile))

class BlastSequencesAllAgainstAll(Job):
    """Take a set of sequences, chunks them up and blasts them.
    """
    def __init__(self, sequenceFileIDs1, finalResultsFileID, blastOptions):
        Job.__init__(self)
        self.sequenceFileIDs1 = sequenceFileIDs1
        self.finalResultsFileID = finalResultsFileID
        self.blastOptions = blastOptions
        blastOptions.roundsOfCoordinateConversion = 1
    
    def getChunks(self, sequenceFiles, chunksDir):
        return [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i %i %s %s" % \
                                                          (getLogLevelString(), 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          chunksDir,
                                                          " ".join(sequenceFiles))).split("\n") if chunk != "" ]
        
    def run(self, fileStore):
        sequenceFiles1 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs1]
        chunks = self.getChunks(sequenceFiles1, makeSubDir(os.path.join(fileStore.getLocalTempDir(), "chunks")))
        logger.info("Broken up the sequence files into individual 'chunk' files")
        chunkIDs = [fileStore.writeGlobalFile(chunk) for chunk in chunks]
        self.addChild(MakeBlastsAllAgainstAll(self.blastOptions, chunkIDs, self.finalResultsFileID))
        
class MakeBlastsAllAgainstAll(Job):
    """Breaks up the inputs into bits and builds a bunch of alignment jobs.
    """
    def __init__(self, blastOptions, chunkIDs, finalResultsFileID):
        Job.__init__(self)
        self.blastOptions = blastOptions
        self.chunkIDs = chunkIDs
        self.finalResultsFileID = finalResultsFileID
        
    def run(self, fileStore):
        #Avoid compression if just one chunk
        self.blastOptions.compressFiles = self.blastOptions.compressFiles and len(self.chunkIDs) > 2
        resultsFileIDs = [fileStore.getEmptyFileStoreID() for i in xrange(len(self.chunkIDs))]

        for chunkID, resID in zip(self.chunkIDs, resultsFileIDs):
            self.addChild(RunSelfBlast(self.blastOptions, chunkID, resID))
        logger.info("Made the list of self blasts")
        #Setup job to make all-against-all blasts
        self.addFollowOn(MakeBlastsAllAgainstAll2(self.blastOptions, self.chunkIDs, resultsFileIDs, self.finalResultsFileID))
    
class MakeBlastsAllAgainstAll2(MakeBlastsAllAgainstAll):
        def __init__(self, blastOptions, chunkIDs, resultsFileIDs, finalResultsFileID):
            MakeBlastsAllAgainstAll.__init__(self, blastOptions, chunkIDs, finalResultsFileID)
            self.resultsFileIDs = resultsFileIDs
           
        def run(self, fileStore):
            #Make the list of blast jobs.
            for i in xrange(0, len(self.chunkIDs)):
                for j in xrange(i+1, len(self.chunkIDs)):
                    resultsFileID = fileStore.getEmptyFileStoreID()
                    self.resultsFileIDs.append(resultsFileID)
                    self.addChild(RunBlast(self.blastOptions, self.chunkIDs[i], self.chunkIDs[j], resultsFileID))
            logger.info("Made the list of all-against-all blasts")
            #Set up the job to collate all the results
            self.addFollowOn(CollateBlasts(self.finalResultsFileID, self.resultsFileIDs))
            
class BlastSequencesAgainstEachOtherWrapper(BlastSequencesAllAgainstAllWrapper):
    def __init__(self, sequenceFiles1, sequenceFiles2, finalResultsFile, blastOptions):
        BlastSequencesAllAgainstAllWrapper.__init__(self, sequenceFiles1, finalResultsFile, blastOptions)
        self.sequenceFiles2 = sequenceFiles2
    def run(self, fileStore):
        sequenceFileIDs1 = [fileStore.writeGlobalFile(path) for path in self.sequenceFiles1]
        sequenceFileIDs2 = [fileStore.writeGlobalFile(path) for path in self.sequenceFiles2]
        finalResultsFileID = fileStore.getEmptyFileStoreID()
        self.addChild(BlastSequencesAgainstEachOther(sequenceFileIDs1, sequenceFileIDs2, finalResultsFileID, self.blastOptions))
        self.addFollowOn(WritePermanentFile(finalResultsFileID, self.finalResultsFile))

class BlastSequencesAgainstEachOther(BlastSequencesAllAgainstAll):
    """Take two sets of sequences, chunks them up and blasts one set against the other.
    """
    def __init__(self, sequenceFileIDs1, sequenceFileIDs2, finalResultsFileID, blastOptions):
        BlastSequencesAllAgainstAll.__init__(self, sequenceFileIDs1, finalResultsFileID, blastOptions)
        self.sequenceFileIDs2 = sequenceFileIDs2
        
    def run(self, fileStore):
        sequenceFiles1 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs1]
        sequenceFiles2 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs2]
        chunks1 = self.getChunks(sequenceFiles1, makeSubDir(os.path.join(fileStore.getLocalTempDir(), "chunks1")))
        chunks2 = self.getChunks(sequenceFiles2, makeSubDir(os.path.join(fileStore.getLocalTempDir(), "chunks2")))
        chunkIDs1 = [fileStore.writeGlobalFile(chunk) for chunk in chunks1]
        chunkIDs2 = [fileStore.writeGlobalFile(chunk) for chunk in chunks2]
        resultsFileIDs = []
        #Make the list of blast jobs.
        for chunkID1 in chunkIDs1:
            for chunkID2 in chunkIDs2:
                resultsFileID = fileStore.getEmptyFileStoreID()
                resultsFileIDs.append(resultsFileID)
                #TODO: Make the compression work
                self.blastOptions.compressFiles = False
                self.addChild(RunBlast(self.blastOptions, chunkID1, chunkID2, resultsFileID))
        logger.info("Made the list of blasts")
        #Set up the job to collate all the results
        self.addFollowOn(CollateBlasts(self.finalResultsFileID, resultsFileIDs))

class BlastIngroupsAndOutgroupsWrapper(Job):
    def __init__(self, blastOptions, ingroupSequenceFiles,
                 outgroupSequenceFiles, finalResultsFile,
                 outgroupFragmentsDir):
        Job.__init__(self, memory = blastOptions.memory)
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 1
        self.ingroupSequenceFiles = ingroupSequenceFiles
        self.outgroupSequenceFiles = outgroupSequenceFiles
        self.finalResultsFile = finalResultsFile
        self.outgroupFragmentsDir = outgroupFragmentsDir
    def run(self, fileStore):
        #move sequence files into the fileStore and then run BlastIngroupsAndOutgroups
        ingroupSequenceFileIDs = [fileStore.writeGlobalFile(path) for path in self.ingroupSequenceFiles]
        outgroupSequenceFileIDs = [fileStore.writeGlobalFile(path) for path in self.outgroupSequenceFiles]
        finalResultsFileID = fileStore.getEmptyFileStoreID()
        outgroupFragmentIDs = [fileStore.getEmptyFileStoreID() for i in xrange(len(outgroupSequenceFileIDs))]
        self.addChild(BlastIngroupsAndOutgroups(self.blastOptions, ingroupSequenceFileIDs, outgroupSequenceFileIDs, finalResultsFileID, outgroupFragmentIDs))
        self.addFollowOn(WritePermanentFile(finalResultsFileID, self.finalResultsFile))

        try:
            os.makedirs(self.outgroupFragmentsDir)
        except os.error:
            pass
        for outgroupID, outgroupPath in zip(outgroupSequenceFileIDs, self.outgroupSequenceFiles):
            self.addFollowOn(WritePermanentFile(outgroupID, os.path.basename(outgroupPath)))


class BlastIngroupsAndOutgroups(Job):
    """Blast ingroup sequences against each other, and against the given
    outgroup sequences in succession. The next outgroup is only
    aligned against the regions that are not found in the previous
    outgroup.
    """
    def __init__(self, blastOptions, ingroupSequenceFileIDs,
                 outgroupSequenceFileIDs, finalResultsFileID,
                 outgroupFragmentIDs):
        Job.__init__(self, memory = blastOptions.memory)
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 1
        self.ingroupSequenceFileIDs = ingroupSequenceFileIDs
        self.outgroupSequenceFileIDs = outgroupSequenceFileIDs
        self.finalResultsFileID = finalResultsFileID
        self.outgroupFragmentIDs = outgroupFragmentIDs

    def run(self, fileStore):
        fileStore.logToMaster("Blasting ingroups vs outgroups to file %s" % (self.finalResultsFileID))
        ingroupResultsID = fileStore.getEmptyFileStoreID()
        self.addChild(BlastSequencesAllAgainstAll(self.ingroupSequenceFileIDs,
                                                        ingroupResultsID,
                                                        self.blastOptions))
        outgroupResultsID = fileStore.getEmptyFileStoreID()
        self.addChild(BlastFirstOutgroup(self.ingroupSequenceFileIDs,
                                               self.ingroupSequenceFileIDs,
                                               self.outgroupSequenceFileIDs,
                                               self.outgroupFragmentIDs,
                                               outgroupResultsID,
                                               self.blastOptions, 1))
        self.addFollowOn(CollateBlasts(self.finalResultsFileID,
                                             [ingroupResultsID, outgroupResultsID]))

class BlastFirstOutgroup(Job):
    """Blast the given sequence(s) against the first of a succession of
    outgroups, only aligning fragments that haven't aligned to the
    previous outgroups. Then recurse on the other outgroups.
    """
    def __init__(self, untrimmedSequenceFileIDs, sequenceFileIDs,
                 outgroupSequenceFileIDs, outgroupFragmentIDs, outputFileID,
                 blastOptions, outgroupNumber):
        Job.__init__(self, memory=blastOptions.memory)
        self.untrimmedSequenceFileIDs = untrimmedSequenceFileIDs
        self.sequenceFileIDs = sequenceFileIDs
        self.outgroupSequenceFileIDs = outgroupSequenceFileIDs
        self.outgroupFragmentIDs = outgroupFragmentIDs
        self.outputFileID = outputFileID
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber

    def run(self, fileStore):
        logger.info("Blasting ingroup sequences %s to outgroup %s" % (self.sequenceFileIDs, self.outgroupSequenceFileIDs[0]))
        blastResultsID = fileStore.getEmptyFileStoreID()
        self.addChild(BlastSequencesAgainstEachOther(self.sequenceFileIDs,
                                                           [self.outgroupSequenceFileIDs[0]],
                                                           blastResultsID,
                                                           self.blastOptions))
        self.addFollowOn(TrimAndRecurseOnOutgroups(self.untrimmedSequenceFileIDs,
                                                         self.sequenceFileIDs,
                                                         self.outgroupSequenceFileIDs,
                                                         self.outgroupFragmentIDs,
                                                         blastResultsID,
                                                         self.outputFileID,
                                                         self.blastOptions,
                                                         self.outgroupNumber))

class TrimAndRecurseOnOutgroups(Job):
    def __init__(self, untrimmedSequenceFileIDs, sequenceFileIDs,
                 outgroupSequenceFileIDs, outgroupFragmentIDs,
                 mostRecentResultsFileID, outputFileID, blastOptions,
                 outgroupNumber):
        Job.__init__(self)
        self.untrimmedSequenceFileIDs = untrimmedSequenceFileIDs
        self.sequenceFileIDs = sequenceFileIDs
        self.outgroupSequenceFileIDs = outgroupSequenceFileIDs
        self.outgroupFragmentIDs = outgroupFragmentIDs
        self.mostRecentResultsFileID = mostRecentResultsFileID
        self.outputFileID = outputFileID
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber

    def run(self, fileStore):
        # Trim outgroup, convert outgroup coordinates, and add to
        # outgroup fragments dir
        outgroupSequenceFiles = [fileStore.readGlobalFile(fileID) for fileID in self.outgroupSequenceFileIDs]
        sequenceFiles = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs]
        untrimmedSequenceFiles = [fileStore.readGlobalFile(fileID) for fileID in self.untrimmedSequenceFileIDs]
        mostRecentResultsFile = fileStore.readGlobalFile(self.mostRecentResultsFileID)

        mostRecentResultsFile = fileStore.readGlobalFile(self.mostRecentResultsFileID)

        trimmedOutgroup = getTempFile(rootDir=fileStore.getLocalTempDir())
        outgroupCoverage = getTempFile(rootDir=fileStore.getLocalTempDir())
        calculateCoverage(outgroupSequenceFiles[0],
                          mostRecentResultsFile, outgroupCoverage)
        # The windowSize and threshold are fixed at 1: anything more
        # and we will run into problems with alignments that aren't
        # covered in a matching trimmed sequence.
        trimGenome(outgroupSequenceFiles[0], outgroupCoverage,
                   trimmedOutgroup, flanking=self.blastOptions.trimOutgroupFlanking,
                   windowSize=1, threshold=1)
        outgroupConvertedResultsFile = getTempFile(rootDir=fileStore.getLocalTempDir())
        system("cactus_upconvertCoordinates.py %s %s 1 > %s" %\
               (trimmedOutgroup, mostRecentResultsFile,
                outgroupConvertedResultsFile))
        
        fileStore.updateGlobalFile(self.outgroupFragmentIDs[0], trimmedOutgroup)

        # Report coverage of the latest outgroup on the trimmed ingroups.
        for trimmedIngroupSequence, ingroupSequence in zip(sequenceFiles, untrimmedSequenceFiles):
            tmpIngroupCoverage = getTempFile(rootDir=fileStore.getLocalTempDir())
            calculateCoverage(trimmedIngroupSequence, mostRecentResultsFile,
                              tmpIngroupCoverage)
            fileStore.logToMaster("Coverage on %s from outgroup #%d, %s: %s%% (current ingroup length %d, untrimmed length %d). Outgroup trimmed to %d bp from %d" % (os.path.basename(ingroupSequence), self.outgroupNumber, os.path.basename(outgroupSequenceFiles[0]), percentCoverage(trimmedIngroupSequence, tmpIngroupCoverage), sequenceLength(trimmedIngroupSequence), sequenceLength(ingroupSequence), sequenceLength(trimmedOutgroup), sequenceLength(outgroupSequenceFiles[0])))


        # Convert the alignments' ingroup coordinates.
        ingroupConvertedResultsFile = getTempFile(rootDir=fileStore.getLocalTempDir())
        if self.sequenceFileIDs == self.untrimmedSequenceFileIDs:
            # No need to convert ingroup coordinates on first run.
            system("cp %s %s" % (outgroupConvertedResultsFile,
                                 ingroupConvertedResultsFile))
        else:
            system("cactus_blast_convertCoordinates --onlyContig1 %s %s 1" % (
                outgroupConvertedResultsFile, ingroupConvertedResultsFile))
        
        # Append the latest results to the accumulated outgroup coverage file
        outputFile = fileStore.readGlobalFile(self.outputFileID)
        with open(ingroupConvertedResultsFile) as results:
            with open(outputFile, 'a') as output:
                output.write(results.read())
        fileStore.updateGlobalFile(self.outputFileID, outputFile)
        os.remove(outgroupConvertedResultsFile)
        os.remove(ingroupConvertedResultsFile)
        os.remove(outgroupCoverage)
        os.remove(trimmedOutgroup)

        # Report coverage of all outgroup alignments so far on the ingroups.
        ingroupCoverageFiles = []
        for ingroupSequence in untrimmedSequenceFiles:
            tmpIngroupCoverage = getTempFile(rootDir=fileStore.getLocalTempDir())
            calculateCoverage(ingroupSequence, outputFile,
                              tmpIngroupCoverage)
            fileStore.logToMaster("Cumulative coverage of %d outgroups on ingroup %s: %s" % (self.outgroupNumber, os.path.basename(ingroupSequence), percentCoverage(ingroupSequence, tmpIngroupCoverage)))
            ingroupCoverageFiles.append(tmpIngroupCoverage)

        # Trim ingroup seqs and recurse on the next outgroup.

        # TODO: Optionally look at coverage on ingroup vs. outgroup,
        # and if coverage is >1 among ingroups but 1 in outgroups,
        # look for it in the next outgroup as well. Would require
        # doing self blast first and sending the alignments here.
        # (Probably needs an extra option in cactus coverage to only
        # count self-alignments, since we need to cut the ingroup
        # sequences in question and not something aligning to both of
        # them.)
        # Could also just ignore the coverage on the outgroup to
        # start, since the fraction of duplicated sequence will be
        # relatively small.
        if len(outgroupSequenceFiles) > 1:
            trimmedSeqs = []
            # Use the accumulated results so far to trim away the
            # aligned parts of the ingroups.
            for i, sequenceFile in enumerate(untrimmedSequenceFiles):
                coverageFile = ingroupCoverageFiles[i]

                trimmed = getTempFile(rootDir=fileStore.getLocalTempDir())
                trimGenome(sequenceFile, coverageFile, trimmed,
                           complement=True, flanking=self.blastOptions.trimFlanking,
                           minSize=self.blastOptions.trimMinSize,
                           threshold=self.blastOptions.trimThreshold,
                           windowSize=self.blastOptions.trimWindowSize)
                trimmedSeqs.append(trimmed)
            for seqID, trimmedSeq in zip(self.sequenceFileIDs, trimmedSeqs):
                fileStore.updateGlobalFile(seqID, trimmedSeq)
            self.addChild(BlastFirstOutgroup(self.untrimmedSequenceFileIDs,
                                                   self.sequenceFileIDs,
                                                   self.outgroupSequenceFileIDs[1:],
                                                   self.outgroupFragmentIDs[1:],
                                                   self.outputFileID,
                                                   self.blastOptions,
                                                   self.outgroupNumber + 1))

def compressFastaFile(fileName):
    """Compress a fasta file.
    """
    system("bzip2 --keep --fast %s" % fileName)
    return fileName + ".bz2"
        
class RunSelfBlast(Job):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFileID, resultsFileID):
        Job.__init__(self, memory=blastOptions.memory)
        self.blastOptions = blastOptions
        self.seqFileID = seqFileID
        self.resultsFileID = resultsFileID
    
    def run(self, fileStore):   
        seqFile = fileStore.readGlobalFile(self.seqFileID)
        tempResultsFile = getTempFile(rootDir=fileStore.getLocalTempDir())
        command = self.blastOptions.selfBlastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE", seqFile)
        system(command)
        resultsFile = getTempFile(rootDir=fileStore.getLocalTempDir())
        system("cactus_blast_convertCoordinates %s %s %i" % (tempResultsFile, resultsFile, self.blastOptions.roundsOfCoordinateConversion))
        if self.blastOptions.compressFiles:
            seqFile = compressFastaFile(seqFile)
            fileStore.updateGlobalFile(self.seqFileID, seqFile)
        fileStore.updateGlobalFile(self.resultsFileID, resultsFile)
        logger.info("Ran the self blast okay")

def decompressFastaFile(fileName, tempFileName):
    """Copies the file from the central dir to a temporary file, returning the temp file name.
    """
    system("bunzip2 --stdout %s > %s" % (fileName, tempFileName))
    return tempFileName
    
class RunBlast(Job):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFileID1, seqFileID2, resultsFileID):
        Job.__init__(self, memory=blastOptions.memory)
        self.blastOptions = blastOptions
        self.seqFileID1 = seqFileID1
        self.seqFileID2 = seqFileID2
        self.resultsFileID = resultsFileID
    
    def run(self, fileStore):
        seqFile1 = fileStore.readGlobalFile(self.seqFileID1)
        seqFile2 = fileStore.readGlobalFile(self.seqFileID2)
        if self.blastOptions.compressFiles:
            seqFile1 = decompressFastaFile(seqFile1, getTempFile(rootDir=fileStore.getLocalTempDir()))
            seqFile2 = decompressFastaFile(seqFile2, getTempFile(rootDir=fileStore.getLocalTempDir()))
        tempResultsFile = getTempFile(rootDir=fileStore.getLocalTempDir())
        command = self.blastOptions.blastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE_1", seqFile1).replace("SEQ_FILE_2", seqFile2)
        system(command)
        resultsFile = getTempFile(rootDir=fileStore.getLocalTempDir())
        system("cactus_blast_convertCoordinates %s %s %i" % (tempResultsFile, resultsFile, self.blastOptions.roundsOfCoordinateConversion))
        fileStore.updateGlobalFile(self.resultsFileID, resultsFile)
        logger.info("Ran the blast okay")

class CollateBlasts(Job):
    """Collates all the blasts into a single alignments file.
    """
    def __init__(self, finalResultsFileID, resultsFileIDs):
        Job.__init__(self)
        self.finalResultsFileID = finalResultsFileID
        self.resultsFileIDs = resultsFileIDs
    
    def run(self, fileStore):
        resultsFiles = [fileStore.readGlobalFile(fileID) for fileID in self.resultsFileIDs]
        finalResultsFile = getTempFile(rootDir=fileStore.getLocalTempDir())
        catFiles(resultsFiles, finalResultsFile)
        fileStore.updateGlobalFile(self.finalResultsFileID, finalResultsFile)
        logger.info("Collated the alignments to the file: %s",  self.finalResultsFileID)
        
class SortCigarAlignmentsInPlace(Job):
    """Sorts an alignment file in place.
    """
    def __init__(self, cigarFile):
        Job.__init__(self)
        self.cigarFile = cigarFile
    
    def run(self, fileStore):
        tempResultsFile = os.path.join(self.getLocalTempDir(), "tempResults.cig")
        system("cactus_blast_sortAlignments %s %s %i" % (getLogLevelString(), self.cigarFile, tempResultsFile))
        logger.info("Sorted the alignments okay")
        system("mv %s %s" % (tempResultsFile, self.cigarFile))

def sequenceLength(sequenceFile):
    """Get the total # of bp from a fasta file."""
    seqLength = 0
    for line in open(sequenceFile):
        line = line.strip()
        if line == '' or line[0] == '>':
            continue
        seqLength += len(line)
    return seqLength

def percentCoverage(sequenceFile, coverageFile):
    """Get the % coverage of a sequence from a coverage file."""
    sequenceLen = sequenceLength(sequenceFile)
    if sequenceLen == 0:
        return 0
    coverage = popenCatch("awk '{ total += $3 - $2 } END { print total }' %s" % coverageFile)
    if coverage.strip() == '': # No coverage lines
        return 0
    return 100*float(coverage)/sequenceLen

def calculateCoverage(sequenceFile, cigarFile, outputFile):
    logger.info("Calculating coverage of cigar file %s on %s, writing to %s" % (
        cigarFile, sequenceFile, outputFile))
    system("cactus_coverage %s %s > %s" % (sequenceFile,
                                           cigarFile,
                                           outputFile))

def trimGenome(sequenceFile, coverageFile, outputFile, complement=False,
               flanking=0, minSize=1, windowSize=10, threshold=1):
    system("cactus_trimSequences.py %s %s %s %s %s %s %s > %s" % (
        nameValue("complement", complement, valueType=bool),
        nameValue("flanking", flanking), nameValue("minSize", minSize),
        nameValue("windowSize", windowSize), nameValue("threshold", threshold),
        sequenceFile, coverageFile, outputFile))

def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
    Job.Runner.addToilOptions(parser)
    blastOptions = BlastOptions()
    
    #output stuff
    parser.add_option("--cigars", dest="cigarFile", 
                      help="File to write cigars in",
                      default="cigarFile.txt")
    
    parser.add_option("--chunkSize", dest="chunkSize", type="int", 
                     help="The size of chunks passed to lastz (must be at least twice as big as overlap)",
                     default=blastOptions.chunkSize)
    
    parser.add_option("--overlapSize", dest="overlapSize", type="int",
                     help="The size of the overlap between the chunks passed to lastz (min size 2)",
                     default=blastOptions.overlapSize)
    
    parser.add_option("--blastString", dest="blastString", type="string",
                     help="The default string used to call the blast program. \
Must contain three strings: SEQ_FILE_1, SEQ_FILE_2 and CIGARS_FILE which will be \
replaced with the two sequence files and the results file, respectively",
                     default=blastOptions.blastString)
    
    parser.add_option("--selfBlastString", dest="selfBlastString", type="string",
                     help="The default string used to call the blast program for self alignment. \
Must contain three strings: SEQ_FILE and CIGARS_FILE which will be \
replaced with the the sequence file and the results file, respectively",
                     default=blastOptions.selfBlastString)
   
    parser.add_option("--compressFiles", dest="compressFiles", action="store_false",
                      help="Turn of bz2 based file compression of sequences for I/O transfer", 
                      default=blastOptions.compressFiles)
    
    parser.add_option("--lastzMemory", dest="memory", type="int",
                      help="Lastz memory (in bytes)", 
                      default=blastOptions.memory)
    
    parser.add_option("--trimFlanking", type=int, help="Amount of flanking sequence to leave on trimmed ingroup sequences", default=blastOptions.trimFlanking)
    parser.add_option("--trimMinSize", type=int, help="Minimum size, before adding flanking sequence, of ingroup sequence to align against the next outgroup", default=blastOptions.trimMinSize)
    parser.add_option("--trimThreshold", type=int, help="Coverage threshold for an ingroup region to not be aligned against the next outgroup", default=blastOptions.trimThreshold)
    parser.add_option("--trimWindowSize", type=int, help="Windowing size to integrate ingroup coverage over", default=blastOptions.trimWindowSize)
    parser.add_option("--trimOutgroupFlanking", type=int, help="Amount of flanking sequence to leave on trimmed outgroup sequences", default=blastOptions.trimOutgroupFlanking)
    

    parser.add_option("--test", dest="test", action="store_true",
                      help="Run doctest unit tests")
    
    parser.add_option("--targetSequenceFiles", dest="targetSequenceFiles", type="string",
                     help="Sequences to align against the input sequences against. If these are not provided then the input sequences are aligned against each other.",
                     default=None)

    parser.add_option("--ingroups", type=str, default=None,
                      help="Ingroups to align (comma-separated) (--outgroups "
                      "must be provided as well")

    parser.add_option("--outgroups", type=str, default=None,
                      help="Outgroups to align (comma-separated) (--ingroups "
                      "must be provided as well")

    parser.add_option("--outgroupFragmentsDir", type=str,
                      default="outgroupFragments/", help= "Directory to "
                      "store outgroup fragments in")

    options, args = parser.parse_args()
    if options.test:
        _test()

    if (options.ingroups is not None) ^ (options.outgroups is not None):
        raise RuntimeError("--ingroups and --outgroups must be provided "
                           "together")
    if options.ingroups:
        firstJob = BlastIngroupsAndOutgroupsWrapper(options,
                                                options.ingroups.split(','),
                                                options.outgroups.split(','),
                                                options.cigarFile,
                                                options.outgroupFragmentsDir)
    elif options.targetSequenceFiles == None:
        firstJob = BlastSequencesAllAgainstAllWrapper(args, options.cigarFile, options)
    else:
        firstJob = BlastSequencesAgainstEachOtherWrapper(args, options.targetSequenceFiles.split(), options.cigarFile, options)
    Job.Runner.startToil(firstJob, options)

def _test():
    import doctest
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.blast.cactus_blast import *
    main()

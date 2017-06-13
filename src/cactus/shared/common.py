#!/usr/bin/env python
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Wrapper functions for assisting in running the various programs of the cactus package.
"""

import os
import sys
import shutil
import subprocess32
import logging
import uuid
import json

from toil.lib.bioio import logger
from toil.lib.bioio import system
from toil.lib.bioio import getLogLevelString

from toil.job import Job

from sonLib.bioio import nameValue
from sonLib.bioio import popenCatch

from cactus.shared.version import cactus_commit

_log = logging.getLogger(__name__)

subprocess32._has_poll = False

def makeURL(path):
    if not (path.startswith("file:") or path.startswith("s3:") or path.startswith("http:")):
        return "file://" + os.path.abspath(path)
    else:
        return path

def catFiles(filesToCat, catFile):
    """Cats a bunch of files into one file. Ensures a no more than maxCat files
    are concatenated at each step.
    """
    if len(filesToCat) == 0: #We must handle this case or the cat call will hang waiting for input
        open(catFile, 'w').close()
        return
    maxCat = 25
    system("cat %s > %s" % (" ".join(filesToCat[:maxCat]), catFile))
    filesToCat = filesToCat[maxCat:]
    while len(filesToCat) > 0:
        system("cat %s >> %s" % (" ".join(filesToCat[:maxCat]), catFile))
        filesToCat = filesToCat[maxCat:]

def nameValue(name, value, valueType=str, quotes=False):
    """Little function to make it easier to make name value strings for commands.
    """
    if valueType == bool:
        if value:
            return "--%s" % name
        return ""
    if value is None:
        return ""
    if quotes:
        return "--%s '%s'" % (name, valueType(value))
    return "--%s %s" % (name, valueType(value))

def cactusRootPath():
    """
    function for finding external location
    """
    import cactus
    i = os.path.abspath(cactus.__file__)
    return os.path.split(i)[0]

def getLogLevelString2(logLevelString):
    """Gets the log level string for the binary
    """
    if logLevelString == None:
        return getLogLevelString()
    return logLevelString

def getOptionalAttrib(node, attribName, typeFn=None, default=None):
    """Get an optional attrib, or default if not set or node is None
    """
    if node != None and node.attrib.has_key(attribName):
        if typeFn != None:
            if typeFn == bool:
                return bool(int(node.attrib[attribName]))
            return typeFn(node.attrib[attribName])
        return node.attrib[attribName]
    return default

def findRequiredNode(configNode, nodeName):
    """Retrieve an xml node, complain if it's not there."""
    nodes = configNode.findall(nodeName)
    if nodes == None:
        raise RuntimeError("Could not find any nodes with name %s in %s node" % (nodeName, configNode))
    assert len(nodes) == 1, "More than 1 node for %s in config XML" % nodeName
    return nodes[0]

#############################################
#############################################
#Following used to gather the names of flowers
#in problems
#############################################
#############################################  

def readFlowerNames(flowerStrings):
    ret = []
    for line in flowerStrings.split("\n"):
        if line == '':
            continue
        flowersAndSizes = line[1:].split()
        numFlowers = flowersAndSizes[0]
        flowers = []
        sizes = []
        currentlyAFlower = True
        for token in flowersAndSizes[1:]:
            if token == 'a' or token == 'b':
                flowers += [token]
            elif currentlyAFlower:
                flowers += [token]
                currentlyAFlower = False
            else:
                sizes += [int(token)]
                currentlyAFlower = True
        assert len(sizes) == int(numFlowers)
        ret += [(bool(int(line[0])), " ".join([numFlowers] + flowers), sizes)]
    return ret

def runCactusGetFlowers(cactusDiskDatabaseString, flowerNames,
                        jobName=None, features=None, fileStore=None,
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, 
                        maxSequenceSizeOfSecondaryFlowerGrouping=-1, 
                        logLevel=None):
    """Gets a list of flowers attached to the given flower. 
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStrings = cactus_call(check_output=True, stdin_string=flowerNames,
                                parameters=["cactus_workflow_getFlowers"],
                                option_string="%s '%s' %i %i %i" % 
                                            (logLevel, cactusDiskDatabaseString,
                                            minSequenceSizeOfFlower,
                                            maxSequenceSizeOfFlowerGrouping,
                                            maxSequenceSizeOfSecondaryFlowerGrouping),
                                job_name=jobName,
                                features=features,
                                fileStore=fileStore)

    l = readFlowerNames(flowerStrings)
    return l

def runCactusExtendFlowers(cactusDiskDatabaseString, flowerNames, 
                        jobName=None, features=None, fileStore=None,
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, 
                        maxSequenceSizeOfSecondaryFlowerGrouping=-1, 
                        logLevel=None):
    """Extends the terminal groups in the cactus and returns the list
    of their child flowers with which to pass to core.
    The order of the flowers is by ascending depth first discovery time.
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStrings = cactus_call(check_output=True, stdin_string=flowerNames,
                                parameters=["cactus_workflow_extendFlowers"],
                                option_string="%s '%s' %i %i %i" %
                                            (logLevel,
                                            cactusDiskDatabaseString,
                                            minSequenceSizeOfFlower,
                                            maxSequenceSizeOfFlowerGrouping,
                                            maxSequenceSizeOfSecondaryFlowerGrouping),
                                job_name=jobName,
                                features=features,
                                fileStore=fileStore)

    l = readFlowerNames(flowerStrings)
    return l

def encodeFlowerNames(flowerNames):
    if len(flowerNames) == 0:
        return "0"
    return "%i %s" % (len(flowerNames), " ".join([ str(flowerNames[0]) ] + [ str(flowerNames[i] - flowerNames[i-1]) for i in xrange(1, len(flowerNames)) ]))
    
def decodeFirstFlowerName(encodedFlowerNames):
    tokens = encodedFlowerNames.split()
    if int(tokens[0]) == 0:
        return None
    if tokens[1] == 'b':
        return int(tokens[2])
    return int(tokens[1])

def runCactusSplitFlowersBySecondaryGrouping(flowerNames):
    """Splits a list of flowers into smaller lists.
    """
    flowerNames = flowerNames.split()
    flowerGroups = []
    stack = []
    overlarge = False
    name = 0
    for i in flowerNames[1:]:
        if i != '':
            if i in ('a', 'b'):
                if len(stack) > 0:
                    flowerGroups.append((overlarge, encodeFlowerNames(stack))) #b indicates the stack is overlarge
                    stack = []
                overlarge = i == 'b'
            else:
                name = int(i) + name
                stack.append(name)
    if len(stack) > 0:
        flowerGroups.append((overlarge, encodeFlowerNames(stack)))
    return flowerGroups

#############################################
#############################################
#All the following provide command line wrappers
#for core programs in the cactus pipeline.
#############################################
#############################################  

def runCactusSetup(cactusDiskDatabaseString, sequences, 
                   newickTreeString, logLevel=None, outgroupEvents=None,
                   makeEventHeadersAlphaNumeric=None):
    logLevel = getLogLevelString2(logLevel)
    outgroupEvents = nameValue("outgroupEvents", outgroupEvents, str, quotes=True)
    makeEventHeadersAlphaNumeric=nameValue("makeEventHeadersAlphaNumeric", makeEventHeadersAlphaNumeric, bool)
    masterMessages = cactus_call(check_output=True,
                                 parameters=["cactus_setup"] + sequences,
                                 option_string="--speciesTree '%s' --cactusDisk '%s' --logLevel %s %s %s" % (newickTreeString, cactusDiskDatabaseString, logLevel, outgroupEvents, makeEventHeadersAlphaNumeric))
    
    logger.info("Ran cactus setup okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]
    


def runConvertAlignmentsToInternalNames(cactusDiskString, alignmentsFile, outputFile, flowerName, isBedFile = False):
    isBedFile = nameValue("bed", isBedFile, bool)
    cactus_call(stdin_string=encodeFlowerNames((flowerName,)),
                option_string="--cactusDisk '%s'" % cactusDiskString,
                parameters=["cactus_convertAlignmentsToInternalNames",
                            alignmentsFile, outputFile,
                            isBedFile])
    
def runStripUniqueIDs(cactusDiskString):
    cactus_call(option_string="--cactusDisk '%s'" % cactusDiskString,
                parameters=["cactus_stripUniqueIDs"])
    

def runCactusCaf(cactusDiskDatabaseString, alignments,
                 flowerNames=encodeFlowerNames((0,)),
                 logLevel=None, 
                 writeDebugFiles=False,
                 annealingRounds=None,
                 deannealingRounds=None,
                 trim=None,
                 minimumTreeCoverage=None,
                 blockTrim=None,
                 minimumBlockDegree=None,
                 minimumIngroupDegree=None,
                 minimumOutgroupDegree=None,
                 alignmentFilter=None,
                 lastzArguments=None,
                 minimumSequenceLengthForBlast=None,
                 maxAdjacencyComponentSizeRatio=None,
                 constraints=None,
                 minLengthForChromosome=None,
                 proportionOfUnalignedBasesForNewChromosome=None, 
                 maximumMedianSequenceLengthBetweenLinkedEnds=None,
                 realign=None,
                 realignArguments=None,
                 phylogenyNumTrees=None,
                 phylogenyScoringMethod=None,
                 phylogenyRootingMethod=None,
                 phylogenyBreakpointScalingFactor=None,
                 phylogenySkipSingleCopyBlocks=None,
                 phylogenyMaxBaseDistance=None,
                 phylogenyMaxBlockDistance=None,
                 phylogenyDebugFile=None,
                 phylogenyKeepSingleDegreeBlocks=None,
                 phylogenyTreeBuildingMethod=None,
                 phylogenyCostPerDupPerBase=None,
                 phylogenyCostPerLossPerBase=None,
                 referenceEventHeader=None,
                 phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce=None,
                 numTreeBuildingThreads=None,
                 doPhylogeny=None,
                 removeLargestBlock=None,
                 phylogenyNucleotideScalingFactor=None,
                 minimumBlockDegreeToCheckSupport=None,
                 minimumBlockHomologySupport=None,
                 removeRecoverableChains=None,
                 minimumNumberOfSpecies=None,
                 maxRecoverableChainsIterations=None,
                 maxRecoverableChainLength=None,
                 phylogenyHomologyUnitType=None,
                 phylogenyDistanceCorrectionMethod=None,
                 features=None,
                 jobName=None,
                 fileStore=None):
    # remove annoying carriage returns in caf command line.
    cactusDiskDatabaseString = cactusDiskDatabaseString.replace('\n', '')

    if alignments:
        alignments = os.path.basename(alignments)

    logLevel = getLogLevelString2(logLevel)
    annealingRounds = nameValue("annealingRounds", annealingRounds, quotes=True)
    deannealingRounds = nameValue("deannealingRounds", deannealingRounds, quotes=True)
    trim = nameValue("trim", trim, quotes=True)
    alignments = nameValue("alignments", alignments)
    lastzArguments = nameValue("lastzArguments", lastzArguments, quotes=True)
    minimumTreeCoverage = nameValue("minimumTreeCoverage", minimumTreeCoverage, float)
    blockTrim = nameValue("blockTrim", blockTrim, int)
    minimumBlockDegree = nameValue("minimumDegree", minimumBlockDegree, int)
    minimumSequenceLengthForBlast = nameValue("minimumSequenceLengthForBlast", minimumSequenceLengthForBlast, int)
    minimumIngroupDegree = nameValue("minimumIngroupDegree", minimumIngroupDegree, int)
    minimumOutgroupDegree = nameValue("minimumOutgroupDegree", minimumOutgroupDegree, int)
    alignmentFilter = nameValue("alignmentFilter", alignmentFilter)
    maxAdjacencyComponentSizeRatio = nameValue("maxAdjacencyComponentSizeRatio", maxAdjacencyComponentSizeRatio, float)
    constraints = nameValue("constraints", constraints)
    realign = nameValue("realign", realign, bool)
    realignArguments = nameValue("realignArguments", realignArguments, quotes=True)
    phylogenyNumTrees = nameValue("phylogenyNumTrees", phylogenyNumTrees, int)
    phylogenyRootingMethod = nameValue("phylogenyRootingMethod", phylogenyRootingMethod, quotes=True)
    phylogenyScoringMethod = nameValue("phylogenyScoringMethod", phylogenyScoringMethod, quotes=True)
    phylogenyBreakpointScalingFactor = nameValue("phylogenyBreakpointScalingFactor", phylogenyBreakpointScalingFactor)
    phylogenySkipSingleCopyBlocks = nameValue("phylogenySkipSingleCopyBlocks", phylogenySkipSingleCopyBlocks, bool)
    phylogenyMaxBaseDistance = nameValue("phylogenyMaxBaseDistance", phylogenyMaxBaseDistance)
    phylogenyMaxBlockDistance = nameValue("phylogenyMaxBlockDistance", phylogenyMaxBlockDistance)
    phylogenyDebugFile = nameValue("phylogenyDebugFile", phylogenyDebugFile)
    phylogenyKeepSingleDegreeBlocks = nameValue("phylogenyKeepSingleDegreeBlocks", phylogenyKeepSingleDegreeBlocks, bool)
    phylogenyTreeBuildingMethod = nameValue("phylogenyTreeBuildingMethod", phylogenyTreeBuildingMethod)
    phylogenyCostPerDupPerBase = nameValue("phylogenyCostPerDupPerBase", phylogenyCostPerDupPerBase)
    phylogenyCostPerLossPerBase = nameValue("phylogenyCostPerLossPerBase", phylogenyCostPerLossPerBase)
    referenceEventHeader = nameValue("referenceEventHeader", referenceEventHeader, quotes=True)
    phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce = nameValue("phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce", phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce)
    numTreeBuildingThreads = nameValue("numTreeBuildingThreads", numTreeBuildingThreads)
    doPhylogeny = nameValue("phylogeny", doPhylogeny, bool)
    minimumBlockDegreeToCheckSupport = nameValue("minimumBlockDegreeToCheckSupport", minimumBlockDegreeToCheckSupport)
    minimumBlockHomologySupport = nameValue("minimumBlockHomologySupport", minimumBlockHomologySupport)
    phylogenyNucleotideScalingFactor = nameValue("phylogenyNucleotideScalingFactor", phylogenyNucleotideScalingFactor)
    removeRecoverableChains = nameValue("removeRecoverableChains", removeRecoverableChains)
    minimumNumberOfSpecies = nameValue("minimumNumberOfSpecies", minimumNumberOfSpecies, int)
    maxRecoverableChainsIterations = nameValue("maxRecoverableChainsIterations", maxRecoverableChainsIterations, int)
    maxRecoverableChainLength = nameValue("maxRecoverableChainLength", maxRecoverableChainLength, int)
    phylogenyHomologyUnitType = nameValue("phylogenyHomologyUnitType", phylogenyHomologyUnitType, quotes=True)
    phylogenyDistanceCorrectionMethod = nameValue("phylogenyDistanceCorrectionMethod", phylogenyDistanceCorrectionMethod, quotes=True)

    minLengthForChromosome = nameValue("minLengthForChromosome", minLengthForChromosome, int)
    proportionOfUnalignedBasesForNewChromosome = nameValue("proportionOfUnalignedBasesForNewChromosome", proportionOfUnalignedBasesForNewChromosome, float)
    maximumMedianSequenceLengthBetweenLinkedEnds = nameValue("maximumMedianSequenceLengthBetweenLinkedEnds", maximumMedianSequenceLengthBetweenLinkedEnds, int)

    masterMessages = cactus_call(stdin_string=flowerNames, check_output=True,
                                 option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                                 parameters=["cactus_caf",
                                             "--logLevel", logLevel, alignments, annealingRounds,
                                             deannealingRounds, 
                                             trim, minimumTreeCoverage, blockTrim, 
                                             minimumBlockDegree, minimumIngroupDegree, minimumOutgroupDegree,  
                                             alignmentFilter, lastzArguments, minimumSequenceLengthForBlast,
                                             maxAdjacencyComponentSizeRatio, constraints,
                                             minLengthForChromosome,
                                             proportionOfUnalignedBasesForNewChromosome,
                                             maximumMedianSequenceLengthBetweenLinkedEnds, realign,
                                             realignArguments, phylogenyNumTrees, phylogenyRootingMethod,
                                             phylogenyScoringMethod, phylogenyBreakpointScalingFactor,
                                             phylogenySkipSingleCopyBlocks, phylogenyMaxBaseDistance,
                                             phylogenyMaxBlockDistance, phylogenyDebugFile,
                                             phylogenyKeepSingleDegreeBlocks, phylogenyTreeBuildingMethod,
                                             phylogenyCostPerDupPerBase, phylogenyCostPerLossPerBase,
                                             referenceEventHeader,
                                             phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce,
                                             numTreeBuildingThreads, doPhylogeny,
                                             minimumBlockDegreeToCheckSupport, minimumBlockHomologySupport,
                                             phylogenyNucleotideScalingFactor, removeRecoverableChains,
                                             minimumNumberOfSpecies, phylogenyHomologyUnitType,
                                             phylogenyDistanceCorrectionMethod,
                                             maxRecoverableChainsIterations, maxRecoverableChainLength],
                                 features=features, job_name=jobName, fileStore=fileStore)
    logger.info("Ran cactus_core okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]
    
def runCactusPhylogeny(cactusDiskDatabaseString,
                       flowerNames=encodeFlowerNames((0,)),
                       logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_phylogeny",
                            "--cactusDisk '%s'" % cactusDiskDatabaseString,
                            "--logLevel", logLevel])
    logger.info("Ran cactus_phylogeny okay")
    
def runCactusAdjacencies(cactusDiskDatabaseString, flowerNames=encodeFlowerNames((0,)), logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_fillAdjacencies",
                            "--cactusDisk '%s'" % cactusDiskDatabaseString,
                            "--logLevel", logLevel])
    logger.info("Ran cactus_fillAdjacencies OK")

def runCactusConvertAlignmentToCactus(cactusDiskDatabaseString, constraintsFile, newConstraintsFile, logLevel=None):
    """Takes a cigar file and makes an equivalent cigar file using the internal coordinate system format of cactus.
    """
    logLevel = getLogLevelString2(logLevel)
    cactus_call(parameters=["cactus_workflow_convertAlignmentCoordinates",
                            logLevel, "'%s'" % cactusDiskDatabaseString,
                            constraintsFile, newConstraintsFile])

def runCactusFlowerStats(cactusDiskDatabaseString, flowerName, logLevel=None):
    """Prints stats for the given flower
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStatsString = cactus_call(check_output=True,
                                    parameters=["cactus_workflow_flowerStats",
                                                logLevel, "'%s'" % cactusDiskDatabaseString, flowerName])
    return flowerStatsString.split("\n")[0]

def runCactusMakeNormal(cactusDiskDatabaseString, flowerNames, maxNumberOfChains=0, logLevel=None):
    """Makes the given flowers normal (see normalisation for the various phases)
    """
    logLevel = getLogLevelString2(logLevel)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_normalisation",
                            "--cactusDisk '%s'" % cactusDiskDatabaseString,
                            "--maxNumberOfChains", maxNumberOfChains,
                            "--logLevel", logLevel])

def runCactusBar(cactusDiskDatabaseString, flowerNames, logLevel=None,
                         spanningTrees=None, maximumLength=None, 
                         gapGamma=None,
                         matchGamma=None,
                         splitMatrixBiggerThanThis=None,
                         anchorMatrixBiggerThanThis=None,
                         repeatMaskMatrixBiggerThanThis=None,
                         diagonalExpansion=None,
                         constraintDiagonalTrim=None,
                         minimumBlockDegree=None,
                         minimumIngroupDegree=None,
                         minimumOutgroupDegree=None,
                         alignAmbiguityCharacters=None,
                         pruneOutStubAlignments=None,
                         useProgressiveMerging=None,
                         calculateWhichEndsToComputeSeparately=None,
                         largeEndSize=None,
                         endAlignmentsToPrecomputeOutputFile=None,
                         precomputedAlignments=None,
                         ingroupCoverageFile=None,
                         minimumSizeToRescue=None,
                         minimumCoverageToRescue=None,
                         minimumNumberOfSpecies=None,
                         jobName=None,
                         fileStore=None,
                         features=None):
    """Runs cactus base aligner."""
    logLevel = getLogLevelString2(logLevel)
    maximumLength = nameValue("maximumLength", maximumLength, int)
    spanningTrees = nameValue("spanningTrees", spanningTrees, int)
    gapGamma = nameValue("gapGamma", gapGamma, float)
    matchGamma = nameValue("matchGamma", matchGamma, float)
    splitMatrixBiggerThanThis=nameValue("splitMatrixBiggerThanThis", splitMatrixBiggerThanThis, int)
    anchorMatrixBiggerThanThis=nameValue("anchorMatrixBiggerThanThis", anchorMatrixBiggerThanThis, int)
    repeatMaskMatrixBiggerThanThis=nameValue("repeatMaskMatrixBiggerThanThis", repeatMaskMatrixBiggerThanThis, int)
    diagonalExpansion=nameValue("diagonalExpansion", diagonalExpansion, int)
    constraintDiagonalTrim = nameValue("constraintDiagonalTrim", constraintDiagonalTrim, int)
    minimumBlockDegree = nameValue("minimumDegree", minimumBlockDegree, int)
    minimumIngroupDegree = nameValue("minimumIngroupDegree", minimumIngroupDegree, int)
    minimumOutgroupDegree = nameValue("minimumOutgroupDegree", minimumOutgroupDegree, int)
    pruneOutStubAlignments = nameValue("pruneOutStubAlignments", pruneOutStubAlignments, bool)
    alignAmbiguityCharacters = nameValue("alignAmbiguityCharacters", alignAmbiguityCharacters, bool)
    useProgressiveMerging=nameValue("useProgressiveMerging", useProgressiveMerging, bool)
    calculateWhichEndsToComputeSeparately=nameValue("calculateWhichEndsToComputeSeparately", calculateWhichEndsToComputeSeparately, bool)
    largeEndSize=nameValue("largeEndSize", largeEndSize, int)
    if endAlignmentsToPrecomputeOutputFile is not None:
        endAlignmentsToPrecomputeOutputFile = os.path.basename(endAlignmentsToPrecomputeOutputFile)
    endAlignmentsToPrecomputeOutputFile=nameValue("endAlignmentsToPrecomputeOutputFile", endAlignmentsToPrecomputeOutputFile, str)
    if precomputedAlignments is not None:
        precomputedAlignments = map(os.path.basename, precomputedAlignments)
        precomputedAlignments = " ".join(precomputedAlignments)
    precomputedAlignments=nameValue("precomputedAlignments", precomputedAlignments, str, quotes=True)
    ingroupCoverageFile = nameValue("ingroupCoverageFile", ingroupCoverageFile, str, quotes=True)
    minimumSizeToRescue = nameValue("minimumSizeToRescue", minimumSizeToRescue, int)
    minimumCoverageToRescue = nameValue("minimumCoverageToRescue", minimumCoverageToRescue, float)
    minimumNumberOfSpecies = nameValue("minimumNumberOfSpecies", minimumNumberOfSpecies, int)

    masterMessages = cactus_call(stdin_string=flowerNames, check_output=True,
                                 option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                                 parameters=["cactus_bar",
                                             "--logLevel", logLevel, spanningTrees, maximumLength,
                                             gapGamma, matchGamma,
                                             splitMatrixBiggerThanThis, anchorMatrixBiggerThanThis,
                                             repeatMaskMatrixBiggerThanThis,
                                             constraintDiagonalTrim, minimumBlockDegree, minimumIngroupDegree,
                                             minimumOutgroupDegree,  
                                             alignAmbiguityCharacters, pruneOutStubAlignments,
                                             diagonalExpansion,
                                             useProgressiveMerging, calculateWhichEndsToComputeSeparately,
                                             largeEndSize, endAlignmentsToPrecomputeOutputFile,
                                             precomputedAlignments, ingroupCoverageFile,
                                             minimumSizeToRescue, minimumCoverageToRescue,
                                             minimumNumberOfSpecies],
                                 job_name=jobName, fileStore=fileStore, features=features)
        
    logger.info("Ran cactus_bar okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]

def runCactusSecondaryDatabase(secondaryDatabaseString, create=True):
    cactus_call(parameters=["cactus_secondaryDatabase",
                secondaryDatabaseString, create])
            
def runCactusReference(cactusDiskDatabaseString, flowerNames, logLevel=None,
                       jobName=None, features=None, fileStore=None,
                       matchingAlgorithm=None, 
                       referenceEventString=None, 
                       permutations=None,
                       useSimulatedAnnealing=None,
                       theta=None,
                       phi=None, 
                       maxWalkForCalculatingZ=None,
                       ignoreUnalignedGaps=None,
                       wiggle=None, 
                       numberOfNs=None,
                       minNumberOfSequencesToSupportAdjacency=None,
                       makeScaffolds=None):
    """Runs cactus reference."""
    logLevel = getLogLevelString2(logLevel)
    matchingAlgorithm = nameValue("matchingAlgorithm", matchingAlgorithm)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    permutations = nameValue("permutations", permutations, int)
    useSimulatedAnnealing = nameValue("useSimulatedAnnealing", useSimulatedAnnealing, bool)
    theta = nameValue("theta", theta, float)
    phi = nameValue("phi", phi, float)
    maxWalkForCalculatingZ = nameValue("maxWalkForCalculatingZ", maxWalkForCalculatingZ, int)
    ignoreUnalignedGaps = nameValue("ignoreUnalignedGaps", ignoreUnalignedGaps, bool)
    wiggle = nameValue("wiggle", wiggle, float)
    numberOfNs = nameValue("numberOfNs", numberOfNs, int)
    minNumberOfSequencesToSupportAdjacency = nameValue("minNumberOfSequencesToSupportAdjacency", minNumberOfSequencesToSupportAdjacency, int)
    makeScaffolds = nameValue("makeScaffolds", makeScaffolds, bool)

    masterMessages = cactus_call(stdin_string=flowerNames, check_output=True,
                option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                parameters=["cactus_reference",
                            "--logLevel", logLevel,
                            matchingAlgorithm, referenceEventString, permutations, 
                            useSimulatedAnnealing, theta, phi, maxWalkForCalculatingZ, ignoreUnalignedGaps,
                            wiggle, numberOfNs, minNumberOfSequencesToSupportAdjacency, makeScaffolds],
                job_name=jobName,
                features=features,
                fileStore=fileStore)
    logger.info("Ran cactus_reference okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]
    
def runCactusAddReferenceCoordinates(cactusDiskDatabaseString, flowerNames,
                                     jobName=None, fileStore=None, features=None,
                                     logLevel=None, referenceEventString=None,
                                     outgroupEventString=None, secondaryDatabaseString=None,
                                     bottomUpPhase=None):
    logLevel = getLogLevelString2(logLevel)
    bottomUpPhase = nameValue("bottomUpPhase", bottomUpPhase, bool)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    outgroupEventString = nameValue("outgroupEventString", outgroupEventString)
    secondaryDatabaseString = nameValue("secondaryDisk", secondaryDatabaseString, quotes=True)
    cactus_call(stdin_string=flowerNames,
                option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                parameters=["cactus_addReferenceCoordinates",
                            secondaryDatabaseString,
                            "--logLevel", logLevel,
                            referenceEventString,
                            outgroupEventString,
                            bottomUpPhase],
                job_name=jobName,
                features=features,
                fileStore=fileStore)

def runCactusCheck(cactusDiskDatabaseString, 
                   flowerNames=encodeFlowerNames((0,)), 
                   logLevel=None, 
                   recursive=None,
                   checkNormalised=None):
    logLevel = getLogLevelString2(logLevel)
    recursive = nameValue("recursive", recursive, bool)
    checkNormalised = nameValue("checkNormalised", checkNormalised, bool)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_check",
                            "--cactusDisk '%s'" % cactusDiskDatabaseString,
                            "--logLevel", logLevel,
                            recursive, checkNormalised])
    logger.info("Ran cactus check")
    
def _fn(toilDir, 
      logLevel=None, retryCount=0, 
      batchSystem="single_machine", 
      rescueJobFrequency=None,
      skipAlignments=False,
      buildAvgs=False, buildReference=False,
      buildHal=False,
      buildFasta=False,
      toilStats=False,
      maxThreads=None,
      maxCpus=None,
      defaultMemory=None,
      logFile=None,
      extraToilArgumentsString=""):
    logLevel = getLogLevelString2(logLevel)
    skipAlignments = nameValue("skipAlignments", skipAlignments, bool)
    buildAvgs = nameValue("buildAvgs", buildAvgs, bool)
    buildReference = nameValue("buildReference", buildReference, bool)
    buildHal = nameValue("buildHal", buildHal, bool)
    buildFasta = nameValue("buildFasta", buildFasta, bool)
    #Jobtree args
    batchSystem = nameValue("batchSystem", batchSystem, str)
    retryCount = nameValue("retryCount", retryCount, int)
    rescueJobFrequency = nameValue("rescueJobsFrequency", rescueJobFrequency, int)
    toilStats = nameValue("stats", toilStats, bool)
    maxThreads = nameValue("maxThreads", maxThreads, int)
    maxCpus = nameValue("maxCpus", maxCpus, int)
    defaultMemory= nameValue("defaultMemory", defaultMemory, int)
    logFile = nameValue("logFile", logFile, str)
    return "%s %s %s %s --logLevel %s %s %s %s %s %s %s %s %s %s %s %s" % (toilDir, skipAlignments, buildAvgs, 
             buildReference, logLevel, buildHal, buildFasta, batchSystem, retryCount, rescueJobFrequency, toilStats, maxThreads, maxCpus, logFile, defaultMemory, extraToilArgumentsString)
     
def runCactusWorkflow(experimentFile,
                      toilDir, 
                      logLevel=None, retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      skipAlignments=False,
                      buildAvgs=False, buildReference=False,
                      buildHal=False,
                      buildFasta=False,
                      toilStats=False,
                      maxThreads=None,
                      maxCpus=None,
                      defaultMemory=None,
                      logFile=None,
                      extraToilArgumentsString=""):
    arguments = ("--experiment %s" % experimentFile) + " " + _fn(toilDir,
                      logLevel, retryCount, batchSystem, rescueJobFrequency, skipAlignments,
                      buildAvgs, buildReference, buildHal, buildFasta, toilStats, maxThreads, maxCpus, defaultMemory, logFile, extraToilArgumentsString=extraToilArgumentsString)

    import cactus.pipeline.cactus_workflow as cactus_workflow
    cactus_workflow.runCactusWorkflow(arguments.split())
    logger.info("Ran the cactus workflow okay")
    
def runCactusCreateMultiCactusProject(experimentFile, outputDir, 
                                      logLevel=None, fixNames=True,
                                      root=None):
    logLevel = getLogLevelString2(logLevel)
    root = nameValue("root", root, str, quotes=True)
    command = "cactus_createMultiCactusProject.py %s %s --fixNames=%s %s" % (experimentFile, outputDir, str(fixNames), root)
    system(command)
    logger.info("Ran the cactus create multi project")
    
def runCactusProgressive(inputDir,
                      toilDir, 
                      logLevel=None, retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      skipAlignments=False,
                      buildHal=None,
                      buildFasta=None,
                      buildAvgs=False, 
                      toilStats=False,
                      maxThreads=None,
                      maxCpus=None,
                      defaultMemory=None,
                      recursive=None,
                      logFile=None,
                      event=None,
                      extraToilArgumentsString="",
                      profileFile=None):
    command = ("cactus_progressive.py --project %s" % inputDir) + " " + _fn(toilDir, 
                      logLevel, retryCount, batchSystem, rescueJobFrequency, skipAlignments,
                      buildAvgs, None,
                      buildHal,
                      buildFasta,
                      toilStats, maxThreads, maxCpus, defaultMemory, logFile, extraToilArgumentsString=extraToilArgumentsString) + \
                      (" %s %s" % (nameValue("recursive", recursive, bool),
                                      nameValue("event", event)))
    if profileFile != None:
        command = "python -m cProfile -o %s %s/bin/%s" % (profileFile, cactusRootPath(), command)
    system(command)                   
    logger.info("Ran the cactus progressive okay")
    
def runCactusHalGenerator(cactusDiskDatabaseString,
                          secondaryDatabaseString, 
                          flowerNames,
                          referenceEventString, 
                          outputFile=None,
                          showOnlySubstitutionsWithRespectToReference=None,
                          logLevel=None,
                          jobName=None,
                          features=None,
                          fileStore=None):
    logLevel = getLogLevelString2(logLevel)
    if outputFile:
        outputFile = os.path.basename(outputFile)
    cactus_call(stdin_string=flowerNames,
                option_string="--cactusDisk '%s' --secondaryDisk '%s'" % (cactusDiskDatabaseString, secondaryDatabaseString),
                parameters=["cactus_halGenerator",
                            "--logLevel", logLevel,
                            nameValue("referenceEventString", referenceEventString),
                            nameValue("outputFile", outputFile),
                            nameValue("showOnlySubstitutionsWithRespectToReference",
                                      showOnlySubstitutionsWithRespectToReference, bool)],
                job_name=jobName, features=features, fileStore=fileStore)
                            
def runCactusFastaGenerator(cactusDiskDatabaseString,
                            flowerName,
                            outputFile,
                            referenceEventString=None, 
                            logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    cactus_call(option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                parameters=["cactus_fastaGenerator",
                            "--flowerName", flowerName,
                            "--outputFile", outputFile,
                            "--logLevel", logLevel,
                            nameValue("referenceEventString", referenceEventString)])
    
def runCactusAnalyseAssembly(sequenceFile):
    return cactus_call(check_output=True,
                parameters=["cactus_analyseAssembly",
                            sequenceFile])[:-1]
    
def runToilStats(toil, outputFile):
    system("toil stats %s --outputFile %s" % (toil, outputFile))
    logger.info("Ran the job-tree stats command apparently okay")
def runToilStatusAndFailIfNotComplete(toilDir):
    command = "toil status %s --failIfNotComplete --verbose" % toilDir
    system(command)

def runLastz(seq1, seq2, alignmentsFile, lastzArguments, work_dir=None, samplingRates=None):
    #Have to specify the work_dir manually for this, since
    #we're adding arguments to the filename
    assert os.path.dirname(seq1) == os.path.dirname(seq2)
    work_dir = os.path.dirname(seq1)
    parameters = ["cPecanLastz",
                            "--format=cigar",
                            "--notrivial",
                            lastzArguments,
                            "%s[multiple][nameparse=darkspace]" % os.path.basename(seq1),
                            "%s[nameparse=darkspace]" % os.path.basename(seq2)]
    if samplingRates:
        parameters += ["--samplingRates", samplingRates]
    cactus_call(work_dir=work_dir, outfile=alignmentsFile,
                parameters=parameters)

def runSelfLastz(seq, alignmentsFile, lastzArguments, work_dir=None, samplingRates=None):
    work_dir = os.path.dirname(seq)
    parameters=["cPecanLastz",
                            "--format=cigar",
                            "--notrivial",
                            lastzArguments,
                            "%s[multiple][nameparse=darkspace]" % os.path.basename(seq),
                            "%s[nameparse=darkspace]" % os.path.basename(seq)]
    if samplingRates:
        parameters += ["--samplingRates", samplingRates]
    cactus_call(work_dir=work_dir, outfile=alignmentsFile,
                parameters=parameters)
    
def runCactusRealign(seq1, seq2, inputAlignmentsFile, outputAlignmentsFile, realignArguments, work_dir=None):
    cactus_call(infile=inputAlignmentsFile, outfile=outputAlignmentsFile, work_dir=work_dir,
                parameters=["cPecanRealign", realignArguments, seq1, seq2])

def runCactusSelfRealign(seq, inputAlignmentsFile, outputAlignmentsFile, realignArguments, work_dir=None):
    cactus_call(infile=inputAlignmentsFile, outfile=outputAlignmentsFile, work_dir=work_dir,
                parameters=["cPecanRealign", realignArguments, seq])

def runCactusCoverage(sequenceFile, alignmentsFile, work_dir=None):
    return cactus_call(check_output=True, work_dir=work_dir,
                parameters=["cactus_coverage", sequenceFile, alignmentsFile])

def runGetChunks(sequenceFiles, chunksDir, chunkSize, overlapSize, work_dir=None):
    return [chunk for chunk in cactus_call(work_dir=work_dir,
                                           check_output=True,
                                           parameters=["cactus_blast_chunkSequences",
                                           getLogLevelString(),
                                           chunkSize,
                                           overlapSize,
                                           chunksDir] + sequenceFiles).split("\n") if chunk != ""]

def pullCactusImage():
    """Ensure that the cactus Docker image is pulled."""
    if os.environ.get('CACTUS_DOCKER_MODE') == "0":
        return
    dockerOrg = getDockerOrg()
    dockerTag = getDockerTag()
    image = "%s/cactus:%s" % (dockerOrg, dockerTag)
    call = ["docker", "pull", image]
    process = subprocess32.Popen(call, stdout=subprocess32.PIPE,
                                 stderr=sys.stderr, bufsize=-1)
    output, _ = process.communicate()
    if process.returncode != 0:
        raise RuntimeError("Command %s failed with output: %s" % (call, output))

def getDockerOrg():
    """Get where we should find the cactus containers."""
    if "CACTUS_DOCKER_ORG" in os.environ:
        return os.environ["CACTUS_DOCKER_ORG"]
    else:
        return "quay.io/comparative-genomics-toolkit"

def getDockerTag():
    """Get what docker tag we should use for the cactus image
    (either forced to be latest or the current cactus commit)."""
    if 'CACTUS_USE_LATEST' in os.environ:
        return "latest"
    else:
        return cactus_commit

def maxMemUsageOfContainer(containerInfo):
    """Return the max RSS usage (in bytes) of a container, or None if something failed."""
    if containerInfo['id'] is None:
        # Try to get the internal container ID from the docker name
        try:
            id = popenCatch("docker inspect -f '{{.Id}}' %s" % containerInfo['name']).strip()
            containerInfo['id'] = id
        except:
            # Not yet running
            return None
    # Try to check for the maximum memory usage ever used by that
    # container, in a few different possible locations depending on
    # the distribution
    possibleLocations = ["/sys/fs/cgroup/memory/docker/%s/memory.max_usage_in_bytes",
                         "/sys/fs/cgroup/memory/system.slice.docker-%s.scope/memory.max_usage_in_bytes"]
    possibleLocations = [s % containerInfo['id'] for s in possibleLocations]
    for location in possibleLocations:
        try:
            with open(location) as f:
                return int(f.read())
        except IOError:
            # Not at this location, or sysfs isn't mounted
            continue
    return None

#TODO: This function is a mess
def cactus_call(tool=None,
                work_dir=None,
                parameters=None,
                rm=True,
                detached=True,
                check_output=False,
                container_name=None,
                mounts=None,
                infile=None,
                outfile=None,
                stdin_string=None,
                option_string="",
                server=False,
                shell=True,
                port=None,
                check_result=False,
                dockstore=None,
                job_name=None,
                features=None,
                fileStore=None):
    if dockstore is None:
        dockstore = getDockerOrg()
    if parameters is None:
        parameters = []

    def moveToWorkDir(work_dir, arg):
        if isinstance(arg, str) and os.path.isfile(arg):
            if not os.path.dirname(arg) == work_dir:
                _log.info('Copying file %s to work dir' % arg)
                shutil.copy(arg, work_dir)

    if work_dir:
        for arg in parameters:
            moveToWorkDir(work_dir, arg)

    parameters = [str(par) for par in parameters]
    if not work_dir:
    #Make sure all the paths we're accessing are in the same directory
        files = [par for par in parameters if os.path.isfile(par)]
        folders = [par for par in parameters if os.path.isdir(par)]
        work_dirs = set([os.path.dirname(fileName) for fileName in files] + [os.path.dirname(folder) for folder in folders])
        _log.info("Work dirs: %s" % work_dirs)
        if len(work_dirs) > 1:
            work_dir = os.path.commonprefix(work_dirs)
        elif len(work_dirs) == 1:
            work_dir = work_dirs.pop()

    #If there are no input files, or if their MRCA is '' (when working
    #with relative paths), just set the current directory as the work
    #dir
    if work_dir is None or work_dir == '':
        work_dir = "."
    _log.info("Docker work dir: %s" % work_dir)

    #We'll mount the work_dir containing the paths as /data in the container,
    #so set all the paths to their basenames. The container will access them at
    #/data/<path>
    def adjustPath(path, wd):
        # Hack to relativize paths that are not provided as a
        # single argument (i.e. multiple paths that are
        # space-separated and quoted)
        if wd != '.':
            if not wd.endswith('/'):
                wd = wd + '/'
            return path.replace(wd, '')
        else:
            return path

    if work_dir and os.environ.get('CACTUS_DOCKER_MODE') != "0":
        parameters = [adjustPath(par, work_dir) for par in parameters]

    base_docker_call = ['docker', 'run',
                        '--interactive',
                        '--net=host',
                        '--log-driver=none',
                        '-u', '%s:%s' % (os.getuid(), os.getgid()),
                        '-e', 'ST_ABORT=1',
                        '-e', 'ST_ABORT_UNCAUGHT=1',
                        '-v', '{}:/data'.format(os.path.abspath(work_dir))]

    if port:
        base_docker_call += ["-p %d:%d" % (port, port)]

    containerInfo = { 'name': str(uuid.uuid4()), 'id': None }
    base_docker_call.extend(['--name', containerInfo['name']])
    if rm:
        base_docker_call.append('--rm')


    parameters = [par for par in parameters if par != '']

    parameters = " ".join(parameters)

    if not tool:
        tool = "cactus"

    docker_tag = getDockerTag()

    if os.environ.get('CACTUS_DOCKER_MODE') == "0":
        _log.info("Calling tool from local cactus installation.")
        call = parameters
    else:
        tool = "%s/%s:%s" % (dockstore, tool, docker_tag)
        call = " ".join(base_docker_call) + " " + tool + " " + parameters
    if option_string:
        call += " " + option_string
    

    if stdin_string:
        _log.info("Input string: %s" % stdin_string)


    stdinFileHandle = None
    stdoutFileHandle = None
    if stdin_string:
        stdinFileHandle = subprocess32.PIPE
    elif infile:
        stdinFileHandle = open(infile, 'r')
    if outfile:
        stdoutFileHandle = open(outfile, 'w')
    if check_output:
        stdoutFileHandle = subprocess32.PIPE


    _log.info("Running the command %s" % call)
    if not shell:
        call = call.split()
    process = subprocess32.Popen(call, shell=shell,
                                 stdin=stdinFileHandle, stdout=stdoutFileHandle,
                                 stderr=sys.stderr, bufsize=-1)

    if server:
        return process

    memUsage = 0
    first_run = True
    while True:
        try:
            # Wait a bit to see if the process is done
            output, nothing = process.communicate(stdin_string if first_run else None, timeout=10)
        except subprocess32.TimeoutExpired:
            # Every so often, check the memory usage of the container
            updatedMemUsage = maxMemUsageOfContainer(containerInfo)
            if updatedMemUsage is not None:
                assert memUsage <= updatedMemUsage, "memory.max_usage_in_bytes should never decrease"
                memUsage = updatedMemUsage
            first_run = False
        else:
            break
    _log.info("Used %s max memory" % memUsage)
    if job_name is not None and features is not None and fileStore is not None:
        # Log a datapoint for the memory usage for these features.
        fileStore.logToMaster("Max memory used for job %s (tool %s) "
                              "on JSON features %s: %s" % (job_name, parameters.split()[0],
                                                           json.dumps(features), memUsage))
    if check_result:
        return process.returncode

    if process.returncode != 0:
        raise RuntimeError("Command %s failed with output: %s" % (call, output))

    if check_output:
        return output

class RunAsFollowOn(Job):
    def __init__(self, job, *args, **kwargs):
        Job.__init__(self, memory=100000000, preemptable=True)
        self._args = args
        self._kwargs = kwargs
        self.job = job
    def run(self, fileStore):
        return self.addFollowOn(self.job(*self._args, **self._kwargs)).rv()
        
class RoundedJob(Job):
    """Thin wrapper around Toil.Job to round up resource requirements.

    Rounding is useful to make Toil's Mesos scheduler more
    efficient--it runs a process that is O(n log n) in the number of
    different resource requirements for every offer received, so
    thousands of slightly different requirements will slow down the
    leader and the workflow.
    """
    # Default rounding amount: 100 MiB
    roundingAmount = 100*1024*1024
    def __init__(self, memory=None, cores=None, disk=None, preemptable=None,
                 unitName=None, checkpoint=False):
        if memory is not None:
            memory = self.roundUp(memory)
        if disk is not None:
            disk = self.roundUp(disk)
        super(RoundedJob, self).__init__(memory=memory, cores=cores, disk=disk,
                                         preemptable=preemptable, unitName=unitName,
                                         checkpoint=checkpoint)

    def roundUp(self, bytesRequirement):
        """
        Round the amount up to the next self.roundingAmount.

        >>> j = RoundedJob()
        >>> j.roundingAmount = 100000000
        >>> j.roundUp(1000)
        10000000
        >>> j.roundUp(200000000)
        200000000
        >>> j.roundUp(200000001)
        300000000
        """
        if bytesRequirement % self.roundingAmount == 0:
            return bytesRequirement
        return (bytesRequirement // self.roundingAmount + 1) * self.roundingAmount

def readGlobalFileWithoutCache(fileStore, jobStoreID):
    """Reads a jobStoreID into a file and returns it, without touching
    the cache.

    Works around toil issue #1532.
    """
    f = fileStore.getLocalTempFile()
    fileStore.jobStore.readFile(jobStoreID, f)
    return f

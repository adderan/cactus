#!/usr/bin/env python

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" This file contains the basic Toil pattern for launching ktservers.
Given a target T we want to effectively do the following:
1) spawn ktserver
2) run T (along with all its children and follow-ons)
3) kill ktserver

Given an input Root target, this is accomplished by means
of creating three new targets:

Launch, Block, and Kill:

Root creates a global temporary file and writes it to the disk.
Root then spawns two child targets: Launch and Block

Launch launches a ktserver on whatever system it is executed on.
Launch then monitors the temporary file.  When this file is
deleted, the server is terminated.

Block polls the temporary (and server log) file indefinitely.  As
soon as it can verify that the server is running, it schedules
T as a child.
Block also schedules Kill as a follow-on

Kill deletes the temporary file, causing the ktserver to terminate.  
"""

import os
import sys
import random
import math
import xml.etree.ElementTree as ET

from sonLib.bioio import getTempFile
from toil.job import Job
from cactus.shared.experimentWrapper import DbElemWrapper
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.pipeline.ktserverControl import runKtserver
from cactus.pipeline.ktserverControl import blockUntilKtserverIsRunnning
from cactus.pipeline.ktserverControl import killKtServer
from cactus.pipeline.ktserverControl import getKtServerReport

###############################################################################
# Launch childTargetClosure (which is a Toil Job with bound parameters)
# as a child of rootTarget, using the ktserver launch pattern described
# above to make sure that a ktserver is running somewhere for the lifespan
# of childTargetClosure using the pattern described above.
#
# rootJob : existing Toil job which will serve as the root
#              of all compuation (which will be added as children)
# rootJobFileStore ; the FileStore of the rootJob, which is used
# to create the kill switch file
# newChild : the child PhaseJob we wish to execute while the
#            ktserver is running.  if isSecondary is true then
#            newChild is a recursion target and not a phase target
# maxMemory : memory (in bytes) for toil to request for ktserver
# maxCpu : number of cpus for jobTree to request for ktserver
# isSecondary : flag whether or not the database is secondary.  If False,
#               then the port and host are written to the regular conf node,
#               if True, then we update the secondaryString (all this within
#               the phase's workflow options)
# createTimeout : seconds to wait for server to be valid after running ktserver
#                 on a new database before giving up
# loadTimeout : seconds to wait for server to be valid after running ktserver
#                 on an existing database before giving up
# blockTimeout : seconds that Block target waits for the ktserver to be
#                launched before aborting with an error
# blockTimestep : length of interval in seconds between polling for the server
# runTimeout : maximum number of seconds for the ktserver to run
# runTimestep : polling interval for above
# killTimeout : amount of time to wait for server to die after deleting
#               the kill switch file before throwing an error
###############################################################################
def addKtserverDependentChild(rootJob, rootJobFileStore, newChild, maxMemory, maxCpu,
                              isSecondary = False,
                              createTimeout = 30, loadTimeout = 10000,
                              blockTimeout=sys.maxint, blockTimestep=10,
                              runTimeout=sys.maxint, runTimestep=10,
                              killTimeout=10000):
    from cactus.pipeline.cactus_workflow import CactusPhasesJob
    from cactus.pipeline.cactus_workflow import CactusRecursionJob
    
    assert isinstance(rootJob, Job)

    if killTimeout < runTimestep * 2:
        killTimeout = runTimestep * 2
    
    killSwitchPath = getTempFile(suffix="_kill.txt",
                                 rootDir=rootJobFileStore.getLocalTempDir())
    killSwitchFile = open(killSwitchPath, "w")
    killSwitchFile.write("init")
    killSwitchFile.close()
    killSwitchFileID = rootJobFileStore.writeGlobalFile(killSwitchPath)

    if isSecondary == False:
        assert isinstance(newChild, CactusPhasesJob)
        wfArgs = newChild.cactusWorkflowArguments
        dbElem = ExperimentWrapper(wfArgs.experimentNode)
    else:
        assert isinstance(newChild, CactusRecursionJob)
        dbString = newChild.getOptionalPhaseAttrib("secondaryDatabaseString")
        assert dbString is not None
        confXML = ET.fromstring(dbString)
        dbElem = DbElemWrapper(confXML)
    
    rootJob.addChild(
        KtserverJobLauncher(dbElem, killSwitchFileID, maxMemory,
                               maxCpu, createTimeout,
                               loadTimeout, runTimeout, runTimestep))
    rootJob.addChild(
        KtserverJobBlocker(killSwitchFileID, newChild, isSecondary,
                                blockTimeout, blockTimestep, killTimeout))


###############################################################################
# Launch the server on whatever node runs this target
###############################################################################
class KtserverJobLauncher(Job):
    def __init__(self, dbElem, killSwitchFileID,
                 maxMemory, maxCpu, createTimeout,
                 loadTimeout, runTimeout, runTimestep):
        Job.__init__(self, memory=maxMemory, cores=maxCpu)
        self.dbElem = dbElem
        self.killSwitchFileID = killSwitchFileID
        self.createTimeout = createTimeout
        self.loadTimeout = loadTimeout
        self.runTimeout = runTimeout
        self.runTimestep = runTimestep
        
    def run(self, fileStore):
        fileStore.logToMaster("Launching ktserver %s with killPath %s" % (
            ET.tostring(self.dbElem.getDbElem()), self.killSwitchFileID))
        runKtserver(self.dbElem, self.killSwitchFileID, fileStore,
                    maxPortsToTry=100, readOnly = False,
                    createTimeout=self.createTimeout,
                    loadTimeout=self.loadTimeout,
                    killTimeout=self.runTimeout,
                    killPingInterval=self.runTimestep)

###############################################################################
# Block until the server's detected.
# Run the child job as a child, and kill the server in a follow-on
###############################################################################
class KtserverJobBlocker(Job):
    def __init__(self, killSwitchFileID, newChild, isSecondary,
                 blockTimeout, blockTimestep, killTimeout):
        Job.__init__(self)
        self.killSwitchFileID = killSwitchFileID
        self.newChild = newChild
        self.isSecondary = isSecondary
        self.blockTimeout = blockTimeout
        self.blockTimestep = blockTimestep
        self.killTimeout = killTimeout
        
    def run(self, fileStore):
        if self.isSecondary == False:
            wfArgs = self.newChild.cactusWorkflowArguments
            experiment = ExperimentWrapper(wfArgs.experimentNode)
            dbElem = experiment
        else:
            dbString = self.newChild.getOptionalPhaseAttrib(
                "secondaryDatabaseString")
            assert dbString is not None
            confXML = ET.fromstring(dbString)
            dbElem = DbElemWrapper(confXML)

        fileStore.logToMaster("Blocking on ktserver %s with killSwitchFileID %s" % (
            ET.tostring(dbElem.getDbElem()), self.killSwitchFileID))
            
        blockUntilKtserverIsRunnning(dbElem, self.killSwitchFileID,
                                     self.blockTimeout, self.blockTimestep)

        if self.isSecondary == False:
            experiment.writeXML(wfArgs.experimentFile)
            wfArgs.cactusDiskDatabaseString = dbElem.getConfString()
        else:
            self.newChild.phaseNode.attrib[
                "secondaryDatabaseString"] = dbElem.getConfString()
            # added on as a hack to get this into the experiment.xml
            etPath = self.newChild.phaseNode.attrib[
                "experimentPath"]
            experiment = ExperimentWrapper(ET.parse(etPath).getroot())
            experiment.setSecondaryDBElem(dbElem)
            experiment.writeXML(etPath)            
        
        self.addChild(self.newChild)
        self.addFollowOn(KtserverJobKiller(dbElem,
                                                    self.killSwitchFileID,
                                                    self.killTimeout))

###############################################################################
# Kill the server by deleting its kill switch file
###############################################################################
class KtserverJobKiller(Job):
    def __init__(self, dbElem, killSwitchFileID, killTimeout):
        Job.__init__(self)
        self.dbElem = dbElem
        self.killSwitchFileID = killSwitchFileID
        self.killTimeout = killTimeout
        
    def run(self, fileStore):
        fileStore.logToMaster("Killing ktserver %s with killSwitchFileID %s" % (
            ET.tostring(self.dbElem.getDbElem()), self.killSwitchFileID))
        report = getKtServerReport(self.dbElem)
        fileStore.logToMaster(report)
        killKtServer(self.dbElem, killSwitchFileID, fileStore,
                     killTimeout=self.killTimeout)
    

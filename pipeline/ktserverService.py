from toil.job import Job
import sys
import subprocess
import os
import socket
from toil.lib.bioio import logger
import threading
from cactus.pipeline.ktserverControl import validateKtserver, isKtServerRunning, isKtServerReorganizing, isKtServerOnTakenPort, getLogPath, getHostName, getKtserverCommand, writeStatusToSwitchFile

def addKtserverDependentChild():
    pass

class KtserverService(Job.Service):
    def __init__(self, dbElem, statusPath,
                              createTimeout = 30, loadTimeout = 10000,
                              blockTimeout=sys.maxint, blockTimestep=10,
                              runTimeout=sys.maxint, runTimestep=10,
                              killTimeout=10000):
        Job.Service.__init__(self)
        self.dbElem = dbElem
        self.dbPathExists = False
        self.loadTimeout = loadTimeout
        self.blockTimeout = blockTimeout
        self.runTimeout = runTimeout
        self.readOnly = False
        self.statusPath = statusPath

        if self.dbElem.getDbInMemory() == False:
            if not os.path.splitext(self.dbElem.getDbName())[1] == ".kch":
                raise RuntimeError("Expected path to end in .kch: %s" %
                        self.dbElem.getDbName())
            self.dbPathExists = os.path.exists(self.dbPath)
        self.logPath = getLogPath(self.dbElem)
        self.logDir = os.path.dirname(self.logPath)
        if os.path.exists(self.logDir) is False:
            os.makedirs(self.logDir)
        assert os.path.isdir(self.logDir)
        self.basePort = self.dbElem.getDbPort()
        self.maxPortsToTry = 100
        host = self.dbElem.getDbHost()
        if host is not None:
            self.dbElem.setDbHost(host)
        else:
            self.dbElem.setDbHost(getHostName())


    def start(self):
        statusFile = open(self.statusPath, "w")
        statusFile.write("init")
        statusFile.close()

        success = False
        self.process = None
        procWaiter = None
        logger.info("Trying ports in range: %s to %s" % (self.basePort, self.basePort + self.maxPortsToTry))
        for port in xrange(self.basePort, self.basePort + self.maxPortsToTry):
            self.dbElem.setDbPort(port)
            logger.info("Trying port %s" % self.dbElem.getDbPort())
            if os.path.exists(self.logPath):
                os.remove(self.logPath)
            if isKtServerOnTakenPort(self.dbElem, pretest=True):
                logger.info("Port taken")
                continue
            cmd = getKtserverCommand(self.dbElem, self.dbPathExists, self.readOnly)
            logger.info("Executing command %s" % cmd.split())
            self.process = subprocess.Popen(cmd.split(), shell=False,
                    stdout=subprocess.PIPE,
                    stderr=sys.stderr, bufsize=-1)
            assert self.process is not None
            procWaiter = ProcessWaiter(self.process)
            procWaiter.start()
            writeStatusToSwitchFile(self.dbElem, self.process.pid, self.statusPath)


            success = validateKtserver(self.process, self.dbElem, self.statusPath, self.createTimeout, self.loadTimeout)
            if success is True:
                logger.info("Launched ktserver")
                break
            else:
                if self.process.returncode is None:
                    logger.info("Killing process")
                    self.process.kill()
                self.process = None
        if success is False:
            raise RuntimeError("Unable to launch ktserver. "+
                    "Server log is: %s" % logPath)
    def stop(self):
        self.process.terminate()

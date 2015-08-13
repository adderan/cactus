#!/usr/bin/env python

"""Script to generate a series of random
"""

from optparse import OptionParser
from toil.src.toil.job import Job
from sonLib.bioio import system, spawnDaemon, setLoggingFromOptions, logger, getLogLevelString
from cactus.shared.common import runToilStatusAndFailIfNotComplete

def getDatabaseConf(options):
    return "<st_kv_database_conf type=\"kyoto_tycoon\"><kyoto_tycoon host=\"%s\" port=\"%s\" database_dir=\"%s\"/></<st_kv_database_conf>" % \
        (options.host, options.port, options.databaseDir)

def runDbTestScript(options, firstKey=0, keyNumber=0, addRecords=False, setRecords=False):
    def fn(stringId, bool):
        if bool:
            return stringId
        return ""
    addRecords = fn("--addRecords", addRecords)
    setRecords = fn("--setRecords", setRecords)
    command = "dbTestScript --databaseConf '%s' --firstKey %s --keyNumber %s %s %s --minRecordSize %s --maxRecordSize %s --logLevel %s" %\
    (getDatabaseConf(options), firstKey, keyNumber, addRecords, setRecords, options.minRecordSize, options.maxRecordSize, getLogLevelString())
    system(command)

class AddKeysPhase(Job):
    def __init__(self, options):
        Job.__init__(self)
        self.options = options
    
    def run(self, fileStore):
        keyIndex = 0
        for jobIndex in xrange(int(self.options.totalJobs)):
            self.addChild(AddKeys(self.options, keyIndex))
            keyIndex += int(self.options.keysPerJob)
        self.addFollowOn(SetKeysPhase(self.options))
    
class AddKeys(Job):
    def __init__(self, options, firstKey):
        Job.__init__(self)
        self.options = options
        self.firstKey = firstKey
        
    def run(self, fileStore):
        runDbTestScript(self.options, self.firstKey, self.options.keysPerJob, addRecords=True)

class SetKeysPhase(AddKeysPhase):
    def run(self, fileStore):
        keyIndex = 0
        for jobIndex in xrange(int(self.options.totalJobs)):
            self.addChild(SetKeys(self.options, keyIndex))
            keyIndex += int(self.options.keysPerJob)

class SetKeys(AddKeys):
    def run(self, fileStore):
        runDbTestScript(self.options, self.firstKey, self.options.keysPerJob, setRecords=True)

def main():
    ##########################################
    #Construct the arguments.
    ##########################################

    parser = OptionParser()
    
    parser.add_option("--host", dest="host")
    parser.add_option("--port", dest="port")
    parser.add_option("--databaseDir", dest="databaseDir")
    parser.add_option("--databaseOptions", dest="databaseOptions")
    parser.add_option("--keysPerJob", dest="keysPerJob")
    parser.add_option("--totalJobs", dest="totalJobs")
    parser.add_option("--minRecordSize", dest="minRecordSize")
    parser.add_option("--maxRecordSize", dest="maxRecordSize")
    parser.add_option("--test", dest="test", action="store_true",
                      help="Run doctest unit tests")
    
    Job.Runner.addToilOptions(parser)

    options, args = parser.parse_args()
    if options.test:
        _test()
    setLoggingFromOptions(options)

    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))
    
    Job.Runner.startToil(AddKeysPhase(options), options)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.dbTest.dbTestScript import *
    main()

#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples

class TestCase(unittest.TestCase):
    def testCactusNormalisation_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=TestStatus.getTestSetup())
        
    def testCactusNormalisation_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True)
    
def main():
    parseCactusSuiteTestOptions() 
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()

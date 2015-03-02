
import unittest
import os
import IORobot


class IsOddTests(unittest.TestCase):
    
    def setUp(self):
        print "Init : Copying files : "
        
        self.testingFolder = "unitTestFolder"
        self.mummerPath = "/Users/kakitlam/Desktop/experimentBench/MUMmer3.23/"
        self.listOfFiles = ["raw_reads.fasta", "contigs.fasta"]
        os.system("rm -rf "+ self.testingFolder)
        
    # Remember to have test* at the beginning

    def testFinisherSCMRubTest(self):
        print "Integration test on FinisherSC:  "
        self.runningTestSet("/Users/kakitlam/Desktop/experimentBench/MRub/" , 1)

        
    def testFinisherSCPhTest(self):
        print "Integration test on FinisherSC:  "
        self.runningTestSet("/Users/kakitlam/Desktop/experimentBench/Ph1/" , 5)
        
    def testFinisherSCEcoliTest(self):
        print "Integration test on FinisherSC:  "
        self.runningTestSet("/Users/kakitlam/Desktop/experimentBench/EcoliTestRun/" , 7)

   
            
    def runningTestSet(self ,myFolderName, ctexpected):
        print "Integration test on FinisherSC:  " + myFolderName
        self.sourceFolder = myFolderName
        os.system("mkdir " + self.testingFolder)
        
        for eachitem in self.listOfFiles:
            os.system("cp "+ self.sourceFolder + eachitem + " " +self.testingFolder)
        
        os.system("python finisherSC.py -par 4 "+ self.testingFolder + " "+ self.mummerPath)
        lenDic = IORobot.obtainLength(self.testingFolder, "/improved3.fasta")
        print lenDic
        assert(len(lenDic) == ctexpected)
        os.system("rm -rf "+ self.testingFolder)
        
    def tearDown(self):
        print "Teardown : Removing used files "
        
        

def main():
    unittest.main()
    
if __name__ == '__main__':
    main()



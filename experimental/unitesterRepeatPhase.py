
import unittest
import os


import polyPhaser
import abunSplitter
import tandemRepeatSolver


from finisherSCCoreLib import IORobot

class repeatPhaserTests(unittest.TestCase):
    
    def setUp(self):
        print "Init : Copying files :"
        # Set directories, etc
        self.testingFolder = "unitTestFolder/"
        self.mummerPath = "/Users/kakitlam/Desktop/experimentBench/MUMmer3.23/"
        self.listOfFiles = ["raw_reads.fasta", "contigs.fasta"]
        
        
    def testXphaser(self):
        print "Integration test XPhaser"
        # polyPhaser.py as starting point
        myFolderName = "/Users/kakitlam/Desktop/Research/Nov16-2014/keyAlgorithms/sampleTest2/"
        ctexpected = 2
        commandList = []
        
        commandList.append("python xPhaser.py " + self.testingFolder + " " + self.mummerPath )
        
        
        self.listOfFiles.append("improved3.fasta")
        self.listOfFiles.append("relatedReads.fasta")
        self.listOfFiles.append("raw_reads.part-*")
        
        matchingContigFile = "improved4.fasta"
        self.runningTestSet(myFolderName, ctexpected, commandList, matchingContigFile)
        
        self.listOfFiles.pop()
        self.listOfFiles.pop()
        self.listOfFiles.pop()
        
    def testTSolver(self):
        print "Integration test TSolver"
        # tandemRepeatSolver.py as starting point
        myFolderName = "/Users/kakitlam/Desktop/Research/Nov16-2014/dataFolders/"
        ctexpected = 1
        commandList = []
        
        commandList.append("python finisherSCCoreLib/finisherSC.py -par 4 "+ self.testingFolder +" " + self.mummerPath)
        commandList.append("python tSolver.py " + self.testingFolder + " " + self.mummerPath )
        
        matchingContigFile = "tademResolved.fasta"
        
        self.runningTestSet(myFolderName, ctexpected, commandList, matchingContigFile)
        
        
    def testASplitter(self):
        print "Integration test ASplitter"
        myFolderName = "/Users/kakitlam/Desktop/metaFinisherSC/src/testdata/"
        ctexpected = 2
        commandList = []
        
        commandList.append("python finisherSCCoreLib/finisherSC.py -par 4 "+ self.testingFolder +" " + self.mummerPath)
        commandList.append("python aSplitter.py -par 4 " + self.testingFolder + " " + self.mummerPath )
        
        matchingContigFile = "abun.fasta"
        
        self.runningTestSet(myFolderName, ctexpected, commandList, matchingContigFile)
    
    def testASplitterParameterCheck(self):
        paraTestList = ["-rp improved2.fasta ", "-ar True ", "-rs 0 "]
        for eachpara in paraTestList:
            self.runningParaterTestSet(eachpara,2)
        
        self.runningParaterTestSet("-rd True ",2)
        
    def runningParaterTestSet(self , parameter,ctexpected):
        myFolderName = "/Users/kakitlam/Desktop/metaFinisherSC/src/testdata/"
        commandList = []
        
        commandList.append("python finisherSCCoreLib/finisherSC.py -par 4 "+ self.testingFolder +" " + self.mummerPath)
        commandList.append("python aSplitter.py -par 4 " + parameter + self.testingFolder + " " + self.mummerPath )
        
        matchingContigFile = "abun.fasta"
        
        self.runningTestSet(myFolderName, ctexpected, commandList, matchingContigFile)
    
        
    def runningTestSet(self ,myFolderName, ctexpected, commandList, matchingContigFile):
        
        print "Integration test on RepeatPhaserMain:  " + myFolderName
        self.sourceFolder = myFolderName
        
        os.system("rm -rf "+ self.testingFolder)
        os.system("mkdir " + self.testingFolder)
        
        for eachitem in self.listOfFiles:
            os.system("cp "+ self.sourceFolder + eachitem + " " +self.testingFolder)
        
        for eachcommand in commandList:
            os.system(eachcommand)
        
        lenDic = IORobot.obtainLength(self.testingFolder,  matchingContigFile)
        
        assert(len(lenDic) == ctexpected)
        os.system("rm -rf "+ self.testingFolder)
    
    
    
    def tearDown(self):
        print "Teardown : Removing used files "
        

def main():
    unittest.main()

    
if __name__ == '__main__':
    main()


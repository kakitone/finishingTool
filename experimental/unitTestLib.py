import unittest
import newPhasing


class IsOddTests(unittest.TestCase):
    def setUp(self):
        print "Init"
        self.folderName = "sampleTest2/"
        self.mummerLink ="MUMmer3.23/" 
    def testgetAllAssociatedReads(self):
        print "1: "
        #newPhasing.getAllAssociatedReads(self.folderName, self.mummerLink)
    
    def testformReadContigStringGraph(self):
        print "2: "
        #newPhasing.formReadContigStringGraph(self.folderName, self.mummerLink)
        
    def testidentifyRepeat(self):
        print "3: "
        #newPhasing.identifyRepeat(self.folderName, self.mummerLink)
    
    def testingdefineRepeatAndFlanking(self):
        print "4: "
        #newPhasing.defineRepeatAndFlanking(self.folderName, self.mummerLink)
    def testingperformPhasing(self):
        print "5: "
        newPhasing.performPhasing(self.folderName, self.mummerLink)
        
    def testCombineTest(self):
        print "0 : "
        #newPhasing.mainFlow(self.folderName, self.mummerLink)

def main():
    unittest.main()
    
if __name__ == '__main__':
    main()

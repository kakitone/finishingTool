import numpy as np 
class parameterRobot(object):
    def __init__(self, defaultFolder = ""):
        self.defaultFolder = defaultFolder
        self.setReadStat()
        self.setGenomeStat()
        self.setThresholdPara()
        self.genome = False
        self.longOnly = False
        
    def setReadStat(self, Nshort= 5000, Nlong= 1250, Lshort=100, Llong=240, p=0.1 , longOnly = False):
        self.Nshort = Nshort
        self.Nlong = Nlong
        self.Lshort = Lshort
        self.Llong = Llong
        self.p = p
        self.longOnly = longOnly
        
        
    def setRealGenome(self ,start1 = 3423513    ,start2 = 3689852, lengthOfRep = 5182):
        self.start1, self.start2, self.lengthOfRep, self.genome = start1, start2, lengthOfRep, True
        
        
    def setGenomeStat(self,G = 10000, lrep=500, lsnp=300, lint=50 ):
        self.G, self.lrep, self.lsnp, self.lint = G, lrep, lsnp, lint
        self.lengthOfRep = lsnp
        
    def setThresholdPara(self,   liid = 30, thresForRandom= 0.5,thresForins =0.4, thresFordel=0.4, insMin=4, delMin=4,thresholdForSupport= 0.15, subthreshold= 9, editsub= -10, editins= -1, editdel= -1, editmatch = 1, lookRange =15):
        
        self.liid, self.thresForRandom,self.thresForins , self.thresFordel, self.insMin, self.delMin, self.editsub, self.editins, self.editdel, self.editmatch = liid, thresForRandom,thresForins , thresFordel, insMin, delMin, editsub, editins, editdel, editmatch
        self.thresholdForSupport = thresholdForSupport
        self.subthreshold = subthreshold
        self.lookRange = lookRange
    def tunePara(self):
        covRatio = self.Nshort*self.Lshort/float(30*self.G)
        self.insMin, self.delMin, self.subthreshold = self.insMin*covRatio, self.delMin*covRatio, self.subthreshold*covRatio
    
    
class voteTable(object):
    def __init__(self, index, eachlongread):
        self.index = index    
        self.longread = eachlongread
        self.SNPlocList = []
        
        self.segmentList = []
        
        self.leftSegList = []
        self.rightSegList = []
        
        
    def initVoteTable(self):
        L = len(self.longread)
        
        self.delCount = [0 for i in range(L)]
        self.confirmCount =[0 for i in range(L)]
        # 1, 2,3 ,4 - > 0, 1, 2, 3 
        self.subCount = [[0 , 0 ,0 ,0 ] for i in range(L)]
        
        self .insCount = [ [] for i in range(L) ]
        self.confirmNoIns= [0 for i in range(L)]
        
    def logSNPloc(self, newSNPList):
        
        self.SNPlocList = self.SNPlocList + newSNPList 
        
    
    def filterSNP(self):
        
        self.SNPlocList = sorted(self.SNPlocList)
        index = 0
        while ( index < len(self.SNPlocList) -1 ):
            if self.SNPlocList[index] == self.SNPlocList[index+1]:
                self.SNPlocList.pop(index)
            else:
                index += 1 
                
    
    def matchedIndex(self, seg1, seg2):
        matchedCount = 0 
        minimumMatchingLength = 5
        
        for test1 in seg1:
            for test2 in seg2:
                for i in range(len(test1) - minimumMatchingLength):
                    if test1[i:len(test1)] == test2[0:len(test1) - i]:
                        matchedCount += 1 

        if matchedCount  ==2 :
            return True
        else:
            return False
    
    def filterSNPNeighbor(self):
        toRemoveList = []
        
        for i in range(len(self.segmentList) -1 ):
            if self.matchedIndex(self.segmentList[i], self.segmentList[i+1]):
                toRemoveList.append(i+1)
        
        for j in toRemoveList[-1:-len(toRemoveList)-1:-1]:
            self.segmentList.pop(j)
            self.SNPlocList.pop(j)
            
    def filterSNPOutOfTrustRange(self,trustRange):
        toRemoveList = []
        
        for i in range(len(self.segmentList)):
            LOfSegment = len(self.segmentList[i][0])
            tmp1 = self.segmentList[i][0][max(0, LOfSegment/2 - trustRange): min(LOfSegment,LOfSegment/2 + trustRange )]
            tmp2 = self.segmentList[i][1][max(0, LOfSegment/2 - trustRange): min(LOfSegment,LOfSegment/2 + trustRange )]
            
            print max(0, LOfSegment/2 - trustRange), min(LOfSegment,LOfSegment/2 + trustRange )
            print LOfSegment
            
            print tmp1 == tmp2
            
            if tmp1 == tmp2:
                toRemoveList.append(i)
                print "toRemoveList", toRemoveList
        
                
        for j in toRemoveList[-1:-len(toRemoveList)-1:-1]:
            
            self.segmentList.pop(j)
            self.SNPlocList.pop(j)
            print "self.SNPlocList", self.SNPlocList
    
    def trimendSNP(self,trustRange):
        toRemoveList= []
        
        for i in range(len(self.segmentList)):
            if len(self.longread) - self.SNPlocList[i]  < trustRange : 
                toRemoveList.append(i)
            elif self.SNPlocList[i] < trustRange:
                toRemoveList.append(i)
                
        for j in toRemoveList[-1:-len(toRemoveList)-1:-1]:
            
            self.segmentList.pop(j)
            self.SNPlocList.pop(j)
            print "self.SNPlocList", self.SNPlocList
            
        
                
        
    def filterSegmentList(self, parameterRobot):
        print "Filter segment List"

        if parameterRobot.Llong >1000:
            trustRange = parameterRobot.liid
        else:
            trustRange = parameterRobot.lookRange/3
            
        
            
        #self.filterSNPOutOfTrustRange(trustRange)
        self.filterSNPNeighbor()
        self.trimendSNP(trustRange)
        
        if len(self.segmentList) >=1 :
            self.leftSegList = self.segmentList[0]
            self.rightSegList = self.segmentList[-1]
        else:
            self.leftSegList = []
            self.rightSegList = []
            
            
    def dist(self):
        if len(self.SNPlocList) == 0:
            return len(self.longread)
        else:
            return max(  self.SNPlocList[0], len(self.longread) - self.SNPlocList[-1]   )
            

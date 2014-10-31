import os
from itertools import groupby
from operator import itemgetter
import sys
import time
import argparse


###################################################### Helper Functions 
def extractMumData(folderName, fileName):
    # "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    f = open(folderName + fileName, 'r')
    dataList = []
    
    for i in range(6):
        tmp = f.readline()

    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr = info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        
        matchLenArr = info[2].split()
    
        matchLen1 = int(matchLenArr[0])
        matchLen2 = int(matchLenArr[1])    
        percentMatch = float(info[3])
        
        
        helperStart, helperEnd = int(firstArr[0]), int(firstArr[1])
        readStart, readEnd = int(filterArr[0]) , int(filterArr[1])
        
        helperName = rdGpArr[0].rstrip().lstrip()
        readName = rdGpArr[1].rstrip().lstrip()
        
        dataList.append([helperStart, helperEnd , readStart, readEnd, matchLen1, matchLen2, percentMatch, helperName, readName ])
    
        tmp = f.readline().rstrip()
                
    f.close()
    
    return dataList


def obtainLength(folderName, fileName):
    f = open(folderName + fileName, 'r')
    tmp = f.readline().rstrip()
    lenDic = {}
    tmplen = 0 
    tmpName = ""
    
    while len(tmp) > 0:
        
        if tmp[0] == '>':
            if tmplen != 0:
                lenDic[tmpName] = tmplen
                tmplen = 0
            tmpName = tmp[1:]
        else:
            tmplen += len(tmp)
        tmp = f.readline().rstrip()
        
    lenDic[tmpName] = tmplen
    

    f.close()
    
    return lenDic
    


def findContigLength(folderName, option):
    
    if option == "contigs":
        print "\n\nfindContigLength(folderName)\n"
        contigLength = {}
        f = open(folderName + "smaller_contigs_Double.fasta", 'r')
        
        
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) > 0:
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()
        
        
    else:
        print "\n\nfindContigLength(folderName)\n"
        contigLength = {}
        f = open(folderName + "improved_Double.fasta", 'r')
        
        
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) > 0:
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()
        
    
    
        f = open(folderName + "relatedReads_Double.fasta", 'r')
        
    
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) > 0:
    
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()

    return contigLength




def putListToFileO(folderName, sourceFileName, targetFileName, myList):
    f = open(folderName + targetFileName + ".txt", 'w')
    for eachitem in myList:
        f.write(eachitem + '\n')
        
    f.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' " + folderName + targetFileName + ".txt " + folderName + sourceFileName + " > " + folderName + targetFileName + ".fasta"
    os.system(command)


def transformCoor(dataList):
    # "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    newList = []
    
    for eachitem in dataList:
        if eachitem[2] < eachitem[3]:
            newList.append(eachitem)
        else:
            tmpitem = eachitem
            tmp = tmpitem[2]
            tmpitem[2] = tmpitem[3]
            tmpitem[3] = tmp
            newList.append(tmpitem)
    
    return newList


def useMummerAlign(mummerLink, folderName, outputName, referenceName, queryName, specialForRaw = False, specialName = ""):
    
    if not specialForRaw:
        if globalFast:
            command = mummerLink + "nucmer -b 50  --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " + folderName + queryName
        else:
            command = mummerLink + "nucmer --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " + folderName + queryName
    else:
        if globalFast:
            command = mummerLink + "nucmer -b 50  --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " +  queryName
        else:
            command = mummerLink + "nucmer --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " +  queryName
            
    os.system(command)
    
    if not specialForRaw:
        command = mummerLink + "show-coords -r " + folderName + outputName + ".delta > " + folderName + outputName + "Out"
    else:
        command = mummerLink + "show-coords -r " + folderName + outputName + ".delta > " + folderName + specialName 
        
    os.system(command)


def writeToFile_Double1(folderName, fileName1, fileName2, option="contig"):

    f2 = open(folderName + fileName2, 'w')
    fOriginal = open(folderName + fileName1, 'r')
    
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) > 0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead + tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)
    
    print "len(readSet)", len(readSet)
    
    fOriginal.close()
    
    if option == "contig":
        header = ">Contig"
    else:
        header = ">Read"
    for eachcontig, dum in zip(readSet, range(len(readSet))):
        f2.write(header + str(dum) + "_p\n")
        f2.write(eachcontig + '\n') 
        f2.write(header + str(dum) + "_d\n")
        f2.write(reverseComplement(eachcontig) + '\n')
        
    f2.close()
    
def reverseComplement(myStr):
    myNewStr = myStr[::-1]
    myNewStr2 = ""
    for i in range(len(myNewStr)):
        if myNewStr[i] == 'A' or myNewStr[i] == 'a':
            myNewStr2 += 'T'
            
        elif  myNewStr[i] == 'T' or myNewStr[i] == 't':
            myNewStr2 += 'A'
            
        elif  myNewStr[i] == 'C' or myNewStr[i] == 'c':
            myNewStr2 += 'G'
            
        elif  myNewStr[i] == 'G' or myNewStr[i] == 'g':
            myNewStr2 += 'C'
        elif myNewStr[i] == 'N' or myNewStr[i] == 'n':
            myNewStr2 += 'N'
        else:
            print myNewStr[i]
            
    return myNewStr2

def readConnectList(folderName, fileName):
    f = open(folderName + fileName , 'r')
    tmp = f.readline().rstrip()
    
    dataList = []
    while len(tmp) > 0:
        index, connector, overlap = tmp.split(',')
        dataList.append([int(index), int(connector), int(overlap)])
        tmp = f.readline().rstrip()
        
    connectorList = [ [-1, -1] for i in range(len(dataList))]    
    for eachitem in dataList:
        connectorList[eachitem[0]] = [eachitem[1], eachitem[2]]
        
    f.close()

    return connectorList


def nameInEdgeList(name, myList):
    
    haveInserted = False
    for eachitem in myList:
        if eachitem[0] == name:
            haveInserted = True
    return haveInserted
                

def removeItem(myList, myname):
    newList = []
    for eachitem in myList:
        if eachitem[0] != myname:
            newList.append(eachitem)
    # print len(newList)- len(myList)
    return newList




def quastEvaluate(folderName, quastLink, originalName, improvedNameList, referenceName):
    
    # ./quast.py ~/git/myFinisher/finishingTool/S_cerivisea/contigs.fasta  ~/git/myFinisher/finishingTool/S_cerivisea/improved.fasta -R ~/git/myFinisher/finishingTool/S_cerivisea/reference.fasta
    bindir = os.path.abspath(os.path.dirname(sys.argv[0]))   
    header = bindir + "/" + quastLink + "quast.py" + " "
    
    originalContigPath = folderName + originalName + " "
    improvedContigPath = "" 
    for eachname in improvedNameList:
        improvedContigPath = improvedContigPath + folderName + eachname + " "
        
    
    referencePath = "-R " + folderName + referenceName + " "
    
    command = header + originalContigPath + improvedContigPath + referencePath 
    
    os.system(command)
    
    bindir = os.path.abspath(os.path.dirname(sys.argv[0]))    
    command = "cp " + bindir + "/quast_results/latest/report.txt " + folderName + "assemblyAssessment.txt"
    os.system(command)

def compareGraphUnitTest(G, G2):
    assert(len(G.graphNodesList) == len(G2.graphNodesList))
    for index in range(len(G.graphNodesList)):
        assert(G.graphNodesList[index].nodeIndex == G2.graphNodesList[index].nodeIndex)
        assert(G.graphNodesList[index].nodeIndexList == G2.graphNodesList[index].nodeIndexList)
        assert(G.graphNodesList[index].overlapList == G2.graphNodesList[index].overlapList)
        assert(G.graphNodesList[index].listOfPrevNodes == G2.graphNodesList[index].listOfPrevNodes)
        assert(G.graphNodesList[index].listOfNextNodes == G2.graphNodesList[index].listOfNextNodes)
        assert(G.graphNodesList[index].visited == G2.graphNodesList[index].visited)


def loadContigsFromFile(folderName, fileName):
    f = open(folderName + fileName, 'r')
    tmp = f.readline().rstrip()
    dataDic = {}
    tmpSeq = ""
    tmpName = ""
    
    while len(tmp) > 0:
        
        if tmp[0] == '>':
            if len(tmpSeq) != 0:
                dataDic[tmpName] = tmpSeq
                tmpSeq = ""
            tmpName = tmp[1:]
        else:
            tmpSeq += tmp
        tmp = f.readline().rstrip()
        
    dataDic[tmpName] = tmpSeq
    

    f.close()
    return dataDic


class seqGraphNode(object):
    def __init__(self, nodeIndex):
        self.nodeIndex = nodeIndex
        self.nodeIndexList = [nodeIndex]
        self.overlapList = []
        self.listOfPrevNodes = []
        self.listOfNextNodes = []
        self.visited = 0


class seqGraph(object):
    def __init__(self, numberOfNodes):
        self.graphNodesList = [seqGraphNode(i) for i in range(numberOfNodes)]
    
    def loadFromFile(self, folderName, fileName):
        
        f = open(folderName + fileName, 'r')
        numberOfNodes = 0 
        tmp = f.readline().rstrip()
        while (len(tmp) > 0):
            tmp = f.readline().rstrip()
            numberOfNodes += 1
        f.close()
        
        self.graphNodesList = [seqGraphNode(i) for i in range(numberOfNodes)]
        
        f = open(folderName + fileName, 'r')
        tmp = f.readline().rstrip()
        runningIndex = 0
        while (len(tmp) > 0):
            dataList = tmp.split(';')
            

            for i in range(6):
                if len(dataList[i]) > 0:
                    myList = dataList[i].split(',')
                else:
                    myList = []

                    
                if i == 0:
                    self.graphNodesList[runningIndex].nodeIndex = int(dataList[i])
                elif i == 1:
                    self.graphNodesList[runningIndex].nodeIndexList = []
                    for eachitem in myList:
                        self.graphNodesList[runningIndex].nodeIndexList.append(int(eachitem))
                        
                elif i == 2:
                    self.graphNodesList[runningIndex].overlapList = []
                    for eachitem in myList:
                        self.graphNodesList[runningIndex].overlapList.append(int(eachitem))
                elif i == 3 or i == 4:
                    for eachitem in myList:
                        mydata = eachitem.split('-')
                        if i == 3 :
                            self.graphNodesList[runningIndex].listOfPrevNodes.append([int(mydata[0]), int(mydata[1])])
                        elif i == 4:
                            self.graphNodesList[runningIndex].listOfNextNodes.append([int(mydata[0]), int(mydata[1])])
                elif i == 5:
                    self.graphNodesList[runningIndex].visited = int(dataList[i])
            
                
            tmp = f.readline().rstrip()
            runningIndex = runningIndex + 1 
        f.close()
    
    def insertEdge(self, i, j, wt):
        if j != -1 and i != -1:
            haveInserted = nameInEdgeList(j, self.graphNodesList[i].listOfNextNodes)

            if not haveInserted:      
                self.graphNodesList[i].listOfNextNodes.append([j, wt])
                
                
            haveInserted = nameInEdgeList(i, self.graphNodesList[j].listOfPrevNodes) 
            
            if not haveInserted:      
                self.graphNodesList[j].listOfPrevNodes.append([i, wt])
                


    def findStartEndList(self):
        # ## A strict method 
        self.myStartList = []
        self.myEndList = []
        self.repeatList = []
        
        for eachnode in self.graphNodesList:
            if len(eachnode.nodeIndexList) > 0:
                if len(eachnode.listOfPrevNodes) == 0 :
                    self.myStartList.append(eachnode.nodeIndex)
                    # self.myStartList= self.myStartList + eachnode.nodeIndexList
   
                if len(eachnode.listOfNextNodes) == 0 :
                    # self.myEndList = self.myEndList + eachnode.nodeIndexList
                    self.myEndList.append(eachnode.nodeIndex)
                    
                if len(eachnode.listOfNextNodes) >= 2 or len(eachnode.listOfPrevNodes) >= 2 :
                    self.repeatList.append(eachnode.nodeIndex)
        
                    
        self.myStartList.sort()
        self.myEndList.sort()
        # print self.myStartList    
        # print self.myEndList
        
    def reportEdge(self):
        for eachnode in self.graphNodesList:
            myName = eachnode.nodeIndex
            if len(eachnode.nodeIndexList) > 0:
                for eachnextnode in eachnode.listOfNextNodes:
                    nextName = eachnextnode[0]
                    wt = eachnextnode[1]
                    if wt > 1000:
                        print str(myName) + "->" + str(nextName) + "{weight:" + str(int(wt / 1000)) + "}"
                    else:
                        print str(myName) + "->" + str(nextName) + "{weight:" + str(0.5) + "}"
                        
        for eachnode in self.graphNodesList:
            myName = eachnode.nodeIndex
            if len(eachnode.nodeIndexList) > 0:
                if len(eachnode.listOfNextNodes) == 0 and len(eachnode.listOfPrevNodes) == 0:
                    print str(myName) + "{color:red}"


    def cleanEdge(self):
        toCleanList = []
        for eachnode in self.graphNodesList:
            if len(eachnode.nodeIndexList) > 0:
                myName = eachnode.nodeIndex
                for eachnextnode in eachnode.listOfNextNodes:
                    nextName = eachnextnode[0]
                    wt = eachnextnode[1]
                    
                    if myName % 2 == 0 :
                        myCom = myName + 1 
                    else:
                        myCom = myName - 1 
                    
                    if nextName % 2 == 0 :
                        nextCom = nextName + 1 
                    else:
                        nextCom = nextName - 1 
                    
                    if not  nameInEdgeList(myCom, self.graphNodesList[nextCom].listOfNextNodes) :
                        toCleanList.append([myName, nextName])
                
        print "len(toCleanList)", len(toCleanList) 
        
        for eachitem in toCleanList:
            myName, nextName = eachitem
            self.removeEdge(myName, nextName)
            
    def removeEdge(self, myName, nextName):
        self.graphNodesList[myName].listOfNextNodes = removeItem(self.graphNodesList[myName].listOfNextNodes, nextName)
        self.graphNodesList[nextName].listOfPrevNodes = removeItem(self.graphNodesList[nextName].listOfPrevNodes, myName)
                

    def condense(self):
        for eachnode in self.graphNodesList:
            myname = eachnode.nodeIndex
            
            for eachnextnode in eachnode.listOfNextNodes:
                nextname = eachnextnode[0]
                wt = eachnextnode[1]
                
                if len(eachnode.listOfNextNodes) == 1 and len(self.graphNodesList[nextname].listOfPrevNodes) == 1:

                    originalPrev = self.graphNodesList[myname].listOfPrevNodes
                    
                    self.graphNodesList[nextname].nodeIndexList = self.graphNodesList[myname].nodeIndexList + self.graphNodesList[nextname].nodeIndexList
                    self.graphNodesList[nextname].overlapList = self.graphNodesList[myname].overlapList + [wt] + self.graphNodesList[nextname].overlapList
                    
                    haveInserted = nameInEdgeList(myname, self.graphNodesList[nextname].listOfPrevNodes)
                    if haveInserted :
                        self.graphNodesList[nextname].listOfPrevNodes = removeItem(self.graphNodesList[nextname].listOfPrevNodes, myname)
                    
                    for eachp in originalPrev:  
                        if not eachp in    self.graphNodesList[nextname].listOfPrevNodes:
                            self.graphNodesList[nextname].listOfPrevNodes.append(eachp)
                         
                    self.graphNodesList[myname].nodeIndexList = [] 
                    self.graphNodesList[myname].overlapList = []
                    self.graphNodesList[myname].listOfNextNodes = []
                    self.graphNodesList[myname].listOfPrevNodes = []
                
                    for eachoriginalprev in originalPrev:
                        prevname = eachoriginalprev[0]
                        # if [myname,eachoriginalprev[1]] in self.graphNodesList[prevname].listOfNextNodes:
                        haveInserted = nameInEdgeList(myname, self.graphNodesList[prevname].listOfNextNodes)
                        if haveInserted:
                            self.graphNodesList[prevname].listOfNextNodes = removeItem(self.graphNodesList[prevname].listOfNextNodes, myname)
                            # self.graphNodesList[prevname].listOfNextNodes.remove([myname,eachoriginalprev[1]])
                        
                        self.graphNodesList[prevname].listOfNextNodes.append([nextname, eachoriginalprev[1]])
                    
                        
    def reportDummyUsefulNode(self):
        countUseful = 0
        countUseless = 0 
        for eachnode in self.graphNodesList:
            # print eachnode.listOfNextNodes
            if len(eachnode.nodeIndexList) > 0:
                countUseful += 1
            else:
                countUseless += 1
                
        print "remained, condensed" , countUseful, countUseless
      
    
    def checkSelfLoops(self):
        ck = True
        for eachnode in self.graphNodesList:
            myIndexList = eachnode.nodeIndexList
           # print "eachnode.listOfNextNodes", len(eachnode.listOfNextNodes)
           # print "eachnode.listOfPrevNodes", len(eachnode.listOfPrevNodes)
            
            for eachnext in eachnode.listOfNextNodes:
                if eachnext[0] in myIndexList:
                    # print  "selfLoop" , eachnode.nodeIndex, eachnode.nodeIndexList
                    ck = False
        
        print "No self-loop existence ?", ck
        
            
    
    
    def checkCompleteness(self):
        originalNumberOfNodes = len(self.graphNodesList)
        totalList = []
        for eachnode in self.graphNodesList:
            myIndexList = eachnode.nodeIndexList
            if len(myIndexList) > 0:
                totalList = totalList + myIndexList
        
        print "nodeIndexList conserved ? ", len(totalList) == len(set(totalList)) and len(set(totalList)) == originalNumberOfNodes
            
        
    def saveToFile(self, folderName, fileName):
        # ## Format : nodeIndex,nodeIndexList, overlapList, listOfPrevNodes, listOfNextNodes, visited
        # 3 ; [5, 3] ; [4757]; [[450, 3887], [123,5678]]  ; []; 0
        # 3 ; 5, 3 ; 4757; 450-3887,123-5678  ; ; 0
        
        f = open(folderName + fileName , 'w')
        for eachnode in self.graphNodesList:
            mystr = ""
            mystr = mystr + str(eachnode.nodeIndex) + ';'
            
            for eachitem, index  in  zip(eachnode.nodeIndexList, range(len(eachnode.nodeIndexList))):
                if index != len(eachnode.nodeIndexList) - 1:
                    mystr = mystr + str(eachitem) + ','
                else: 
                    mystr = mystr + str(eachitem)
                    
            mystr = mystr + ';'
            
            for eachitem, index in  zip(eachnode.overlapList, range((len(eachnode.overlapList)))):
                if index != len(eachnode.overlapList) - 1:
                    mystr = mystr + str(eachitem) + ','
                else:
                    mystr = mystr + str(eachitem)
                
            mystr = mystr + ';'
            
            for eachitem, index in  zip(eachnode.listOfPrevNodes, range(len(eachnode.listOfPrevNodes))):
                if index != len(eachnode.listOfPrevNodes) - 1:
                    mystr = mystr + str(eachitem[0]) + '-' + str(eachitem[1]) + ','
                else:
                    
                    mystr = mystr + str(eachitem[0]) + '-' + str(eachitem[1]) 
            mystr = mystr + ';'
            
            for eachitem, index in  zip(eachnode.listOfNextNodes, range(len(eachnode.listOfNextNodes))):
                if index != len(eachnode.listOfNextNodes) - 1:
                    mystr = mystr + str(eachitem[0]) + '-' + str(eachitem[1]) + ','
                else:
                    
                    mystr = mystr + str(eachitem[0]) + '-' + str(eachitem[1]) 
            mystr = mystr + ';'
                

            mystr = mystr + str(eachnode.visited) + '\n'
            f.write(mystr)
            
        f.close()


    def MBResolve(self):
        # ## Algorithm for resolving the repeats in string graph by MB 
        # 1 ) Log down the 1 successor / 1 predecessor nodes[cf: edge] and their associated weights
        # 2 ) Sort according to edge weights
        # 3 ) Make commitment and  remove edges from both sides
        # 4 ) Condense graph 
        # 5 ) Go back to 1 unless no such nodes exist any more 
        
        
        oneSucList, onePreList = [-1] , [-1]
        
        while len(oneSucList) > 0 or len(onePreList) > 0:
            oneSucList, onePreList = [] , []
            
            for eachnode in self.graphNodesList:
                if len(eachnode.nodeIndexList) > 0 :
                    if len(eachnode.listOfNextNodes) == 1 :
                        nextNodeIndex = eachnode.listOfNextNodes[0][0]
                        if len(self.graphNodesList[nextNodeIndex].listOfPrevNodes) <= 2:
                            oneSucList.append([eachnode.nodeIndex, eachnode.listOfNextNodes[0][0], eachnode.listOfNextNodes[0][1] , 1 ])
                    if len(eachnode.listOfPrevNodes) == 1: 
                        prevNodeIndex = eachnode.listOfPrevNodes[0][0]
                        if len(self.graphNodesList[prevNodeIndex].listOfNextNodes) <= 2:
                            onePreList.append([eachnode.listOfPrevNodes[0][0], eachnode.nodeIndex, eachnode.listOfPrevNodes[0][1], 0 ])
                            
            
             
            print "oneSucList", oneSucList 
            print "onePreList", onePreList
            
            
            # ## Do them together when sort ? 
            combinedList = oneSucList + onePreList 
            combinedList.sort(key=itemgetter(2), reverse=True)
           
            
            for eachitem in combinedList:
                i, j, wt, myType = eachitem[0], eachitem[1] , eachitem[2], eachitem[3]
                
                if myType == 1:
                    if nameInEdgeList(i , self.graphNodesList[j].listOfPrevNodes):
                        removeList = []
                        for eachprev in self.graphNodesList[j].listOfPrevNodes:
                            if eachprev[0] != i :
                                removeList.append(eachprev[0])
                        
                        for eachToRemove in removeList:
                            self.removeEdge(eachToRemove, j)
    
                elif myType == 0:
                    
                    if nameInEdgeList(j , self.graphNodesList[i].listOfNextNodes):
                        removeList = []
                        for eachnext in self.graphNodesList[i].listOfNextNodes:
                            if eachnext[0] != j :
                                removeList.append(eachnext[0])
                        
                        for eachToRemove in removeList:
                            self.removeEdge(i, eachToRemove)
    
                
            self.condense()
        

def loadOpenList(folderName):
    
    f = open(folderName + 'openZone.txt', 'r')
    tmp = f.readline().rstrip()
    N = int(tmp)
    usableJunction = [[ 0, 0 ] for i in range(N) ]
    
    tmp = f.readline().rstrip()
    while len(tmp) > 0:
        myInfo = tmp.split(',')
        contigIndex = int(myInfo[0][5:])
        if myInfo[1] == 'noprev':
            usableJunction[contigIndex][0] = 1
        if myInfo[1] == 'nonext' :
            usableJunction[contigIndex][1] = 1
             
        tmp = f.readline().rstrip()
    
    f.close() 
    
    print "usableJunction", usableJunction
    return usableJunction




def writeToFile(f2, runningIndex, seq):
    f2.write(">Seg_" + str(runningIndex))
    f2.write('\n')
    f2.write(seq)
    f2.write('\n')
   



def formRelatedReadsFile(folderName, mummerLink):    
    # Find associated read and extract into a file associatedReads.fasta
    # Input: contigs.fasta, cleaned_Reads.fasta 
    # Output: relatedReads.fasta

    # ## Extract heads of the contigs
    print ">formRelatedReadsFile"
    

    f = open(folderName + "improved.fasta", 'r')
    f2 = open(folderName + "improvedTrunc.fasta", 'w')
    temp = f.readline()
    tempContig = ""
    thres = 400
    runningIndex = 0
    endThres = 10 
    
    while len(temp) > 0:
        if temp[-1] == '\n':
            temp = temp[0:-1]
        
        
        if temp[0] == '>':

            if len(tempContig) > 0:
                writeToFile(f2, runningIndex, tempContig[0:thres])
                runningIndex = runningIndex + 1
                
                writeToFile(f2, runningIndex, tempContig[-thres:])
                runningIndex = runningIndex + 1 
                
                                
                writeToFile(f2, runningIndex, reverseComplement(tempContig[0:thres]))
                runningIndex = runningIndex + 1
                
                writeToFile(f2, runningIndex, reverseComplement(tempContig[-thres:]))
                runningIndex = runningIndex + 1
                
                tempContig = ""
        else:
            tempContig = tempContig + temp
        
        temp = f.readline()

    writeToFile(f2, runningIndex, tempContig[0:thres])
    runningIndex = runningIndex + 1
    
    writeToFile(f2, runningIndex, tempContig[-thres:])
    runningIndex = runningIndex + 1
                  
    writeToFile(f2, runningIndex, reverseComplement(tempContig[0:thres]))
    runningIndex = runningIndex + 1
    
    writeToFile(f2, runningIndex, reverseComplement(tempContig[-thres:]))
    runningIndex = runningIndex + 1
    
    
    f2.close()
    f.close()
    
    # ## Write double stranded reads
    writeToFile_Double1(folderName, "improved.fasta", "improved_Double.fasta", "contig")
    # writeToFile_Double1(folderName, "raw_reads.fasta", "raw_reads_Double.fasta","read")
    
    # ## Apply MUMMER on them using cleanedReads against them
    assoiatedReadIndex = []
    nameList = []
    
    numberOfFiles = 20
    if True:
        bindir = os.path.abspath(os.path.dirname(sys.argv[0]))
        command = bindir + "/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + folderName + "raw_reads.fasta"
        os.system(command)
    
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
        
        if True:
            useMummerAlign(mummerLink, folderName, "out", "improvedTrunc.fasta", "raw_reads.part-" + indexOfMum + ".fasta", True, "fromMum" + indexOfMum )
            
            '''
            command = mummerLink + "nucmer --maxmatch --nosimplify -p " + folderName + "out " + folderName + "improvedTrunc.fasta raw_reads.part-" + indexOfMum + ".fasta"
            os.system(command)
    
            command = mummerLink + "show-coords -r " + folderName + "out.delta > " + folderName + "fromMum" + indexOfMum
            os.system(command)
            '''
            
        f = open(folderName + "fromMum" + indexOfMum, 'r')
    
        for i in range(6):
            tmp = f.readline()
        
        while len(tmp) > 0:
            infoArr = tmp.split('|')
            myArr = infoArr[-1].split('\t')
            rdGpArr = infoArr[-1].split('\t')
            contigName = rdGpArr[0].rstrip().lstrip()
            readName = rdGpArr[1].rstrip().lstrip()
            
            endSegArr = infoArr[0].split(" ")
            pos = []
            for eachitem in endSegArr:
                if len(eachitem) > 0:
                    pos.append(int(eachitem))
                    
            startPos = pos[0]
            endPos = pos[1]
            if startPos < endThres and endPos > thres - endThres:
                assoiatedReadIndex.append(myArr[1])
                nameList.append([int(contigName.split('_')[1]), readName])
            tmp = f.readline()
        
        f.close()
    
    
    nameList.sort()

    assoiatedReadIndex.sort()
    
    # print "assoiatedReadIndex", assoiatedReadIndex
    
    ckIndex = 0
    f = open(folderName + "associatedNames.txt", 'w')
    oneItem = 0
    keyFound = []
    for key, items in groupby(assoiatedReadIndex):
        
        countItem = 0
        for eachitem in items:
            countItem += 1
            
        if countItem == 1:
            
            oneItem += 1
        else:
            key = key.rstrip()
            if not key in keyFound:
                f.write(key + '\n')
                keyFound.append(key)

        ckIndex += 1
    
    print "ckIndex,oneItem: ", ckIndex, oneItem
    f.close()

    fFilter = open(folderName + "associatedNames.txt", 'r')
    
    fout = open(folderName + "associatedNames2.txt", 'w') 
    
    maxCount = 12000
    mytmpDum = fFilter.readline() 
    i = 0
    while i < maxCount and len(mytmpDum) > 0:
        fout.write(mytmpDum)  
        mytmpDum = fFilter.readline() 
        i = i + 1
        
    fout.close()   
    fFilter.close()

    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' " + folderName + "associatedNames2.txt " + folderName + "raw_reads.fasta > " + folderName + "relatedReads.fasta"
    os.system(command)
    
    writeToFile_Double1(folderName, "relatedReads.fasta", "relatedReads_Double.fasta", "read")
    
    numberOfFiles = 20
    if True:
        bindir = os.path.abspath(os.path.dirname(sys.argv[0]))
        command = bindir + "/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + folderName + "relatedReads_Double.fasta"
        os.system(command)
    

def extractEdgeSet(folderName, mummerLink, option="nopolish"):
    # Tasks: reconstruct the string  graph
    
    # Input : relatedReads_Double.fasta, conig_Double.fasta
    # Intermediate files: fromMum_overlap , fromMum_overlap
    # Output: connectivity of eachNode: InList, OutList [critical]
    #         connectivity of eachNode: arrow representatiowith size [optional]
    
    
    # ## Perform MUMMER alignment
    print ">Extract Edge set"
    # lengthDic = obtainLength(folderName, "improved_Double.fasta")
    # print lengthDic
    lengthDic = findContigLength(folderName, "improved")
    
    numberOfContig = 0
    f = open(folderName + "improved_Double.fasta", 'r')
    tmp = f.readline()
    tmp = tmp.rstrip()
    while (len(tmp) > 0):
        numberOfContig += 1
        tmp = f.readline()
        tmp = tmp.rstrip()
    
    numberOfContig = numberOfContig / 2
    print "numberOfContig", numberOfContig

    f.close()
    K = 400
    
    thres = 5
    # Nodes are contigs, 
    nodes = [i for i in range(numberOfContig)]
    dataSet = []
    
    # ## Apply MUMMER on them using cleanedReads against them
    fmyFile = open(folderName + "improved_Double.fasta", 'r')
    fSmaller = open(folderName + "smaller_improvedContig.fasta", 'w')

    tmp = fmyFile.readline().rstrip()
    maxSize = 25000
    
    dummySeq = ""
    for i in range(0):
        dummySeq = dummySeq + "A"
        
    
    myName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            fSmaller.write(tmp + '\n')
            myName = tmp[1:]
        else:
            component = tmp[0:min(len(tmp), maxSize)] 
            countComp = len(component)
            fSmaller.write(component)
            
            fSmaller.write(dummySeq)
            countComp = countComp + len(dummySeq)
            
            component = tmp[max(0, len(tmp) - maxSize):len(tmp)]
            fSmaller.write(component)
            countComp = countComp + len(component)
            
            lengthDic[myName] = countComp 
            print "DebugName", myName, countComp
            fSmaller.write('\n')

        tmp = fmyFile.readline().rstrip()

    fSmaller.close()
    fmyFile.close()
    
    numberOfFiles = 20
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)

        if True:
            useMummerAlign(mummerLink, folderName, "outRefine", "smaller_improvedContig.fasta", "relatedReads_Double.part-" + indexOfMum + ".fasta", True,  "fromMumRefine" + indexOfMum)
            
            '''
            command = mummerLink + "nucmer --maxmatch --simplify -p " + folderName + "outRefine " + folderName + "smaller_improvedContig.fasta " + "relatedReads_Double.part-" + indexOfMum + ".fasta"
            os.system(command)
            
            command = mummerLink + "show-coords -r " + folderName + "outRefine.delta > " + folderName + "fromMumRefine" + indexOfMum
            os.system(command)
            '''
            
        f = open(folderName + "fromMumRefine" + indexOfMum, 'r')
        for i in range(6):
            tmp = f.readline()
    
        while len(tmp) > 0:
            info = tmp.split('|')
            filterArr = info[1].split()
            rdGpArr = info[-1].split('\t')
            firstArr = info[0].split()
            matchLenArr = info[2].split()
           
        
            matchLen = int(matchLenArr[1])    
            contigStart, contigEnd = int(firstArr[0]), int(firstArr[1])
            readStart, readEnd = int(filterArr[0]) , int(filterArr[1])
            
            contigName = rdGpArr[0].rstrip().lstrip()
            readName = rdGpArr[1].rstrip().lstrip()
                    
            if readStart < readEnd and matchLen > K and min(contigStart, readStart) < thres and min(lengthDic[contigName] - contigEnd , lengthDic[readName] - readEnd) + 1 < thres:
                conditionForMatch = True
            else:
                conditionForMatch = False
    
            if conditionForMatch :
                if contigStart < thres:
                    dataSet.append((readName, contigName, 'L', matchLen))
                
                if lengthDic[contigName] - contigEnd + 1 < thres :
                    dataSet.append((readName, contigName, 'R', matchLen))
            
            tmp = f.readline()
            
        f.close()  
        
    
    # ## repeat aware
    usableJunction = loadOpenList(folderName)
    dataSet, blockedSet = filterRepeatEnd(dataSet, usableJunction)
    # ## repeat aware end
    
    dataSet.sort()
    matchPair = formMatchPairFromReadInfo(dataSet)
    

    
    # print matchPair
    keyFound = []
    bestMatchPair = []
    rawReadList = []
    
    for key, items in groupby(matchPair, itemgetter(0, 1)):
        maxvalue = -1
        maxLenPair = []
        for eachitem in items:
            if eachitem[2] > maxvalue:
                maxvalue = eachitem[2]
                maxLenPair = [eachitem[3], eachitem[4], eachitem[5]]
        bestMatchPair.append([key[0], key[1], maxvalue, maxLenPair[0], maxLenPair[1], maxLenPair[2]])
    
    bestMatchPair.sort(key=lambda tup: tup[2], reverse=True)
    
    print "bestMatchPair", bestMatchPair


    leftConnect = [[-1, -1, -1] for i in range(numberOfContig)]
    rightConnect = [[-1, -1, -1] for i in range(numberOfContig)]
    
    for eachitem in bestMatchPair:

        prefixContig = eachitem[0]
        suffixContig = eachitem[1]
        
        if leftConnect[suffixContig][0] == -1 and rightConnect[prefixContig][0] == -1:
            leftConnect[suffixContig][0] = prefixContig 
            leftConnect[suffixContig][1] = eachitem[4]
            leftConnect[suffixContig][2] = eachitem[5]
            
            rightConnect[prefixContig][0] = suffixContig
            rightConnect[prefixContig][1] = eachitem[3]
            rightConnect[prefixContig][2] = eachitem[5]
            
            rawReadList.append(eachitem[5])
    
    startList = []
    print "leftConnect", leftConnect
    for i in range(len(leftConnect)):
        if leftConnect[i][0] == -1:
            startList.append(i)
    
    print "startList", startList
    checkLoopList = [False for i in range(len(leftConnect))]
        
    for eachitem in startList:
        tmp = eachitem
        checkLoopList[tmp] = True
        while rightConnect[tmp][0] != -1:
            tmp = rightConnect[tmp][0]
            if tmp != -1:
                checkLoopList[tmp] = True
    
    for dumdumi in range(len(leftConnect)):
        if checkLoopList[dumdumi] == False:
            startList.append(dumdumi)
            tmp = dumdumi
            
            while checkLoopList[tmp] == False :
                if tmp != -1:
                    checkLoopList[tmp] = True
                    
                tmp = rightConnect[tmp][0]


    contigList = []
    print "startList", startList
    for eachitem in startList:
        tmp = eachitem
        myList = [tmp]
        mystart = tmp
        while rightConnect[tmp][0] != -1 and rightConnect[tmp][0] != mystart:
            tmp = rightConnect[tmp][0]
            if tmp != -1:
                myList.append(tmp)
        contigList.append(myList)
    
    print "contigList", contigList

    
    # ## repeat aware logging
    myExtraLinkList = loggingReadsToRepeat(blockedSet + dataSet, contigList)
    # print "myExtraLinkList", myExtraLinkList
    # ## end repeat aware logging
    
    
    i = 0
    fOriginal = open(folderName + "improved.fasta", 'r')
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) > 0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead + tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)  
    fOriginal.close()
    
    # ## Put the needed rawReads into the RAM using Dictionary
    
    fAppendRaw = open(folderName + "appendRaw.txt", 'w')
    for eachraw in rawReadList:
        fAppendRaw.write(eachraw)
        fAppendRaw.write('\n')
    fAppendRaw.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' " + folderName + "appendRaw.txt " + folderName + "relatedReads_Double.fasta > " + folderName + "rawToAppend.fasta"
    os.system(command)
    
    rawRead = {}
    
    fOriginal = open(folderName + "rawToAppend.fasta", 'r')
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    tmpName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            
            if len(tmpRead) > 0:
                rawRead[tmpName] = tmpRead
                tmpRead = ""
                
            tmpName = tmp[1:]
        else:
            tmpRead = tmpRead + tmp
            
        tmp = fOriginal.readline().rstrip()
        
    rawRead[tmpName] = tmpRead
    # ## End
    
    seqToPrint = []
    contigUsed = [False for i in range(numberOfContig / 2)]
    storedStrand = [[-1, 'n'] for i in range(numberOfContig)]
    
    fAllInOne = open(folderName + "allInOne2.fasta", 'w')
    finalList = []
    for eachContig, i in zip(contigList, range(len(contigList))):
        tmpList = []
        for eachitem in eachContig:

            readNum = eachitem / 2
            if contigUsed[readNum] == False:
                seqToPrint.append(eachitem)
                tmpList.append(eachitem)
                contigUsed[readNum] = True
                # ## mark ouput strandinfo
                storedStrand[eachitem] = [len(finalList), 'p']
                
                
        if len(tmpList) > 0:
            finalList.append(tmpList)
    
    
    for kkk in range(len(storedStrand)):
        if storedStrand[kkk][1] == 'n':
            if kkk % 2 == 0:
                storedStrand[kkk][0] = storedStrand[kkk + 1][0]
                storedStrand[kkk][1] = 'd'
            else:
                storedStrand[kkk][0] = storedStrand[kkk - 1][0]
                storedStrand[kkk][1] = 'd'
    
    # ## begin stored output blocked pairs
    blockExtraStored(storedStrand, myExtraLinkList, folderName)
    # ## end output blocked pairs stored
    
    
    fImproved = open(folderName + "improved2.fasta", 'w')
    
    for eachcontig, dummyIndex in zip(finalList, range(len(finalList))):
        fImproved.write(">Segkk" + str(dummyIndex) + '\n')
        tmpStore = -1997
        tmpStore2 = -1998
        tmpStore3 = -1999
        
        for eachseg, hidum in zip(eachcontig, range(len(eachcontig))):
            readNum = eachseg / 2
            orientation = eachseg % 2
            
            newStart = 0 
            
            x , y , l = tmpStore, leftConnect[eachseg][1], tmpStore2
            extraRead = ""
            if hidum == 0:
                newStart = 0
            else:
                
                if l < x + y:
                    newStart = x + y - l
                else:
                    newStart = 0
                    if option == 'polish':
                        print "Missing polish"
                        extraRead = tmpStore3[x:l - y]
                        # extraRead = performPolishing(leftConnect[eachseg][0], eachseg, tmpStore3[x:l-y],  dataSet, folderName)
                    else:
                        extraRead = tmpStore3[x:l - y]
                    
    
            print extraRead[0:10], len(extraRead)
            
            
            fImproved.write(extraRead)
            
            if orientation == 0:
                fImproved.write(readSet[readNum][newStart:])   

            else:
                fImproved.write(reverseComplement(readSet[readNum])[newStart:])

            
            if rightConnect[eachseg][1] != -1:
                tmpStore = rightConnect[eachseg][1]
                tmpStore2 = len(rawRead[rightConnect[eachseg][2]])
                tmpStore3 = rawRead[rightConnect[eachseg][2]]
                
        fImproved.write('\n')
        
    fImproved.close()
            


    for eachContigIndex, dum in zip(seqToPrint, range(len(seqToPrint))):

        readNum = eachContigIndex / 2
        orientation = eachContigIndex % 2
        
        fAllInOne.write(">Seg_" + str(dum) + '\n')
        if orientation == 0:
            fAllInOne.write(readSet[readNum])
            fAllInOne.write('\n')
        else:
            fAllInOne.write(reverseComplement(readSet[readNum]))
            fAllInOne.write('\n')
                    
    fAllInOne.close()
    


    for eachcontig, i in zip(finalList, range(len(finalList))):
        fout = open(folderName + "improvedContig2_" + str(i) + ".fasta", 'w')
        seqToPrint = []

        for eachitem in eachcontig:

            readNum = eachitem / 2
            seqToPrint.append(eachitem)

        print "ImprovedContig ", i 
        for eachhaha in seqToPrint:
            print len(readSet[eachhaha / 2])
            
        for eachContigIndex, dum in zip(seqToPrint, range(len(seqToPrint))):

            readNum = eachContigIndex / 2
            orientation = eachContigIndex % 2
            
            fout.write(">Seg_" + str(dum) + '\n')
            if orientation == 0:
                fout.write(readSet[readNum])
                fout.write('\n')
            else:
                fout.write(reverseComplement(readSet[readNum]))
                fout.write('\n')


        i += 1
        fout.close()


def filterRepeatEnd(dataSet, usableJunction):
    newDataSet = []
    blockedSet = []
    # Format : dataSet.append((readName, contigName, 'R',matchLen)) 
    for eachitem in dataSet:
        readName, contigName, seekDir, matchLen = eachitem
        # print contigName
        myInfo = contigName[6:].split('_')
        index, orient = int(myInfo[0]), myInfo[1]
        
        # print "index, orient", index, orient
        check = False
        if orient == 'p':
            if usableJunction[index][0] == 1 and seekDir == 'L':
                check = True
            elif usableJunction[index][1] == 1 and seekDir == 'R':
                check = True
                
        elif orient == 'd':
            if usableJunction[index][1] == 1 and seekDir == 'L':
                check = True
            elif usableJunction[index][0] == 1 and seekDir == 'R':
                check = True

        if check:
            newDataSet.append(eachitem)
        else:
            blockedSet.append(eachitem)
    


    return newDataSet, blockedSet

def obtainLinkInfo(folderName, mummerLink, inputFile, mummerFile):
    thres = 5
    minLen = 400
    # thres = 10
    # minLen = 200
    
    
    writeToFile_Double1(folderName, inputFile + ".fasta", inputFile + "_Double.fasta", "contig")
    
    fmyFile = open(folderName + inputFile + "_Double.fasta", 'r')
    fSmaller = open(folderName + inputFile + "_contigs_Double.fasta", 'w')

    tmp = fmyFile.readline().rstrip()
    maxSize = 50000

    myName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            fSmaller.write(tmp + '\n')
            myName = tmp[1:]
        else:
            component = tmp[0:min(len(tmp), maxSize)] 
            countComp = len(component)
            fSmaller.write(component)
            
            component = tmp[max(0, len(tmp) - maxSize):len(tmp)]
            fSmaller.write(component)
            countComp = countComp + len(component)
            

            print "DebugName", myName, countComp
            fSmaller.write('\n')

        tmp = fmyFile.readline().rstrip()

    fSmaller.close()
    fmyFile.close()
    
    if True:
        useMummerAlign(mummerLink, folderName, mummerFile, inputFile + "_contigs_Double.fasta", inputFile + "_contigs_Double.fasta")
        
        
    lengthDic = obtainLength(folderName, inputFile + "_contigs_Double.fasta") 
    
    dataSetRaw = extractMumData(folderName, mummerFile + "Out")
    
    # ## Format [ helperStart, helperEnd , readStart, readEnd,matchLen1,matchLen2,percentMatch,helperName,readName]
    
    
    dataSet = []
    
    for eachitem in dataSetRaw: 
        helperStart, helperEnd , readStart, readEnd, matchLen1, matchLen2, percentMatch, helperName, readName = eachitem 
        
        detailHelper = helperName.split('_')
        detailRead = readName.split('_')
        

        if detailHelper[0] != detailRead[0] and  helperName != readName and max(matchLen1, matchLen2) > minLen and readStart < readEnd  and min(helperStart, readStart) < thres and min(lengthDic[helperName] - helperEnd, lengthDic[readName] - readEnd) + 1 < thres:
            conditionForMatch = True
        else:
            conditionForMatch = False

        if conditionForMatch :
            if helperStart < thres:
                
                dataSet.append((max(matchLen1, matchLen2), readName, helperName))
    
    dataSet.sort(reverse=True)
    
    numberOfContig = len(lengthDic)
    
    return numberOfContig, dataSet


def loggingReadsToRepeat(blockedSet, contigList):
    
    # for eachitem in blockedSet:
    #    print eachitem
    
    matchPair = formMatchPairFromReadInfo(blockedSet)
    matchPair.sort()
    
    print "loggingReadsToRepeat" 
    
    startList = []
    endList = []
    
    returnList = []
    
    for eachitem in contigList:
        startList.append(eachitem[0])
        endList.append(eachitem[-1])
    
    # print startList, endList
    for eachitem in matchPair: 
        mystart = eachitem[0]
        myend = eachitem[1]
        len1, len2 = eachitem[3], eachitem[4]
        
        
        if myend in startList and mystart in endList :
            returnList.append([mystart, myend, len1, len2])
            # print eachitem
    
    returnList.sort()

    return returnList

        
def blockExtraStored(storedStrand, myExtraLinkList, folderName):
    print "blockExtraStored"
    myExtraLinkList.sort()
    f = open(folderName + 'extraConnect.txt', 'w')
    for key, items in groupby(myExtraLinkList, itemgetter(0, 1)):
        mystart = key[0]
        myend = key[1]
        maxLen = -1
        storedLenPair = [-1, -1]
        count = 0
        for eachsub in items:
            print "\t", eachsub
            count += 1
            len1, len2 = eachsub[2], eachsub[3]
            if min(len1, len2) > maxLen:
                maxLen = min(len1, len2)
                storedLenPair = [len1, len2]
        print "Copy count", count, maxLen, storedLenPair
        # if count > 1:
        f.write(str(storedStrand[mystart][0]) + '_' + storedStrand[mystart][1] + ';' + str(storedStrand[myend][0]) + '_' + storedStrand[myend][1] + ';' + str(storedLenPair[0]) + ';' + str(storedLenPair[1]) + '\n')
            
    f.close()
    # print storedStrand 
    # print myExtraLinkList
        

def loadEdgeFromBlockedReads(folderName):
    extraEdges = []
    f = open(folderName + 'extraConnect.txt', 'r')
    tmp = f.readline().rstrip()
    while len(tmp) > 0 :
        myInfo = tmp.split(';')
        inNode = myInfo[0].split('_')
        outNode = myInfo[1].split('_')
        inWt = int(myInfo[2])
        outWt = int(myInfo[3])
        
        myInIndex , myOutIndex = 1997, 1997 
        if inNode[1] == 'p':
            myInIndex = 2 * int(inNode[0])
        else:
            myInIndex = 2 * int(inNode[0]) + 1 
        
        if outNode[1] == 'p':
            myOutIndex = 2 * int(outNode[0])
        else:
            myOutIndex = 2 * int(outNode[0]) + 1 
        
        extraEdges.append([myInIndex, myOutIndex, inWt, outWt])
        tmp = f.readline().rstrip()
    
    f.close()
    print extraEdges
    return extraEdges



def formMatchPairFromReadInfo(dataSet):
    dataSet.sort()
    
    
    matchPair = []
    for key, items in groupby(dataSet, itemgetter(0)):
        print "key", key
        left = []
        right = []
        
        for subitem in items:

            myArr = subitem[1].split('_')
            orientation = myArr[1]

            if orientation == 'p' :
                contigNum = int(myArr[0][6:]) * 2
            else:
                contigNum = int(myArr[0][6:]) * 2 + 1
            
            if subitem[2] == 'L':
                left.append([contigNum, subitem[3]])
            else:
                right.append([contigNum, subitem[3]])

        for eachleft in left:
            for eachright in right:
                leftIndex , rightIndex = eachleft[0], eachright[0]
                leftLen, rightLen = eachleft[1], eachright[1]
                
                if leftIndex != rightIndex and leftIndex / 2 != rightIndex / 2:
                    matchPair.append([rightIndex, leftIndex, min(leftLen, rightLen), rightLen, leftLen, key])

    matchPair.sort()
    return matchPair


###################################################### Key functions

# ## 0) Preprocess by removing embedded contigs (I: contigs.fasta ; O : noEmbed.fasta)


def removeEmbedded(folderName , mummerLink):
    print "removeEmbedded"
    thres = 10
    os.system("sed -e 's/|//g' " + folderName + "contigs.fasta  > " + folderName + "contigs2.fasta")

    os.system("cp " + folderName + "contigs2.fasta " + folderName + "contigs.fasta") 

    if True:
        useMummerAlign(mummerLink, folderName, "self", "contigs.fasta", "contigs.fasta")
    
    dataList = extractMumData(folderName, "selfOut")
    
    dataList = transformCoor(dataList)
    
    lenDic = obtainLength(folderName, 'contigs.fasta')
    
    removeList = []
    for eachitem in dataList:
        match1, match2, name1, name2 = eachitem[4], eachitem[5], eachitem[7], eachitem[8]
        
        if name1 != name2:
            l1, l2 = lenDic[name1], lenDic[name2]
            
            if abs(l1 - match1) < thres and abs(l2 - match2) > thres:
                removeList.append(name1)
            elif abs(l1 - match1) > thres and abs(l2 - match2) < thres:
                removeList.append(name2)
            elif abs(l1 - match1) < thres and abs(l2 - match2) < thres:
                print "Both shortembedd", eachitem
                
    
    
    nameList = []
    for eachitem in lenDic:
        nameList.append(eachitem)

    print len(nameList)
    
    for eachitem in removeList:
        if eachitem in nameList:
            nameList.remove(eachitem)
    print len(nameList)
    
    putListToFileO(folderName, "contigs.fasta", "noEmbed", nameList)


# ## 1) Fetch best successor and best predecessor for each contigs (I: noEmbed.fasta ;  O:  left_connect , right_connect )
def fetchSuccessor(folderName , mummerLink): 
    
    print "fetchSuccessor"
    left_connect, right_connect = [], [] 
        
    print "Direct greedy"
    numberOfContig, dataSet = obtainLinkInfo(folderName, mummerLink, "noEmbed", "greedy")
    # [next_item, overlap_length]
    
    leftConnect = [[-1, -1] for i in range(numberOfContig)]
    rightConnect = [[-1, -1] for i in range(numberOfContig)]
    
    dataSet.sort(reverse=True, key=itemgetter(1))
    
    for key, items in groupby(dataSet, itemgetter(1)):
        # if key == "Contig217_d":
        #    print "dddd"
        maxVal = -1
        myName = key
        connectorName = "" 
        for eachsubitem in items:
            if eachsubitem[0] > maxVal:
                maxVal = eachsubitem[0]
                connectorName = eachsubitem[2]
        

        prefix = myName.split('_')
        suffix = connectorName.split('_')
        lengthOfOverlap = maxVal
        
        if prefix[1] == 'p':
            prefixContig = int(prefix[0][6:]) * 2 
        else:
            prefixContig = int(prefix[0][6:]) * 2 + 1
        
        if suffix[1] == 'p':
            suffixContig = int(suffix[0][6:]) * 2 
        else:
            suffixContig = int(suffix[0][6:]) * 2 + 1
            
        assert(rightConnect[prefixContig][0] == -1)
        rightConnect[prefixContig][0] = suffixContig
        rightConnect[prefixContig][1] = lengthOfOverlap
        

    dataSet.sort(reverse=True, key=itemgetter(2))
    
    for key, items in groupby(dataSet, itemgetter(2)):

        maxVal = -1
        myName = key
        connectorName = "" 
        for eachsubitem in items:
            if eachsubitem[0] > maxVal:
                maxVal = eachsubitem[0]
                connectorName = eachsubitem[1]
        

        prefix = connectorName.split('_')
        suffix = myName.split('_')
        lengthOfOverlap = maxVal
        
        if prefix[1] == 'p':
            prefixContig = int(prefix[0][6:]) * 2 
        else:
            prefixContig = int(prefix[0][6:]) * 2 + 1
        
        if suffix[1] == 'p':
            suffixContig = int(suffix[0][6:]) * 2 
        else:
            suffixContig = int(suffix[0][6:]) * 2 + 1
            
        assert(leftConnect[suffixContig][0] == -1)
        leftConnect[suffixContig][0] = prefixContig 
        leftConnect[suffixContig][1] = lengthOfOverlap
    
    
    # ## Write to file: 
    f = open(folderName + 'rightConnect.txt', 'w')
    for eachitem, dummyIndex in zip(rightConnect, range(len(rightConnect))):
        f.write(str(dummyIndex) + ',' + str(eachitem[0]) + ',' + str(eachitem[1]) + '\n')
        
    f.close()
    
    f = open(folderName + 'leftConnect.txt', 'w')
    for eachitem, dummyIndex in zip(leftConnect, range(len(leftConnect))):
        f.write(str(dummyIndex) + ',' + str(eachitem[0]) + ',' + str(eachitem[1]) + '\n')
        
    f.close()

# ## 2) Form seqGraph (I: left_connect, right_connect ; O: startList, graphNodes )
def formSeqGraph(folderName , mummerLink):
    print "formSeqGraph" 
    startList, graphNodes = [], []
    
    rightConnect = readConnectList(folderName, "rightConnect.txt")
    leftConnect = readConnectList(folderName, "leftConnect.txt")
    
    numberOfNodes = len(rightConnect)
    
    G = seqGraph(numberOfNodes)
        
    for eachitem, i  in zip(rightConnect, range(len(rightConnect))):
        index = i
        connector, weight = eachitem
        G.insertEdge(index, connector, weight)
    
    for eachitem, i  in zip(leftConnect, range(len(leftConnect))):
        index = i
        connector, weight = eachitem
        G.insertEdge(connector, index, weight)
    

    G.cleanEdge()
    G.condense()
    G.saveToFile(folderName, "condensedGraph.txt")
    G.checkSelfLoops()
    G.checkCompleteness()
    
    G2 = seqGraph(0)
    G2.loadFromFile(folderName, "condensedGraph.txt")
    
    compareGraphUnitTest(G, G2)
    G.reportDummyUsefulNode()
    G.reportEdge()
    
    graphFileName = "condensedGraph.txt"
    contigFile = "noEmbed_Double.fasta"
    outContigFile = "improved.fasta"
    outOpenList = "openZone.txt"
    
    readContigOut(folderName, mummerLink, graphFileName, contigFile, outContigFile, outOpenList)
   
    

# ## 3) X-phased seqGraph (I: startList, graphNodes; O: startList, graphNodes )
def xPhased(folderName , mummerLink):
    # ## Repeat resolution  [Proxy for MB]
    # 1. Re-form the contig string graph with ALL connections from contigs only V
    # 2. Log down the reads and associated blocked contigs V 
    # 3. Use reads to connect;
    # 4. Transform graph by identifying 1 successor/predecessor case ; Condense(important);
    # 5. Read out contigs
    
    numberOfContig, dataSet = obtainLinkInfo(folderName, mummerLink, "improved2", "mb")
    
    lenDic = obtainLength(folderName, "improved2_Double.fasta")
    
    confidenLenThres = 0 
    
    G = seqGraph(numberOfContig)
    extraEdges = loadEdgeFromBlockedReads(folderName)
    
    for eachitem in dataSet:
        # print eachitem
        wt, myin, myout = eachitem
        myInData = myin[6:].split('_')
        myOutData = myout[6:].split('_')
        
        if myInData[1] == 'p':
            offsetin = 0
        else:
            offsetin = 1
        
        if myOutData[1] == 'p':
            offsetout = 0
        else:
            offsetout = 1
            
        i = int(myInData[0]) * 2 + offsetin
        j = int(myOutData[0]) * 2 + offsetout
        
        ck = False
        
        for eachedge in extraEdges:
            mystart, myend, len1, len2 = eachedge[0], eachedge[1], eachedge[2] , eachedge[3]
            if [i, j] == [mystart, myend] and min(len1, len2) >= wt and lenDic[myin] >= confidenLenThres and lenDic[myout] >= confidenLenThres:
                ck = True
                
        if ck:
            G.insertEdge(i, j, wt)
    
    
    # G.reportEdge()
    G.MBResolve()
    G.reportEdge()
        
    G.saveToFile(folderName, "condensedGraphMB.txt")
    graphFileName = "condensedGraphMB.txt"
    contigFile = "improved2_Double.fasta"
    outContigFile = "improved3.fasta"
    outOpenList = "openZoneMB.txt"
    
    readContigOut(folderName, mummerLink, graphFileName, contigFile, outContigFile, outOpenList)
    
    
    # ## Repeat resolution  [Proxy for phasing step]
    # 6. Find out the repeat region by MSA
    # 7. Find out the location of SNPs and extend across repeat 
    # [short cut : use contig creator : your job here is to get data into the correct formats]
    
    
    
    
    print "xPhased"



# ## 4) EC reduction (I: startList, graphNodes ; O: startList, graphNodes )
def ECReduction(folderName , mummerLink):
    print "ECReduction" 
    

# ## 5) Read the contigs out (I: startList, graphNodes, ; O:improved.fasta, openZone.txt)
def readContigOut(folderName, mummerLink, graphFileName, contigFile, outContigFile, outOpenList):
    
    print "readContigOut"
    
    G = seqGraph(0)
    G.loadFromFile(folderName, graphFileName)
    G.findStartEndList()
    
    myContigsDic = loadContigsFromFile(folderName, contigFile)
    
    contigUsed = [False for i in range(len(G.graphNodesList) / 2)]
     
    seqToPrint = []
    openList = []
    
    noForRevMismatch = True
    
    for eachnode in G.graphNodesList:

        if len(eachnode.nodeIndexList) > 0:
            tmpSeq = ""
            # ## debug consistency of t/f
            ckList = []
            for dummy in eachnode.nodeIndexList:
                indexToAdd = dummy
                readNum = indexToAdd / 2
                ckList.append(contigUsed[readNum])
                
            if (len(ckList) > 0 and not all(ckList) and any(ckList)):
                noForRevMismatch = False
                
            
            # ## end debug 
            
            for i in range(len(eachnode.nodeIndexList)):
                
                indexToAdd = eachnode.nodeIndexList[i]
                readNum = indexToAdd / 2
                orientation = indexToAdd % 2 
            
                    
                if contigUsed[readNum] == False:

                    if i != len(eachnode.nodeIndexList) - 1:

                        overlapLen = eachnode.overlapList[i]

                        if orientation == 0:
                            tmpSeq = tmpSeq + myContigsDic['Contig' + str(readNum) + '_' + 'p'][0:-overlapLen]
                        else:
                            tmpSeq = tmpSeq + myContigsDic['Contig' + str(readNum) + '_' + 'd'][0:-overlapLen]
                    else:
                        if orientation == 0:
                            tmpSeq = tmpSeq + myContigsDic['Contig' + str(readNum) + '_' + 'p']
                        else:
                            tmpSeq = tmpSeq + myContigsDic['Contig' + str(readNum) + '_' + 'd']
                            
                    contigUsed[readNum] = True
                    
            if len(tmpSeq) > 0:
                if eachnode.nodeIndex in G.myStartList:
                    openList.append('Segkk' + str(len(seqToPrint)) + ',noprev')
                if eachnode.nodeIndex in G.myEndList:
                    openList.append('Segkk' + str(len(seqToPrint)) + ',nonext')
                
                # ## Debug
                if eachnode.nodeIndex == 444:
                    print 439, len(seqToPrint)
                
                if eachnode.nodeIndex == 67 :
                    print 67, len(seqToPrint)
                    
                    
                # ## End Debug
                seqToPrint.append(tmpSeq)
    
    print "No forward/reverse mismatch ?", noForRevMismatch
    fImproved = open(folderName + outContigFile, 'w')
    for eachcontig, dummyIndex in zip(seqToPrint, range(len(seqToPrint))):
        fImproved.write(">Segkk" + str(dummyIndex) + '\n')
        fImproved.write(eachcontig + '\n')
         
    fImproved.close()
    
    print "All contigs used? ", all(contigUsed)
    print "NContig", len(seqToPrint)

    f = open(folderName + outOpenList, 'w')
    f.write(str(len(seqToPrint)) + '\n')
    for eachitem in openList:
        f.write(str(eachitem) + str('\n'))
    f.close()


# ## 6) Fill gap(I: improved.fasta, openZone,txt ; O: improved2.fasta )
def fillGap(folderName , mummerLink):
    print "fillGap"
    # ## Load an extra field V
    # ## migrate the existing code to here  V 
    # ## add it the functionality of checking filtered list V 
    # ## testing 
    
    
    
    formRelatedReadsFile(folderName, mummerLink)
    extractEdgeSet(folderName, mummerLink)
    
    
    
    # os.system("cp raw_reads.part* "+ folderName)
    os.system("rm raw_reads.part*")
    
    # os.system("cp relatedReads_Double.part* "+ folderName)
    os.system("rm relatedReads_Double.part*")


    


# ## 7) Compare with reference (I: improved.fasta, improved2.fasta, reference.fasta ; O : assembly assessment report )
def compareWithReference(folderName , mummerLink):
    print "compareWithReference"
    
    quastEvaluate(folderName, "quast-2.3/", originalName="contigs.fasta", improvedNameList=["noEmbed.fasta", "improved.fasta", "improved2.fasta", "improved3.fasta"] , referenceName="reference.fasta")
    # quastEvaluate(folderName, "quast-2.3/", originalName = "contigs.fasta", improvedNameList= ["noEmbed.fasta", "improved.fasta"] , referenceName= "reference.fasta" )
    


def performMapping(folderName, mummerLink, mapcontigsname):
    print "performMapping"
    info = mapcontigsname.split("_")
    
    oldRef, newRef = info[0], info[1]
    print oldRef, newRef
    command = mummerLink + "nucmer --maxmatch --nosimplify  -p " + folderName + "outMapper " + folderName + newRef + " " + folderName + oldRef
    os.system(command)
    
    command = mummerLink + "show-tiling -v 50 -g 50000 -c " + folderName + "outMapper.delta > " + folderName + "mappingResults.txt"
    os.system(command)
    
    command = "more " + folderName + "mappingResults.txt"
    os.system(command)

    

     

###################################################### Starting point
def mainFlow(folderName , mummerLink, pickupname, mapcontigsname):      
    print "Go Bears! ! !" 
    
    print "pickupname, mapcontigsname", pickupname, mapcontigsname
    
    
    if not pickupname in ["noEmbed.fasta", "improved.fasta", "improved2.fasta"]:
        removeEmbedded(folderName , mummerLink)
    
    if not pickupname in ["improved.fasta", "improved2.fasta"]:
        fetchSuccessor(folderName , mummerLink)
        formSeqGraph(folderName , mummerLink)
        
    if not pickupname in ["improved2.fasta"]:
        fillGap(folderName , mummerLink)

    xPhased(folderName , mummerLink)

    
    
    # ECReduction(folderName , mummerLink )
    # compareWithReference(folderName , mummerLink)
    
    if mapcontigsname != None:
        performMapping(folderName, mummerLink, mapcontigsname)
        
    print "<3 Do cool things that matter <3"

def checkingPath(folderName, mummerLink):
    
    pathExists = True
    newFolderName, newMummerLink = "", ""
    
    if folderName[-1] == "/":
        newFolderName = folderName
    else:
        newFolderName = folderName + "/" 
    
    if mummerLink[-1] == "/":
        newMummerLink = mummerLink
    else:
        newMummerLink = mummerLink + "/" 
    
    if not os.path.isdir(newFolderName):
        pathExists = False
        print "Not exists : " + newFolderName
    
    if not os.path.isdir(newMummerLink):
        pathExists = False
        print "Not exists : " + newMummerLink
    
    if not os.path.exists(newFolderName + "contigs.fasta"):
        pathExists = False
        print "Not exists : " + newFolderName + "contigs.fasta"
    
    if not os.path.exists(newFolderName + "raw_reads.fasta"):
        pathExists = False
        print "Not exists : " + newFolderName + "raw_reads.fasta"
    
    
    return pathExists, newFolderName, newMummerLink
     
# folderName = "S_cerivisea/"
# mummerLink = "MUMmer3.23/"

t0 = time.time()

parser = argparse.ArgumentParser(description='FinisherSC : a repeat-aware tool to upgrade de-novo assembly with long reads')
parser.add_argument('folderName')
parser.add_argument('mummerLink')
parser.add_argument('-p', '--pickup', help='Picks up existing work (input is noEmbed.fasta, improved.fasta or improved2.fasta)', required=False)
parser.add_argument('-o', '--mapcontigs', help='Maps new contigs to old contigs(input is of the format of contigs.fasta_improved3.fasta which means improved3.fasta will be mapped back to contigs.fasta; Output can be found in mappingResults.txt in the destinedFolder;)', required=False)
parser.add_argument('-f', '--fast', help= 'Fast aligns contigs (input is True)', required=False)

args = vars(parser.parse_args())

print "args", args
pathExists, newFolderName, newMummerLink = checkingPath(args['folderName'], args['mummerLink'])

if args['fast'] == "True":
    globalFast = True
else:
    globalFast = False
    

if pathExists:
    mainFlow(newFolderName, newMummerLink, args['pickup'], args['mapcontigs'])
else:
    print "Sorry. The above folders or files are missing. If you continue to have problems, please contact me(Ka-Kit Lam) at kklam@eecs.berkeley.edu"

print  "Time", time.time() - t0

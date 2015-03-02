from itertools import groupby
import abunHouseKeeper

from finisherSCCoreLib import IORobot
from finisherSCCoreLib import graphLib



class seqGraphNodeWt(graphLib.seqGraphNode):
    def __init__(self, nodeIndex):
        graphLib.seqGraphNode.__init__(self, nodeIndex)
        self.nodeWt = 0 
            
        
    def getNextNodesFromIndices(self):
        returnList = []
        for eachitem in self.listOfNextNodes:
            returnList.append(eachitem[0])
        return returnList
    
    def getPrevNodesFromIndices(self):
        returnList = []
        for eachitem in self.listOfPrevNodes:
            returnList.append(eachitem[0])
        return returnList

        
class seqGraphWt(graphLib.seqGraph):
    def __init__(self, numberOfNodes):
        self.graphNodesList = [seqGraphNodeWt(i) for i in range(numberOfNodes)] 
    def formReportName(self, index):
        return str(index) + "_" +str(self.graphNodesList[index].nodeWt)

    def loadFromFile(self, folderName, fileName):
        
        f = open(folderName + fileName, 'r')
        numberOfNodes = 0 
        tmp = f.readline().rstrip()
        while (len(tmp) > 0):
            tmp = f.readline().rstrip()
            numberOfNodes += 1
        f.close()
        
        self.graphNodesList = [seqGraphNodeWt(i) for i in range(numberOfNodes)]
        
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
        
    def findConnectedComponents(self):
        for eachitem in self.graphNodesList:
            eachitem.visited = False 
        connectedList = []
        
        for eachnode in self.graphNodesList:
            if len(eachnode.nodeIndexList) > 0:
                tmpList = []
                
                if eachnode.visited == False:
                    stack = [eachnode]
                    tmpList = [eachnode.nodeIndex]
                    while len(stack) > 0:
                        currenttmp = stack.pop(0)
                        currenttmp.visited = True
                        
                        for eachsubIndex in currenttmp.listOfNextNodes:
                            eachsub = self.graphNodesList[eachsubIndex[0]]
                            if eachsub.visited == False:
                                stack.append(eachsub)
                                tmpList.append(eachsub.nodeIndex)
                    connectedList.append(tmpList)
                                
        
        return connectedList
    
        


def findAllReachable(i, N1, G):
    for eachnode in G.graphNodesList:
        eachnode.visited = False
        
    myList = []
    tmpList = [G.graphNodesList[i]]
    neighborList = []
    while len(tmpList) > 0:
        current = tmpList.pop(0)
        for eachChildIndex in current.listOfNextNodes:
            eachChild = G.graphNodesList[eachChildIndex[0]]
            
            if eachChild.nodeIndex >= N1:
                if eachChild.visited == False:
                    eachChild.visited = True
                    tmpList.append(eachChild)
            else:
                neighborList.append(eachChild)
    
    connectedNodeIndexList = []
    
    for eachitem in neighborList:
        connectedNodeIndexList.append(eachitem.nodeIndex)
    
    connectedNodeIndexList.sort()

    for key, items in groupby(connectedNodeIndexList):
        myList.append(key)
    
    return myList 



def filterEdge(adjacencyList, folderName, contigFilename):
    lenDic = IORobot.obtainLength(folderName, contigFilename + "_Double.fasta")
    thresFoPhase = 2000
    smallList, largeList = [], []
    for eachitem in lenDic:
        id = abunHouseKeeper.parseEdgeNameToID(eachitem, 'C')
        if lenDic[eachitem] < thresFoPhase:
            smallList.append(id)
        else:
            largeList.append(id)
    
    newAdjacencyList = [[] for i in range(len(adjacencyList))]
    
    for i in largeList:
        for eachitem in adjacencyList[i]:
######## IMPORTANT:
            if  eachitem in largeList and eachitem / 2 != i / 2:
######## NEED TO REMOVE IN PRODUCTION if True
                newAdjacencyList[i].append(eachitem)
    
    
    print "len(smallList)  , len(largeList): ", len(smallList)  , len(largeList)
    print "lenDic: ", lenDic
    
    for eachitem in newAdjacencyList:
        print "newAdjacencyList :", eachitem 
        
    return newAdjacencyList
        

def addIndicesToReachable(x, G, N1):
    
    
    tmp = int(x.split('_')[0])

    findAllReachable(tmp , N1, G)
    
    for eachnode in G.graphNodesList:
        if eachnode.visited == True and eachnode.nodeIndex >= N1:
            eachnode.visitLabelList.append(x)


def markReachableIndices(G, Grev, rIn, rOut, N1):
    print "markReachableIndices"
    for eachitem in G.graphNodesList:
        eachitem.visitLabelList = []
        eachitem.visited = False
    
    for eachitem in Grev.graphNodesList:
        eachitem.visitLabelList = []
        eachitem.visited = False
    
    for x in rIn:
        addIndicesToReachable(x, G, N1)
    
    # # Bug : wrong indexing 

    for y in rOut:
        addIndicesToReachable(y, Grev, N1)
        
    for i in range(len(Grev.graphNodesList)):
        if len(Grev.graphNodesList[i].visitLabelList) > 0:
            G.graphNodesList[i].visitLabelList = G.graphNodesList[i].visitLabelList + Grev.graphNodesList[i].visitLabelList  
    
    
        
        

def formReverseGraph(G):
    nNode = len(G.graphNodesList)
    Grev = seqGraphWt(nNode)
    for i in range(nNode):
        for j in range(nNode):
            haveInserted = graphLib.nameInEdgeList(j, G.graphNodesList[i].listOfNextNodes) 
            if haveInserted:      
                Grev.insertEdge(j, i, 100)
    return Grev


def markInsideNodes(G, kkIn, kkOut):

    print "markInsideNodes"
    singleMissList = [[int(kkIn[0].split('_')[0])] , [int(kkIn[1].split('_')[0])], [int(kkOut[0].split('_')[0])], [int(kkOut[1].split('_')[0])]]
    allPassList = []
    
    for eachitem in G.graphNodesList:
        result = checkMiss(eachitem, kkIn, kkOut)
        if result == -1: 
            allPassList.append(eachitem.nodeIndex)
        elif result == 0:
            singleMissList[1].append(eachitem.nodeIndex)
        elif result == 1:
            singleMissList[0].append(eachitem.nodeIndex)
        elif result == 2:
            singleMissList[3].append(eachitem.nodeIndex)
        elif result == 3:
            singleMissList[2].append(eachitem.nodeIndex)
        
        # elif 0 <= result <= 3:
        #    singleMissList[result].append(eachitem.nodeIndex)
        
    return singleMissList, allPassList



def checkMiss(myNode, rIn, rOut):
    index = 1997
    myList = myNode.visitLabelList

    combineCheckList = rIn + rOut
    if set(myList) == set(combineCheckList):
        return -1
    for i in range(4):
        if set(myList) == set(combineCheckList[0:i] + combineCheckList[i + 1:]):
            return i
            break
    return index


 
def markStartEndNodes(G, rIn, rOut, singleMissList, allPassList):
    print "markStartEndNodes"
    myStartIndex, myEndIndex = -1, -1
    
    print singleMissList 
    print allPassList
    for i in allPassList:
        nextList = G.graphNodesList[i].getNextNodesFromIndices()
        prevList = G.graphNodesList[i].getPrevNodesFromIndices()
        
        counter = 0 
        for j in range(2):
            if len(set(prevList).intersection(set(singleMissList[j]))) > 0:
                counter = counter + 1
        
        
        for j in range(2, 4):
            if len(set(nextList).intersection(set(singleMissList[j]))) > 0:
                counter = counter + 10
        
        if counter / 10 == 2 and counter % 10 == 0:
            myEndIndex = i
        
        if counter / 10 == 0 and counter % 10 == 2:
            myStartIndex = i
            
        if counter == 22:
            myStartIndex = i
            myEndIndex = i
            break
    
    return myStartIndex, myEndIndex
            

def BFS(x, y, G, N1):
    for eachnode in G.graphNodesList:
        eachnode.visited = False
        eachnode.backPtr = -1
        
    tmpList = [G.graphNodesList[x]]
    starter = True
    while len(tmpList) > 0:
        current = tmpList.pop(0)
        if current.nodeIndex >= N1 or starter:
            starter = False
            for i in current.getNextNodesFromIndices():
                if G.graphNodesList[i].visited == False:
                    G.graphNodesList[i].visited = True
                    G.graphNodesList[i].backPtr = current.nodeIndex
                    tmpList.append(G.graphNodesList[i])
            
        
    # print "end BFS", y ,x
    z = y
    tmpList.append(y)           
    while z != x:
        z = G.graphNodesList[z].backPtr
        tmpList.append(z)
    
    # print "end Back Ptr"

    return tmpList[::-1]

def markInterior(G , myStartIndex, myEndIndex, N1):
    print "markInterior"
    myPathway = BFS(myStartIndex, myEndIndex , G, N1)
    return myPathway

   
def markFlankingRegion(G, rIn, rOut, myStartIndex, myEndIndex, N1):
    
    print "markFlankingRegion"
    print myStartIndex, myEndIndex
    myPathwayList = [[] for i in range(4)]
    
    for i in range(2):
        myPathwayList[i] = BFS(rIn[i], myStartIndex, G, N1)
        myPathwayList[i + 2] = BFS(myEndIndex, rOut[i], G, N1)
        
    return myPathwayList

def DFS(G, x, N1, startIndex, endIndex, mypath):
    
    if x >= N1 or x == startIndex:
        if G.graphNodesList[x].visited == False:
            G.graphNodesList[x].visited = True
            for eachChild in G.graphNodesList[x].listOfNextNodes:
                DFS(G, eachChild[0], N1, startIndex, endIndex, mypath + [x])
    elif x == endIndex:
        print mypath + [x]
        
        
def recCheck(G, i , N1, counter, contigList):
    if counter > 0:
        for eachitem in G.graphNodesList[i].listOfNextNodes:
            if G.graphNodesList[eachitem[0]].nodeIndex < N1:
                contigList.append(eachitem[0])
                continue
            recCheck(G, eachitem[0] , N1, counter - 1, contigList)
            
def checkFourHoppers(G, i , N1):
    contigList = []
    recCheck(G, i , N1, 3, contigList)
    
    returnList = []
    contigList.sort()
    for key, items in groupby(contigList):
        returnList.append(key)
    return returnList
                        
def debugGraphPath(startIndex, endIndex, G, N1):
    for eachnode in G.graphNodesList:
        eachnode.visited = False
        
    print startIndex , checkFourHoppers(G, startIndex , N1)

    
def markAssociatedReads(G, singleMissList, allPassList):
    print "markAssociatedReads"
    flankingList = [[] for i in range(4)]
    repeatList = []
    for eachitem in G.graphNodesList:
        nextList = eachitem.getNextNodesFromIndices()
        prevList = eachitem.getPrevNodesFromIndices()
        
        
        if len(set(nextList).intersection(set(allPassList))) > 0 :
            nextInt = True
        else:
            nextInt = False
        
        if len(set(prevList).intersection(set(allPassList))) > 0 :
            prevInt = True
        else:
            prevInt = False
        
        
        missArr = [False for i in range(4)]
        for j in range(2):
            if len(set(prevList).intersection(set(singleMissList[j]))) > 0 :
                missArr[j] = True
            if len(set(nextList).intersection(set(singleMissList[j + 2]))) > 0:
                missArr[j + 2] = True
        
        for j in range(2):
            if missArr[j] and prevInt:
                flankingList[j].append(eachitem.nodeIndex)
            elif missArr[j + 2] and nextInt:
                flankingList[j + 2].append(eachitem.nodeIndex)
            elif prevInt and nextInt:
                repeatList.append(eachitem.nodeIndex)
        
    return flankingList, repeatList
            

def checkPathLength(path, G, N1, folderName):
    
    lenDicRR = IORobot.obtainLength(folderName, "phasingSeedName_Double.fasta")
    sumLength = 0
    overlapLength = 0
    for index, i in zip(path, range(len(path))):
        header = "Read" + str((index - N1) / 2) + "_"
        if (index - N1) % 2 == 0:
            header = header + "p"
        else:
            header = header + "d"
        print "lenDicRR[header], ", lenDicRR[header], header 
        print (index - N1) * 2 + 1, (index - N1) * 2 + 2
        sumLength = sumLength + lenDicRR[header]
        
        if i != len(path) - 1:
            for eachnext in G.graphNodesList[index].listOfNextNodes:
                if eachnext[0] == path[i + 1]:
                    overlapLength = overlapLength + eachnext[1]
                    break 
    print sumLength, overlapLength, sumLength - overlapLength
                


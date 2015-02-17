from operator import itemgetter


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
            
    return newList


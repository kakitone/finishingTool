import commonLib
from itertools import groupby
from operator import itemgetter
import os
import json
import cleaner
import extender
import common 
import sys
import time

# ## Debug check 
def checkGraphLength(G, N1, lenDicRR):
    for eachitem in G.graphNodesList:
        for eachnext in eachitem.listOfNextNodes:
            index, wt = eachnext[0], eachnext[1]
            if index >= N1 and eachitem.nodeIndex >= N1:
                header = "Read" + str((index - N1) / 2) + "_"
                if (index - N1) % 2 == 0:
                    header = header + "p"
                else:
                    header = header + "d"
                
                header2 = "Read" + str((eachitem.nodeIndex - N1) / 2) + "_"
                    
                if (eachitem.nodeIndex - N1) % 2 == 0:
                    header2 = header2 + "p"
                else:
                    header2 = header2 + "d"
                    
                if lenDicRR[header] < wt:
                    print lenDicRR[header], lenDicRR[header2], wt
                    print header, header2
                    
def checkPathLength(path, G, N1, folderName):
    
    lenDicRR = commonLib.obtainLength(folderName, "phasingSeedName_Double.fasta")
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
                
        
# ## Helper Functions
def outputResults(folderName, mummerLink, toPhaseList, N1, G):
    '''    
    Algorithm :
    a) Write as contigs 
    b) Add back reverse complement 
    c) Create G2 as the readOut part 
    d) Output the contigs by a function call

    '''
    # a) 
    combinedName = "contigAndRead_Double.fasta"
    os.system("cp " + folderName + "improved3_Double.fasta " + folderName + combinedName)
    
    fout = open(folderName + combinedName, 'a')
    fin = open(folderName + "phasingSeedName_Double.fasta", 'r')

    tmp = fin.readline().rstrip()
    while len(tmp) > 0:
        if tmp[0] != ">":
            fout.write(tmp + "\n")
        else:
            infoArr = tmp[5:].split("_")
            fout.write(">Contig" + str(int(infoArr[0]) + N1 / 2))
            fout.write("_" + infoArr[1] + "\n")
        tmp = fin.readline().rstrip()
        
    fin.close()
    fout.close()

    # b)
    '''
    [28], [[2, 690, 28], [6, 126, 28], [28, 212, 0], [28, 216, 4]], 1
    
    [2 , 690, 28, 212, 0]
    '''
    
    completePhaseList = []
    for eachitem in toPhaseList:
        repeat = eachitem[-3]
        flanking = eachitem[-2]
        result = eachitem[-1]
        
        
        revrepeat = []
        for eachsub in eachitem[-3][-1::-1]:
            revrepeat.append(eachsub + pow(-1, eachsub))
            
        revflanking = [[] for i in range(4)] 
        
        for j in range(2):
            for eachsub in eachitem[-2][j + 2][-1::-1]:
                revflanking[j].append(eachsub + pow(-1, eachsub))
            for eachsub in eachitem[-2][j][-1::-1]:
                revflanking[j + 2].append(eachsub + pow(-1, eachsub))
            
        revresult = eachitem[-1]
        
        completePhaseList.append([repeat, flanking, result])
        completePhaseList.append([revrepeat, revflanking, revresult])
    
    print "completePhaseList", completePhaseList
    # c) 
    G2 = commonLib.seqGraph(N1)
    nameDic = {}
    for i in range(N1):
        nameDic[i] = i
        
    for eachitem in completePhaseList:
        repeat, flanking, result = eachitem[0] , eachitem[1] , eachitem[2]
        path = [[], []]
        
        if result == 0:
            path[0] = flanking[0][0:-1] + repeat + flanking[2][1:]
            path[1] = flanking[1][0:-1] + repeat + flanking[3][1:]
        else:
            path[0] = flanking[0][0:-1] + repeat + flanking[3][1:]
            path[1] = flanking[1][0:-1] + repeat + flanking[2][1:]
        
        print path[0] , path[1]
        for i  in range(2):
            eachpath = path[i]
            currentNode = G2.graphNodesList[eachpath[0]]
            
            for nextNodeIndex, ctr in zip(eachpath[1:], range(len(eachpath[1:]))):
                if ctr != len(eachpath[1:]) - 1:
                    myindex = len(G2.graphNodesList)
                    nameDic[myindex] = nextNodeIndex
                    
                    newNode = commonLib.seqGraphNode(myindex)
                    G2.graphNodesList.append(newNode)
                else:
                    newNode = G2.graphNodesList[nextNodeIndex]
                    
                wt = 0
                for eachck in G.graphNodesList[nameDic[currentNode.nodeIndex]].listOfNextNodes:
                    if eachck[0] == nextNodeIndex:
                        wt = eachck[1]
                        break
                    
                newNode.listOfPrevNodes.append([currentNode.nodeIndex, wt])
                currentNode.listOfNextNodes.append([newNode.nodeIndex, wt])
                
                currentNode = newNode
                
    graphFileName = "phaseGraphFinal"
    G2.condense()
    G2.saveToFile(folderName, graphFileName)
    
    commonLib.readContigOut(folderName, mummerLink, graphFileName, combinedName, "improved4.fasta", "outOpenListphaing", nameDic)


def filterReverseComp(loadData, N1):
    newLoadData = []
    isUsedSeekRight = [False for i in range(N1)]
    isUsedSeekLeft = [False for i in range(N1)]
    
    for eachitem in loadData:
        s = [-1, -1]
        e = [-1, -1]
        s[0] , s[1], e[0], e[1] = eachitem[-1][0][0], eachitem[-1][1][0], eachitem[-1][2][-1], eachitem[-1][3][-1]
        check = True
        print s, e 
        for eachsub in s:
            if isUsedSeekRight[eachsub] == True:
                check = False
        for eachsub in e:
            if isUsedSeekLeft[eachsub] == True:
                check = False
        if check:
            newLoadData.append(eachitem)
            for eachsub in s:
                isUsedSeekRight[eachsub] = True
                isUsedSeekLeft[eachsub + pow(-1, eachsub % 2)] = True
                
            for eachsub in e:
                isUsedSeekLeft[eachsub] = True
                isUsedSeekRight[eachsub + pow(-1, eachsub % 2)] = True
        
   
    return newLoadData  

def formMap(i, j, wt, lenDicCR, dicFromOriginal, dicToOriginal, N1):
    # shortToLongMap : indexlong    indexshort    longstart    longend    shortstart    shortend   
    # print lenDicCR
    header = "Read" + str((i - N1) / 2) + "_"
    if (i - N1) % 2 == 0 :
        header = header + 'p'
    else:
        header = header + 'd' 
    
    indexlong, indexshort = dicFromOriginal[i], dicFromOriginal[j]
    longstart, longend, shortstart, shortend = lenDicCR[header] - wt, lenDicCR[header] - 1 , 0, wt - 1
    shortToLongMap = [ indexlong   , indexshort  , longstart   , longend    , shortstart  , shortend  ] 
    altMap = [  indexshort  , indexlong   , shortstart  , shortend, longstart   , longend      ] 
    return [shortToLongMap, altMap ]


def transformReadFormat(inputRead):
    # May need to change dependig on what cleaner is looking for... 
    # need to change cleaner here as well
    # ACGT = 1234 
    
    itemIndices = ['A', 'C', 'G', 'T', 'N']
    outputRead = []
    for x in inputRead:
        for i in range(len(itemIndices)):
            if x == itemIndices[i]:
                outputRead.append(i + 1)

    return outputRead

def reformatNoisyReads(folderName, flankingList, repeatList, N1):
    noisyReads, tranDic = [], [] 
    noisyReadsIndiceList = []
    combinedIndices = repeatList 
    
    f = open(folderName + "phasingSeedName_Double.fasta", 'r')
    tmp = f.readline().rstrip()
    readData = []
    while len(tmp) > 0:
        if tmp[0] != '>':
            readData.append(tmp)
        tmp = f.readline().rstrip()
        
    f.close()
    
    
    for eachitem in flankingList:
        combinedIndices = combinedIndices + eachitem
    
    combinedIndices.sort()
    
    for key, items in groupby(combinedIndices):
        if key >= N1:
            noisyReadsIndiceList.append(key)
    
    dicToOriginal, dicFromOriginal = {}, {}
    
    counter = 0
    for i in noisyReadsIndiceList:
        # print i, len(readData)
        noisyReads.append(transformReadFormat(readData[i - N1]))
        dicToOriginal[counter] = i 
        dicFromOriginal[i] = counter  
        counter = counter + 1
        
    return noisyReads, dicToOriginal, dicFromOriginal
    
def reformatToProcessList(folderName , flankingList, repeatList, dicFromOriginal, N1):
    assert(1 == 1)
    toProcessList = [[] for i in range(5)]
    for i in range(len(flankingList)):
        
        for eachsubitem in flankingList[i]:
            if eachsubitem >= N1:
                toProcessList[i].append(dicFromOriginal[eachsubitem])
    
    for eachsubitem in repeatList:
        if eachsubitem >= N1:
            toProcessList[4].append(dicFromOriginal[eachsubitem])
    
    return toProcessList 
    
    
def formShortToLongMapping(folderName, G, toProcessList, dicFromOriginal, dicToOriginal, lenDicCR, N1):
    combinedList = []
    for eachitem in toProcessList:
        combinedList = combinedList + eachitem
    
    combinedList.sort()
    distinctCombined = []
    mappingList = []
    
    for key, items in groupby(combinedList):
        distinctCombined.append(key)

    for i in distinctCombined:
        id = dicToOriginal[i]
        for eachprev in G.graphNodesList[id].listOfPrevNodes:
            if eachprev[0] in dicFromOriginal :
                mappingList = mappingList + formMap(eachprev[0], id, eachprev[1], lenDicCR, dicFromOriginal, dicToOriginal, N1)
        for eachnext in G.graphNodesList[id].listOfNextNodes:
            if eachnext[0] in dicFromOriginal:
                mappingList = mappingList + formMap(id, eachnext[0], eachnext[1], lenDicCR, dicFromOriginal, dicToOriginal, N1)
    
    return mappingList
    
def createIndelRobot(folderName):
    indelRobot = common.parameterRobot()
    indelRobot.defaultFolder = folderName
    # indelRobot.setReadStat( Nshort= parameterRobot.N, Nlong=  parameterRobot.N, Lshort= parameterRobot.L, Llong= parameterRobot.L, p= parameterRobot.p , longOnly = True)
    # indelRobot.setGenomeStat(G = parameterRobot.G, lrep=500, lsnp=200, lint=50 )
    indelRobot.setThresholdPara(liid=30, thresForRandom=0.5, thresForins=0.4, thresFordel=0.4, insMin=4, delMin=4, thresholdForSupport=0.15, subthreshold=9, editsub=-10, editins=-1, editdel=-1, editmatch=1, lookRange=15)
    indelRobot.tunePara()
    indelRobot.snprate = 0.01
    
    return indelRobot
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
            
    
def markFlankingRegion(G, rIn, rOut, myStartIndex, myEndIndex, N1):
    
    print "markFlankingRegion"
    print myStartIndex, myEndIndex
    myPathwayList = [[] for i in range(4)]
    
    for i in range(2):
        myPathwayList[i] = BFS(rIn[i], myStartIndex, G, N1)
        myPathwayList[i + 2] = BFS(myEndIndex, rOut[i], G, N1)
        
    return myPathwayList

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
            
        #    myStartIndex = 
        #    if G.graphNodesList[i].listOfPrevNodes in :
        #    myEndIndex = 
        


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
    Grev = commonLib.seqGraph(nNode)
    for i in range(nNode):
        for j in range(nNode):
            haveInserted = commonLib.nameInEdgeList(j, G.graphNodesList[i].listOfNextNodes) 
            if haveInserted:      
                Grev.insertEdge(j, i, 100)
    return Grev
                
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


# ## Debug code

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

    

# ## end debug code 


def filterEdge(adjacencyList, folderName, contigFilename):
    lenDic = commonLib.obtainLength(folderName, contigFilename + "_Double.fasta")
    thresFoPhase = 2000
    smallList, largeList = [], []
    for eachitem in lenDic:
        id = parseEdgeNameToID(eachitem, 'C')
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
        

def addDataToList(dataList, G, startIndex1, startIndex2, type1, type2):
    threshold = 40
    for eachitem in dataList:
        wt = min(eachitem[4] , eachitem[5])
        
        if eachitem[0] < threshold:

            j = parseEdgeNameToID(eachitem[-2], type1) + startIndex1
            i = parseEdgeNameToID(eachitem[-1], type2) + startIndex2
        else:
            j = parseEdgeNameToID(eachitem[-1], type2) + startIndex2
            i = parseEdgeNameToID(eachitem[-2], type1) + startIndex1
        
        if i == 5 and j == 3279 : 
            print i, j , eachitem
        if i == 3279 and j == 4 :
            print i, j , eachitem
        
        G.insertEdge(i, j, wt) 

def parseEdgeNameToID(name, mytype):

    if mytype == 'C':
        dataInfo = name[6:].split('_')
    elif mytype == 'R':
        dataInfo = name[4:].split('_')
        
    contigNum = int(dataInfo[0])
    
    if dataInfo[1] == 'p':
        id = contigNum * 2
    elif dataInfo[1] == 'd':
        id = contigNum * 2 + 1
        
    return id 

def checkSatisfy(eachitem, lenDic):
    # "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    threshold = 40
     
    seedName = eachitem[-1]
    detectedName = eachitem[-2]
    
    # seedIndices = seedName.split("/")[-1].split("_")
    # detectIndices = detectedName.split("/")[-1].split("_")
    # lenSeed = int(seedIndices[1]) - int(seedIndices[0])
    # lenDetect = int(detectIndices[1]) - int(detectIndices[0])
    
    lenSeed = lenDic[seedName]
    lenDetect = lenDic[detectedName]
    
    checkSeed = False
    checkDetect = False 
    
    if eachitem[0] < threshold or eachitem[1] > lenSeed - threshold: 
        checkSeed = True
    
    if min(eachitem[2], eachitem[3]) < threshold or max(eachitem[2], eachitem[3]) > lenDetect - threshold:
        checkDetect = True
    
    if checkSeed and checkDetect: 
        return True
    else:
        return False

def headTailMatch(eachitem, lenDic):
    start1, end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]
    l1, l2 = lenDic[eachitem[-2]], lenDic[eachitem[-1]]
    thres = 40 
    
    
    diffContig , forwardStrand, headtailoverlap = False, False , False
    
    if eachitem[-2] == eachitem[-1]:
        diffContig = False
    else:
        diffContig = True

    if start2 > end2 : 
        forwardStrand = False
    else:
        forwardStrand = True
    
    
    if (start1 <= thres and end2 >= l2 - thres) or (end1 >= l1 - thres and start2 <= thres): 
        headtailoverlap = True
    else:
        headtailoverlap = False
        
    
    if diffContig and forwardStrand and headtailoverlap:
        return True
    else:
        return False
    
    
    return True

def filterData(dataList, lenDic):
    newDataList = []
    for eachitem in dataList:
        if headTailMatch(eachitem, lenDic):
            newDataList.append(eachitem)

    return newDataList

# ## New Phasing Step
# ## Pre-processing of the phasing 

def getAllAssociatedReads(folderName, mummerLink,forFastaName):
    '''
    Input : relatedReads.fasta, raw_reads.fasta 
    Output : all_associated_reads.fasta
    
     Algorithm : 
        a) Get all the associated reads
        b) Loop for N=1 times : ==> this correspond 4 reads to link between the bridge in total
            i) Align the raws and tmp_seedReads
            ii) Put the new reads into the SeedReads
    '''
    header, referenceFile, queryFile = "seedReads", forFastaName + ".fasta" , "raw_reads.fasta"
    command = "cp " + folderName + "relatedReads.fasta " + folderName + referenceFile
    os.system(command)
    N = 1
    
    for trial in range(N):
        print "trial", trial
        if True:
            command = mummerLink + "nucmer --maxmatch --nosimplify -p " + folderName + header + " " + folderName + referenceFile + " " + folderName + queryFile
            os.system(command)
            
            command = mummerLink + "show-coords -r " + folderName + header + ".delta > " + folderName + header + "Out"
            os.system(command)
        
        dataList = commonLib.extractMumData(folderName, header + "Out")
        filterList = []
        
        lenDicRR = commonLib.obtainLength(folderName, queryFile)
        
        print "len(dataList)", len(dataList)
        for eachitem in dataList:
            if checkSatisfy(eachitem, lenDicRR):
                filterList.append(eachitem)
            
        filterList.sort(key=itemgetter(-1))
        newReads = []
        
        for key, items in groupby(filterList, itemgetter(-1)):
            newReads.append(key)
                                    
        
        f = open(folderName + forFastaName + ".txt", 'w')
        
        for eachitem in newReads:
            f.write(eachitem + "\n")
        f.close()
            
        command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' " + folderName + forFastaName + ".txt " + folderName + "raw_reads.fasta > " + folderName + forFastaName + ".fasta"
        os.system(command)
        
    
def formReadContigStringGraph(folderName, mummerLink, contigFilename, readsetFilename, optTypeFileHeader, graphName):
    
    '''
    Input : all_associated_reads.fasta, improved3.fasta
    Output : (G) String Graph linking the reads and contigs
    Algorithm: 
        a) Form double reads and contigs                            V
        b) Mummer the data and extract dataList three times         V
        c) Use the subroutine to output a graph                     V
        d) Output the graph to a file phasing_String_graph.graph    V
    '''

    G = []

    commonLib.writeToFile_Double1(folderName, contigFilename + ".fasta", contigFilename + "_Double.fasta", "contig")
    commonLib.writeToFile_Double1(folderName, readsetFilename + ".fasta", readsetFilename + "_Double.fasta", "reads")
    
    
    header, referenceFile, queryFile = optTypeFileHeader + "CC", contigFilename + "_Double.fasta" , contigFilename + "_Double.fasta"
    if True:
        commonLib.useMummerAlign(mummerLink, folderName, header, referenceFile, queryFile)

    lenDicCC = commonLib.obtainLength(folderName, contigFilename + "_Double.fasta")
    dataListCC = commonLib.extractMumData(folderName, header + "Out")
    dataListCC = filterData(dataListCC, lenDicCC)
    
    header, referenceFile, queryFile = optTypeFileHeader + "RR", readsetFilename + "_Double.fasta" , readsetFilename + "_Double.fasta"
    if True:
        commonLib.useMummerAlign(mummerLink, folderName, header, referenceFile, queryFile)
    
    lenDicRR = commonLib.obtainLength(folderName, readsetFilename + "_Double.fasta")

    dataListRR = commonLib.extractMumData(folderName, header + "Out")
    dataListRR = filterData(dataListRR, lenDicRR)

    header, referenceFile, queryFile = optTypeFileHeader + "CR", contigFilename + "_Double.fasta" , readsetFilename + "_Double.fasta"
    if True:
        commonLib.useMummerAlign(mummerLink, folderName, header, referenceFile, queryFile)
    
    lenDicCR = dict(lenDicCC.items() + lenDicRR.items())
    dataListCR = commonLib.extractMumData(folderName, header + "Out")
    dataListCR = filterData(dataListCR, lenDicCR)
            
    numberOfNodes = len(lenDicCR) 
    G = commonLib.seqGraph(numberOfNodes)
    N1, N2 = len(lenDicCC), len(lenDicRR)
    print "N1, N2, numberOfNodes: ", N1, N2, numberOfNodes
    
    '''
    e.g. of dataListCC[0], dataListRR[0], dataListCR[0]
    
    [1, 520, 2913194, 2913716, 520, 523, 99.05, 'Contig0_d', 'Contig2_d']
    [1, 1383, 1253, 2603, 1383, 1351, 82.39, 'Read0_d', 'Read1705_p']
    [1, 718, 4334, 5074, 718, 741, 91.91, 'Contig0_d', 'Read1018_d']
    
    '''
    
    # print dataListCC[0]
    # print dataListRR[0]
    # print dataListCR[0]
    
    # for eachitem in dataListCC:
    #    print eachitem
    addDataToList(dataListCC, G, 0, 0, 'C', 'C')
    # for eachitem in dataListRR[0:10]:
    #    print eachitem , lenDicRR[eachitem[-2]], lenDicRR[eachitem[-1]]
    
    
    addDataToList(dataListRR, G, N1, N1, 'R', 'R')
    
    addDataToList(dataListCR, G, 0, N1, 'C', 'R')
    # G.reportEdge()
    G.saveToFile(folderName, graphName)
    
    checkGraphLength(G, N1, lenDicRR)
    
    # print len(G.graphNodesList[0].listOfPrevNodes), len(G.graphNodesList[0].listOfNextNodes)
    print "len(G.graphNodesList)", len(G.graphNodesList)
    
    
def getDistinct(myList):
    newList = [] 
    myList.sort()
    
    for key, items in groupby(myList):
        newList.append(key)
    
    return newList
def identifyRepeat(folderName, mummerLink,contigFilename,contigReadGraph, repeatFilename, optionToRun  ):
    '''
    Input : Graph --- phaseStringGraph1
    Output: repeat pairs { [ (1,2), (3,4) ] , [(5,6),(7,8)] } 
    Algorithm: 
        a) Reachability test on the graph to find the partners
        b) Form Bipartite graph
        c) Find connected component in the bipartite and define as repeat pairs

    '''
    
    # ## (a) reachability test to find partners 
    G = commonLib.seqGraph(0)
    G.loadFromFile(folderName, contigReadGraph)
    # G.reportEdge()
    lenDicCC = commonLib.obtainLength(folderName, contigFilename+"_Double.fasta")
    
    adjacencyList = [[] for i in range(len(lenDicCC))]
    
    N1 = len(lenDicCC)
    
    
    # # Debug
    # for i in range(14):
    #    debugGraphPath(i, 2, G, N1)
    # # End Debug
    
    for i in range(len(lenDicCC)):
        adjacencyList[i] = findAllReachable(i, N1, G) 
        print "i, adjacencyList[i] : ", i , adjacencyList[i]
    
    # ## (b) formation of bipartite graph
    if optionToRun == "tandem" :
        newAdjacencyList = adjacencyList
    elif optionToRun == "xphase": 
        newAdjacencyList = filterEdge(adjacencyList, folderName, contigFilename)
    
    G2 = commonLib.seqGraph(N1 * 2)
    for i in range(N1):
        for j in newAdjacencyList[i]:
            G2.insertEdge(2 * i, 2 * j + 1, 1)
            G2.insertEdge(2 * j + 1, 2 * i, 1)

    clusters = G2.findConnectedComponents()
    
    repeatList = []
    for eachitem in clusters:
        leftList, rightList = [], []
        for eachsubitem in eachitem:
            if eachsubitem % 2 == 0 :
                leftList.append(eachsubitem)
            else:
                rightList.append(eachsubitem)
                
        
        repeatList.append([getDistinct(leftList), getDistinct(rightList)])
           
    with open(folderName + repeatFilename, 'w') as outfile:
        json.dump(repeatList, outfile)

    
    json_data = open(folderName + repeatFilename, 'r')
    loadData = json.load(json_data)
    
    
    assert(loadData == repeatList)
    
def defineRepeatAndFlanking(folderName, mummerLink,contigFilename,contigReadGraph,repeatFilename,repeatSpec ):
    '''
    Input : 
V        a) String graph : G                
V        b) Repeat Pairing : repeatList     
        
    Output : 
V        a) chain of repeat indices (e.g. [S= R1, R33, R45, R24= E]) 
V        b) chain of flanking region indices for in1/2 out1/2 middle (e.g. [C1, R2, R4] )
V        c) in1/2 out1/2 and middle reads per repeat (e.g. [R1, R33, R45, R24])  
        
    Algorithm : 
V        1. Find repeat by graph operations
V        2. Find flanking region by graph operations
V        3. Find associated reads by graph operations
    '''
    
    print "defineRepeatAndFlanking: "


    
    
    # 0. Load previous data
    G = commonLib.seqGraph(0)
    G.loadFromFile(folderName, contigReadGraph)
    Grev = formReverseGraph(G)
    
    json_data = open(folderName + repeatFilename, 'r')
    repeatList = json.load(json_data)
    
    lenDicCC = commonLib.obtainLength(folderName, contigFilename+"_Double.fasta")
    N1 = len(lenDicCC)
    
    
    print "repeatList: ", repeatList
    print "len(G.graphNodesList)", len(G.graphNodesList)
     
    bigDumpList = []
    
    print "len(repeatList)", len(repeatList) , repeatList
    for r in repeatList:
        rIn, rOut = [], []
        for eachitem in r[0]:
            rIn.append(eachitem / 2)
        for eachitem in r[1]:
            rOut.append(eachitem / 2)
        
        if ( len(rIn) == 2 and len(rOut) == 2) or (len(rIn) == 1 and len(rOut) == 1):
            print rIn, rOut
            if  (len(rIn) == 1 and len(rOut) == 1):
                rIn = [rIn[0], rIn[0]]
                rOut = [rOut[0], rOut[0]]
            
            # 1. Records reachable indices
            kkIn , kkOut = [], []
            for eachkk in rIn:
                kkIn.append(str(eachkk)+"_"+"in")
            
            for eachkk in rOut:
                kkOut.append(str(eachkk)+"_"+"out")
                
            
            markReachableIndices(G, Grev, kkIn, kkOut, N1)
            
            # 2. Marks inside nodes
            singleMissList, allPassList = markInsideNodes(G, kkIn, kkOut)
            for i in range(4): 
                print "len(singleMissList[i]), len(allPassList)", len(singleMissList[i]), len(allPassList)

            # 3. Finds start/end of repeat
            myStartIndex, myEndIndex = markStartEndNodes(G, rIn, rOut, singleMissList, allPassList)
            print myStartIndex, myEndIndex
            
            # 4. Find repeat interior by shortest path joining S/E
            repeatPathway = markInterior(G , myStartIndex, myEndIndex, N1)
            print "repeatPathway", repeatPathway
            #checkPathLength(repeatPathway, G, N1, folderName)
            
            # 5. Find flanking region by shortest path search again
            flankingPathsList = markFlankingRegion(G, rIn, rOut, myStartIndex, myEndIndex, N1)
            print flankingPathsList
            
            # 6. Find associated reads by graph node query
            flankingList, repeatList = markAssociatedReads(G, singleMissList, allPassList)
            
            # ## Experimental
            repeatList = allPassList
            
            # ## End Experimental
            for eachlist in flankingList:
                print len(eachlist), len(repeatList)
            
            bigDumpList.append([flankingList, repeatList, repeatPathway, flankingPathsList])
        

     


    # 7. Format return and move on to the phasing 
    with open(folderName + repeatSpec, 'w') as outfile:
        json.dump(bigDumpList, outfile)

    
def performPhasing(folderName, mummerLink):
    print "performPhasing"
    '''
    1. Interface from alignmentBridge.py : 
        shortToLongMap = formRelatedMap(f2, noisyReads, currentNode, indelRobot, toProcessList)
        cleaner.cleaning([noisyReads,noisyReads] ,shortToLongMap, toProcessList,indelRobot, "init")
        in1List, in2List, out1List, out2List, commonList, longReadToUse  = cleaner.cleaning([noisyReads, noisyReads],shortToLongMap, toProcessList,indelRobot, "vote")
        extendResult = extender.readExtender(in1List, in2List, out1List, out2List, commonList,indelRobot,longReadToUse, True)
    
    2. Format of input data data : 
        bigDumpList.append([flankingList, repeatList, repeatPathway, flankingPathsList])
    
    3. IO : 
        a) Input :
            repeatSpecification.txt, phasingSeedName_Double.fasta, graph G 
        b) Output :
            improved4.fasta
            
    3. Algorithm: 
        a) reformatNoisyReads 
        b) reformatToProcessList
        c) formShortToLongMapping
    
    '''

    json_data = open(folderName + 'repeatSpecification.txt', 'r')
    loadData = json.load(json_data)
    
    G = commonLib.seqGraph(0)
    G.loadFromFile(folderName, "phaseStringGraph1")
    
    lenDicRR = commonLib.obtainLength(folderName, "phasingSeedName_Double.fasta")
    
    lenDicCC = commonLib.obtainLength(folderName, "improved3_Double.fasta")
    N1 = len(lenDicCC)
    
    lenDicCR = dict(lenDicCC.items() + lenDicRR.items())
    
    loadData = filterReverseComp(loadData, N1)
    
    toPhaseList = []
    
    if True:
        for eachitem in loadData:
            # print eachitem
            flankingList, repeatList, repeatPathway, flankingPathsList = eachitem[0], eachitem[1], eachitem[2], eachitem[3] 
            
            noisyReads, dicToOriginal, dicFromOriginal = reformatNoisyReads(folderName, flankingList, repeatList, N1)
            
            toProcessList = reformatToProcessList(folderName , flankingList, repeatList, dicFromOriginal, N1)
    
            shortToLongMap = formShortToLongMapping(folderName, G, toProcessList, dicFromOriginal, dicToOriginal, lenDicCR, N1)
            
            indelRobot = createIndelRobot(folderName)
            
            cleaner.cleaning([noisyReads, noisyReads] , shortToLongMap, toProcessList, indelRobot, "init")
            in1List, in2List, out1List, out2List, commonList, longReadToUse = cleaner.cleaning([noisyReads, noisyReads], shortToLongMap, toProcessList, indelRobot, "vote")
            extendResult = extender.readExtender(in1List, in2List, out1List, out2List, commonList, indelRobot, longReadToUse, True)
            
            if extendResult != -1:
                print "extendResult: ", extendResult
                toPhaseList.append(eachitem + [extendResult])
            
        with open(folderName + 'toPhaseList.txt', 'w') as outfile:
            json.dump(toPhaseList, outfile)

    json_data = open(folderName + 'toPhaseList.txt', 'r')
    toPhaseList = json.load(json_data)
    
    outputResults(folderName, mummerLink, toPhaseList, N1, G)
    
def mainFlow(folderName="SampleTest2/", mummerLink="MUMmer3.23/"):
    print "Hello world"
    
    contigFilename = "improved3"
    readsetFilename = "phasingSeedName"
    optTypeFileHeader = "phaseString"
    contigReadGraph = "phaseStringGraph1"
    repeatFilename = "phaseRepeat.txt"
    repeatSpec = "repeatSpecification.txt"
    
    optionToRun = "xphase"
    
    getAllAssociatedReads(folderName, mummerLink,readsetFilename)
    formReadContigStringGraph(folderName, mummerLink,contigFilename, readsetFilename, optTypeFileHeader , contigReadGraph )
    identifyRepeat(folderName, mummerLink,contigFilename,contigReadGraph, repeatFilename, optionToRun )
    defineRepeatAndFlanking(folderName, mummerLink,contigFilename,contigReadGraph,repeatFilename,repeatSpec )
    performPhasing(folderName, mummerLink)




# mainFlow()  

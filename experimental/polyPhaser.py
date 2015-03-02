import json
import os
from itertools import groupby

import associatedReadFinder
import readContigGraphFormer
import repeatFinder
import repeatFlankingDefiner

from finisherSCCoreLib import IORobot
from finisherSCCoreLib import graphLib


from xPhaserLib import cleaner
from xPhaserLib import extender
from xPhaserLib import common


def mainFlow(folderName="SampleTest2/", mummerLink="MUMmer3.23/"):
    print "Hello world"
    
    contigFilename = "improved3"
    readsetFilename = "phasingSeedName"
    optTypeFileHeader = "phaseString"
    contigReadGraph = "phaseStringGraph1"
    repeatFilename = "phaseRepeat.txt"
    repeatSpec = "repeatSpecification.txt"
    
    optionToRun = "xphase"
    
    if True:
        associatedReadFinder.getAllAssociatedReads(folderName, mummerLink,readsetFilename)
    
        readContigGraphFormer.formReadContigStringGraph(folderName, mummerLink,contigFilename, readsetFilename, optTypeFileHeader , contigReadGraph )
    
        repeatFinder.identifyRepeat(folderName, mummerLink,contigFilename,contigReadGraph, repeatFilename, optionToRun )
    
    repeatFlankingDefiner.defineRepeatAndFlanking(folderName, mummerLink,contigFilename,contigReadGraph,repeatFilename,repeatSpec )    
    
    performPhasing(folderName, mummerLink)



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
    
    G = graphLib.seqGraph(0)
    G.loadFromFile(folderName, "phaseStringGraph1")
    
    lenDicRR = IORobot.obtainLength(folderName, "phasingSeedName_Double.fasta")
    
    lenDicCC = IORobot.obtainLength(folderName, "improved3_Double.fasta")
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
    G2 = graphLib.seqGraph(N1)
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
                    
                    newNode = graphLib.seqGraphNode(myindex)
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
    
    IORobot.readContigOut(folderName, mummerLink, graphFileName, combinedName, "improved4.fasta", "outOpenListphaing", nameDic)


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



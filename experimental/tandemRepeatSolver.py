

import json
from itertools import groupby
from operator import itemgetter

from finisherSCCoreLib import IORobot
from finisherSCCoreLib import graphLib
from finisherSCCoreLib import alignerRobot
import associatedReadFinder
import readContigGraphFormer
import repeatFinder

import repeatFlankingDefiner
import abunHouseKeeper


'''
Interface : 
    Input : contigs.fasta, raw_reads.fasta
    Output : tandemSolved.fasta
    Algorithm : 
        1) Color nodes

        2) Identify skeleton of tandem repeats and associated redundant contigs

        For each repeat group formed
        3) If 1 in and 1 out for that repeat, 
            a) Form a repeat template
            b) Align raw reads
            c) Do count 
            d) Make join and remove redundant contigs

'''

    
def colorNodes(folderName, mummerPath,sourceFilename, contigFilename, readsetFilename):
    print "colorNodes"
    lenDic = IORobot.obtainLength(folderName, sourceFilename+".fasta")
    print lenDic
    thresForShort = 15000
    shortList = []
    longList = []
    for eachitem in lenDic:
        if lenDic[eachitem] > thresForShort:
            longList.append(eachitem)
        else:
            shortList.append(eachitem)
    
    IORobot.putListToFileO(folderName, sourceFilename+".fasta", contigFilename, longList)
    IORobot.putListToFileO(folderName, sourceFilename+".fasta", readsetFilename, shortList)
    
    
def skeletonIdentification(folderName, mummerLink, contigFilename, readsetFilename, optTypeFileHeader, contigReadGraph , repeatFilename, repeatSpec, optionToRun):
    print "skeletonIdentification"
    associatedReadFinder.getAllAssociatedReads(folderName, mummerLink,readsetFilename)
    readContigGraphFormer.formReadContigStringGraph(folderName, mummerLink, contigFilename, readsetFilename, optTypeFileHeader , contigReadGraph)
    repeatFinder.identifyRepeat(folderName, mummerLink, contigFilename, contigReadGraph, repeatFilename, optionToRun)
    repeatFlankingDefiner.defineRepeatAndFlanking(folderName, mummerLink, contigFilename, contigReadGraph, repeatFilename, repeatSpec)

def findCoverageFromRawData(folderName):
    contigLenDic  = IORobot.obtainLength(folderName, "contigs.fasta")
    readLenDic = IORobot.obtainLength(folderName, "raw_reads.fasta")

    G = 0 
    NL = 0 
    for eachitem in contigLenDic:
        G = G+ contigLenDic[eachitem]
        
    for eachitem in readLenDic:
        NL = NL+ readLenDic[eachitem]
    
    c = (NL*1.0)/G 
    print c 
    return c 
def DFSwithPath(G, x , pathList, N1, isTerminate):
    if not isTerminate:
        if x.visited == 0:
            x.visited = 1
            returnPathList = pathList
            for eachnext in x.listOfNextNodes:
                i = eachnext[0]
                if i >= N1:
                    isTerminate, returnPathList = DFSwithPath(G, G.graphNodesList[i], pathList + [i], N1,isTerminate)
                    if isTerminate:
                        return isTerminate, returnPathList
            x.visited = -1
            return isTerminate, returnPathList
        elif x.visited == 1:
            return True, pathList
        elif x.visited == -1:
            return False, pathList

        
    
def resolvingTandem(folderName, mummerPath, contigReadGraph,contigFilename, readsetFilename, optTypeFileHeader, repeatSpec):
    print "resolvingTandem"
    '''
    Input : repeat info 
    Output : count, join. 
    
    Algorithm: 
    1. Find loops
    2. Form repeat
    3. Form chain of repeat copies back to back
    4. Align reads
    5. Calculate extra bases beyond flanking region
    6. Calculate count
    7. Join the contigs
    '''
    # 0 ) Load all the data
    thres = 5 
    
    G = graphLib.seqGraph(0)
    G.loadFromFile(folderName, contigReadGraph)
    lenDicCC = IORobot.obtainLength(folderName, contigFilename+"_Double.fasta")
    N1 = len(lenDicCC)

    maxDuplicate = 10
    repeatTempFilename = "tandemRepeatTemplate.fasta"
    mummerFile = "myTandemRepeatTemplate"
    


    myContigsDic = IORobot.loadContigsFromFile(folderName, readsetFilename+"_Double.fasta")    
    lenDicRR = IORobot.obtainLength(folderName, readsetFilename + "_Double.fasta")
    
    header = optTypeFileHeader + "RR"
    dataListRR = alignerRobot.extractMumData(folderName, header + "Out")
    dataListRR = abunHouseKeeper.filterData(dataListRR, lenDicRR)
    dataListRRDic = {}
    for eachitem in dataListRR: 
        if eachitem[2] < thres:
            dataListRRDic[eachitem[-2] +";"+eachitem[-1]] = eachitem[4]

    header = optTypeFileHeader + "CR"
    lenDicCC = IORobot.obtainLength(folderName, contigFilename + "_Double.fasta")
    lenDicCR = dict(lenDicCC.items() + lenDicRR.items())
    
    dataListCR = alignerRobot.extractMumData(folderName, header + "Out")
    dataListCR = abunHouseKeeper.filterData(dataListCR, lenDicCR)
    dataListCRDic = {}
    for eachitem in dataListCR: 
        if eachitem[2] < thres:
            dataListCRDic[eachitem[-2] +";"+eachitem[-1]] = eachitem[4]

    print dataListCRDic



    json_data = open(folderName + repeatSpec, 'r')
    loadData = json.load(json_data)
    
    contigsTmp = IORobot.loadContigsFromFile(folderName, contigFilename+"_Double.fasta")
    readTmp = IORobot.loadContigsFromFile(folderName, readsetFilename + "_Double.fasta")

    happyTandemList = {}
    
    
    
    for eachrepProfile in loadData:
        # 1) 
        startContig = eachrepProfile[-1][0][0]
        isTerminate, returnPathList = DFSwithPath(G, G.graphNodesList[startContig], [startContig], N1, False)
       
        # 2) 
        if isTerminate:
            v = returnPathList[-1]
            i =0 
            tandemPath = []
            while i < len(returnPathList):
                if returnPathList[i] == v:
                    tandemPath = returnPathList[i:]
                    i = len(returnPathList)
                i = i +1
                
            print returnPathList
            print tandemPath
        # 3) [fix it when have time later ; to just use graph; bug at the min thing]
        
        repeatContent = ""
    
        for kk in range(len(tandemPath[0:-1])): 
            eachitem = tandemPath[kk]- N1
            nextitem = tandemPath[kk+1] - N1
            readName = "Read" + str(eachitem/2) + "_"
            nextReadName = "Read" + str(nextitem/2) + "_"
            if eachitem %2 ==0 :
                readName = readName + "p"
            elif eachitem %2 ==1:
                readName = readName + "d"
            
            if nextitem %2 ==0 :
                nextReadName = nextReadName + "p"
            elif nextitem %2 ==1:
                nextReadName = nextReadName + "d"
            
            overlap = dataListRRDic[readName + ";" + nextReadName]
            print overlap
            repeatContent = repeatContent +  myContigsDic[readName][0:-overlap]
            
        print "len(repeatContent)", len(repeatContent)
        
        fout = open(folderName + repeatTempFilename, 'w')
        fout.write(">RepeatSegment\n")
        repeatContentLarge = ""
        
        for i in range(maxDuplicate):
            fout.write(repeatContent)
            repeatContentLarge= repeatContentLarge + repeatContent
        fout.close()
        
        # 4)
        repeatReadList =  eachrepProfile[1]
        
        myList= []
        for eachitem in repeatReadList:
            
            readName = "Read" + str((eachitem- N1)/2) + "_"
    
            if eachitem %2 ==0 :
                readName = readName + "p"
            elif eachitem %2 ==1:
                readName = readName + "d"
            myList.append(readName)
            
        IORobot.putListToFileO(folderName, readsetFilename+"_Double.fasta", "toAlignReads", myList)
        
        if True:
            alignerRobot.useMummerAlign(mummerPath, folderName,mummerFile , repeatTempFilename, "toAlignReads.fasta")
        
        dataList = alignerRobot.extractMumData(folderName, mummerFile+"Out")
        
        
        # 5)
        totalBasesMatch = 0
        lrepeat = len(repeatContent)
        c = findCoverageFromRawData(folderName)
        
        
        # print "dataList[0]", dataList[0]
        dataList.sort(key = itemgetter(-1))
        for key, values in  groupby(dataList,itemgetter(-1)):
            maxValue = -1
            for eachsub in values:
                if eachsub[5] > maxValue:
                    maxValue = eachsub[5]
    
            #print key, maxValue
            totalBasesMatch = totalBasesMatch + maxValue
        
    
        print c, lrepeat, totalBasesMatch
        ct = totalBasesMatch*1.0/(c*lrepeat)
        print "BIG NUMBER of THE DAY: ", ct
    
        # 6) 
        # a) find the starting point 
        startContig = eachrepProfile[-1][0][0]
        firstRead = eachrepProfile[-1][0][1]-N1

        contigName = "Contig"+ str(startContig/2)
        if startContig %2 == 0:
            contigName = contigName + "_p"
        elif startContig%2 ==1:
            contigName = contigName + "_d"
        
        readName = "Read"+ str(firstRead/2)
        if firstRead %2 == 0:
            readName = readName + "_p"
        elif firstRead%2 ==1:
            readName = readName + "_d"
        
        overlapFirst = dataListCRDic[contigName+";"+readName]
        tmpCombine = contigsTmp[contigName][0:-overlapFirst] + readTmp[readName]
        
        f1 = open(folderName + "firstOverlap.fasta", 'w')
        f1.write(">combined\n")
        f1.write(tmpCombine)
        f1.close()
        
        if True:
            alignerRobot.useMummerAlign(mummerPath, folderName,"myFirstOverlap" , repeatTempFilename, "firstOverlap.fasta")
        
        dataList = alignerRobot.extractMumData(folderName, "myFirstOverlap"+"Out")
        
        dataList.sort(key = itemgetter(0))
        maxVal = -1
        maxItm = []
        for eachi in dataList:
            if eachi[5] > maxVal:
                maxVal = eachi[5]
                maxItm = eachi
        
        print maxItm
        if len(maxItm) > 0 :
            repeatStart = maxItm[0]
            contigEnd = maxItm[2]
        else:
            repeatStart = 0
            contigEnd = -1
        # b) format return : prepare the repeat template 
        print "ct*lrepeat", int(repeatStart + ct*lrepeat)
        print "repeatStart", repeatStart
        happyTandemList[contigName]= repeatContentLarge[repeatStart:int(repeatStart + ct*lrepeat)]
        contigsTmp[contigName] = tmpCombine[0:contigEnd]
        print "len(contigsTmp[contigName])", len(contigsTmp[contigName])
        print "len(happyTandemList[contigName])", len(happyTandemList[contigName])
        
    # 7) Combine all the repeat information and do the join
    
    leaderList = [i for i in range(len(contigsTmp))]
    for eachrepProfile in loadData:
        startContig = eachrepProfile[-1][0][0]
        endContig = eachrepProfile[-1][-1][-1]
        leaderContig = leaderList[startContig]
        
        leaderName = parseIDToName(leaderContig)
        endName = parseIDToName(endContig)
        startName = parseIDToName(startContig)
        
        contigsTmp[leaderName] = contigsTmp[leaderName] + happyTandemList[startName]
        
        if endContig != leaderContig:
            contigsTmp[leaderName] = contigsTmp[leaderName] + contigsTmp[endName]
            contigsTmp[endName] = ""
            leaderList[endContig] = leaderContig
        
    
    leaderAgg = [[] for i in range(len(leaderList))]
    for i in range(len(leaderList)):
        leaderAgg[leaderList[i]].append(i) 
    
    checkingList = [False for i in range(N1)]
    
    fout = open(folderName + "tademResolved.fasta", 'w')
    
    counter = 0
    for eachcontig in contigsTmp:
        id = abunHouseKeeper.parseEdgeNameToID(eachcontig, 'C')
        if checkingList[id/2] == False:
        
            fout.write(">Segkk"+str(counter)+ "\n")
            
            fout.write(contigsTmp[eachcontig])
            counter = counter + 1    
            for eachkk in leaderAgg[leaderList[id]]:
                checkingList[eachkk/2] = True
    
    fout.close()
    
def parseIDToName(k):
    name = "Contig"+ str(k/2)
    if k%2 ==0 :
        name = name + "_p"
    else:
        name = name + "_d"
    return name
    
def mainFlowForTandemResolve(folderName, mummerPath):
    
    sourceFilename = "contigs"
    contigFilename = "tandemLongContigs"
    shortContigFilename = "tandemShortContigs"
    readsetFilename = "tandemReads"
    
    optTypeFileHeader = "tandemString"
    contigReadGraph = "tandemStringGraph1"
    repeatFilename = "tandemRepeat.txt"
    repeatSpec = "tandemRepeatSpecification.txt"
    
    optionToRun = "tandem"

    if True:
        colorNodes(folderName, mummerPath, sourceFilename, contigFilename, shortContigFilename)
        skeletonIdentification(folderName, mummerPath, contigFilename, readsetFilename, optTypeFileHeader, contigReadGraph , repeatFilename, repeatSpec,optionToRun)
    
    resolvingTandem(folderName, mummerPath,contigReadGraph,contigFilename,readsetFilename, optTypeFileHeader, repeatSpec)
    


    

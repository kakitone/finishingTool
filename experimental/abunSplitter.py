import json
import numpy as np
import os
from itertools import groupby
from operator import itemgetter
import copy
import sys

from finisherSCCoreLib import alignerRobot
from finisherSCCoreLib import graphLib
from finisherSCCoreLib import IORobot
from finisherSCCoreLib import houseKeeper

import associatedReadFinder
import readContigGraphFormer
import repeatFinder
import abunHouseKeeper


### Abundance split and its subroutines
def obtainNonEmpty(repeatPairs):
    newRepeatPairs = []
    for eachitem in repeatPairs:
        if len(eachitem[0]) > 0 and len(eachitem[1]) > 0:
            newRepeatPairs.append(eachitem)
    
    return newRepeatPairs

def convert4to1base(i):
    return i/4

def convert4to2base(i):
    return i/2

def getCt(inList ,myCountDic):
    newInList = []
    for i in inList: 
        tmp1 = convert4to1base(i)
        tmp2 = convert4to2base(i)
        newInList.append([tmp2 , myCountDic["Segkk"+str(tmp1)]])
    return newInList 

def satisfyMatch(initem ,newOutList, sd):
    found = -1
    inIndex , inCt = initem[0], initem[1]
    # First check pt
    targetItem = -1
    for outitem in newOutList:
        outIndex, outCt = outitem[0] , outitem[1]
        if abs(outCt - inCt) < 0.5*sd and inIndex != outIndex:
            targetItem = outIndex
            
    # Second check pt
    rejection = False
    for outitem in newOutList:
        outIndex, outCt = outitem[0] , outitem[1]
        if outIndex != targetItem :
            if abs(outCt - inCt) <= 2.01*sd :
                rejection = True
                
    # Combined check 
    if not rejection and targetItem != -1:
        found = targetItem 
    else:
        found = -1
        
    return found
     
def determindMatch(inList, outList, myCountDic):
    resolvedList = []
    newInList , newOutList = [], []
    
    newInList = getCt(inList ,myCountDic)
    newOutList = getCt(outList, myCountDic)
    
    sizeList = []
    for eachitem in myCountDic:
        sizeList.append(myCountDic[eachitem])
        
    sd = np.std(sizeList)
    
    for eachitem in newInList:
        found = satisfyMatch(eachitem ,newOutList, sd)
        if found != -1 :
            resolvedList.append([eachitem[0], found])
        
            
    return resolvedList

def addEdges(G, resolvedList):
    for eachitem in resolvedList:
        G.insertEdge(eachitem[0], eachitem[1], 1997)
    


def abunSplit(folderName, mummerLink, myCountDic):
    
    '''
    Input : repeatSpecification.txt , myCountDic.json, improved3.fasta, raw_reads.fasta
    Output : abunsplit.fasta
    
    Algorithm : 
    
    1. Load data from various sources [various json files]
    
    2. For each repeat interior:
        a) identify the abundances associated with in/out contigs
        b) perform a split and record the split
    
    3. Use split results to generate contigs [may already exist in newPhasing.py ] 
        a) use a graph to capture the split results 
        b) use reads to fill in any gaps 
        c) read out the contigs 
    
    '''
    json_data = open(folderName + "phaseRepeat.txt", 'r')
    repeatPairs = json.load(json_data)
    repeatPairs = obtainNonEmpty(repeatPairs)
    
    N1 = len(myCountDic)*2
    
    G = graphLib.seqGraph(N1)
    
    for eachitem in repeatPairs:
        inList, outList = eachitem[0], eachitem[1]
        resolvedList = determindMatch(inList, outList, myCountDic)
        addEdges(G, resolvedList)
        
    G.condense()
    IORobot.extractGraphToContigs(G, folderName, mummerLink, "abun.fasta", "improved3_Double.fasta")
    
    

def splitter(folderName, mummerLink):

    '''
    Input : repeatSpecification.txt , myCountDic.json, improved3.fasta, raw_reads.fasta
    Output : abunsplit.fasta
    
    Algorithm : 
    
    1. Load data from various sources [various json files]
    
    2. For each repeat interior:
        a) identify the abundances associated with in/out contigs
        b) perform a split and record the split
    
    3. Use split results to generate contigs [may already exist in newPhasing.py ] 
        a) use a graph to capture the split results 
        b) use reads to fill in any gaps 
        c) read out the contigs 
    
    '''
    
    
    with open(folderName + 'myCountDic.json') as f:
        myCountDic = json.load(f)
        
    abunSplit(folderName, mummerLink, myCountDic)    
    


'''
Input : string_graph_3, improved3.fasta, raw_reads.fasta
Output : string_graph_4 with weights [need a data structure to store the weight on node]

Algorithm : 
1. Find your favorite mappers to map read back
    a. MUMmer, Bowtie, bbmap, any that works V 
    b. And then write a short parser to parse the results V 
2. Calculate count on the abundances 
    a. Aggregate by taking average [put weights on bin along contigs]
    b. Inheritance and a subclass 
3. Find your favorite graphical tool to display 
    a. Use a javascript library [halfviz should just work ! put weight on edge ]

'''

def evaluateCoverage(dataList, lenDic, readLenDic, folderName,mummerLink, continueFilter):
    
    myCountDic = {}
    for eachitem in lenDic:
        myCountDic[eachitem] = 0
            
    dataList.sort(key = itemgetter(-1)) 
    
    ctkk, ctbase = 0, 0
    toAddBackDic = copy.deepcopy(readLenDic)
    
    for key, items in groupby(dataList, itemgetter(-1)):
        maxMatch = -1
        bestname = ""
        
        for eachitem in items:
            ct = eachitem[6]/100.0 * eachitem[4]
            if ct > maxMatch:
                maxMatch = ct 
                bestname = eachitem[-2]
        myCountDic[bestname] += readLenDic[key] 
        
        ctkk = ctkk + 1 
        ctbase = ctbase + readLenDic[key]
        toAddBackDic[key] = -1
    
    cttot = 0
    for eachitem in readLenDic:
        cttot = cttot + readLenDic[eachitem]
        
    print "Missed coverage  ", (cttot - ctbase)/(4.7*pow(10, 6))
    print "percentage miss read", (len(readLenDic) - ctkk)/(1.0*len(readLenDic)) 
    
    toAddReadList = []
    for eachitem in toAddBackDic:
        if toAddBackDic[eachitem] >= 0 :
            toAddReadList.append(eachitem)
    
    '''
    This part need the most parallelism because it is most intense with -l 10 
    split V, workerList V , combine 
    '''
    
    if continueFilter:
        numberOfFiles= 20
        
        IORobot.putListToFileO(folderName, "raw_reads.fasta" , "selected_raw", toAddReadList)
        
        bindir =  os.path.abspath(os.path.dirname(sys.argv[0]))   
        command = bindir + "/finisherSCCoreLib/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + folderName + "selected_raw.fasta"
        os.system(command)
        
        workerList = []
        
        for dummyI in range(1, numberOfFiles + 1):
            indexOfMum = ""
            if dummyI < 10:
                indexOfMum = "0" + str(dummyI)
            else:
                indexOfMum = str(dummyI)
           
            outputName, referenceName, queryName, specialName= "outAbunRefine"+indexOfMum, "improved3.fasta", "selected_raw.part-"+ indexOfMum + ".fasta",  "abunMissOut" + indexOfMum
            workerList.append([outputName, referenceName, queryName, specialName])
            
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, workerList, houseKeeper.globalParallel ,specialForRaw = True, refinedVersion = True)
        alignerRobot.combineMultipleCoorMum( True, mummerLink, folderName, "outAbunRefine", "abunMissOut", numberOfFiles)
        

        
    for i in range(len(myCountDic)):
        eachitem = "Segkk"+str(i)
        print eachitem , myCountDic[eachitem]/(1.0*lenDic[eachitem])
        myCountDic[eachitem] = myCountDic[eachitem]/(1.0*lenDic[eachitem])
        
    return myCountDic
    

def generateAbundanceGraph(folderName, mummerLink):
    
    
    print "generateAbundanceGraph"
    
    '''
    1. Find your favorite mappers to map read back
        a. MUMmer, Bowtie, bbmap, any that works V 
        b. And then write a short parser to parse the results V 
    '''
    numberOfFiles = 20
    workerList = []
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
        
        '''
        "outGapFillRefine"+indexOfMum , "smaller_improvedContig.fasta",  "relatedReads_Double.part-" + indexOfMum + ".fasta",  "fromMumRefine" + indexOfMum
        '''
        outputName, referenceName, queryName, specialName= "outAbun"+indexOfMum, "improved3.fasta", "raw_reads.part-"+ indexOfMum + ".fasta",  "outAbun" + indexOfMum
        workerList.append([outputName, referenceName, queryName, specialName])
    
    if True:
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, workerList, houseKeeper.globalParallel ,False)
        '''
        command = mummerLink + "nucmer --maxmatch --nosimplify -p " + folderName + "out " + folderName + "improved3.fasta "+folderName+"raw_reads.part-" + indexOfMum + ".fasta"
        os.system(command)
    
        command = mummerLink + "show-coords -r " + folderName + "out.delta > " + folderName + "fromMumAbun" + indexOfMum
        os.system(command)
        '''
        
    dataList = []
    
    for i in range(1, 1+numberOfFiles): 
        if i < 10:
            indexOfMum = "0" + str(i)
        else:
            indexOfMum = str(i)
        dataList = dataList+ alignerRobot.extractMumData(folderName, "outAbun"+ str(indexOfMum)+"Out")
    

    '''
    2. Calculate count on the abundances 
        a. Aggregate by taking average [put weights on bin along contigs]
        b. Inheritance and a subclass 
    '''
         
    lenDic = IORobot.obtainLength(folderName, "improved3.fasta")
    readLenDic = IORobot.obtainLength(folderName , "raw_reads.fasta")
    

    myCountDic = {}
    for eachitem in lenDic:
        myCountDic[eachitem] = [0 for i in range(lenDic[eachitem])]

    thres = 30
    lenSum = 0
    extraDataList= []
    
    
    print "len(dataList)", len(dataList)
    
    if not abunHouseKeeper.abunGlobalAvoidrefine: 
        myCountDic =  evaluateCoverage(dataList, lenDic, readLenDic, folderName, mummerLink,  True)
        extraDataList = alignerRobot.extractMumData(folderName, "abunMissOut" )
    else:
        extraDataList = []
        
    dataList = dataList + extraDataList
    myCountDic = evaluateCoverage(dataList, lenDic, readLenDic, folderName, mummerLink,False)
    
    with open(folderName + 'myCountDic.json', 'w') as f:
        json.dump(myCountDic, f)

    
    return myCountDic

def mainFlow(folderName, mummerLink):
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
                
    if True:
        myCountDic = generateAbundanceGraph(folderName, mummerLink)
        
    if True :
        splitter(folderName, mummerLink)
        
    os.system("cp selected_raw.part-* "+ folderName )
    os.system("rm selected_raw.part-*")
        
        

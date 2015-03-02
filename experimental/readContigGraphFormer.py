import os
import sys 

import abunHouseKeeper

from finisherSCCoreLib import IORobot
from finisherSCCoreLib import alignerRobot
from finisherSCCoreLib import graphLib
from finisherSCCoreLib import houseKeeper


def addDataToList(dataList, G, startIndex1, startIndex2, type1, type2):

    threshold = 20
    
    for eachitem in dataList:
        wt = min(eachitem[4] , eachitem[5])
        
        if eachitem[0] < threshold:

            j = abunHouseKeeper.parseEdgeNameToID(eachitem[-2], type1) + startIndex1
            i = abunHouseKeeper.parseEdgeNameToID(eachitem[-1], type2) + startIndex2
        else:
            j = abunHouseKeeper.parseEdgeNameToID(eachitem[-1], type2) + startIndex2
            i = abunHouseKeeper.parseEdgeNameToID(eachitem[-2], type1) + startIndex1
        
        
        G.insertEdge(i, j, wt) 

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
   
def alignerSubRoutine(folderName ,referenceFile,  queryFile, mummerLink, header ):   
    #alignerRobot.useMummerAlign(mummerLink, folderName, header, referenceFile, queryFile)
    numberOfFiles = 20
    bindir =  os.path.abspath(os.path.dirname(sys.argv[0]))   
    command = bindir + "/finisherSCCoreLib/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + folderName + queryFile
    os.system(command)
    
    os.system("cp *.fasta " + folderName )
    os.system("rm *.fasta ")
    
    workerList = []
    
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
       
        outputName, referenceName, queryName, specialName= header+indexOfMum, referenceFile,queryFile[0:-6]+".part-"+ indexOfMum + ".fasta" ,  header + indexOfMum
        workerList.append([outputName, referenceName, queryName, specialName])
        
    alignerRobot.useMummerAlignBatch(mummerLink, folderName, workerList, houseKeeper.globalParallel ,specialForRaw = False, refinedVersion = False)
    alignerRobot.combineMultipleCoorMum( True, mummerLink, folderName, header,header +"Out", numberOfFiles)
                    

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

    IORobot.writeToFile_Double1(folderName, contigFilename + ".fasta", contigFilename + "_Double.fasta", "contig")
    IORobot.writeToFile_Double1(folderName, readsetFilename + ".fasta", readsetFilename + "_Double.fasta", "reads")
    
    
    header, referenceFile, queryFile = optTypeFileHeader + "CC", contigFilename + "_Double.fasta" , contigFilename + "_Double.fasta"
    if True:
        alignerRobot.useMummerAlign(mummerLink, folderName, header, referenceFile, queryFile)

    lenDicCC = IORobot.obtainLength(folderName, contigFilename + "_Double.fasta")
    dataListCC = alignerRobot.extractMumData(folderName, header + "Out")
    dataListCC = abunHouseKeeper.filterData(dataListCC, lenDicCC)
    
    header, referenceFile, queryFile = optTypeFileHeader + "RR", readsetFilename + "_Double.fasta" , readsetFilename + "_Double.fasta"
    
    
    lenDicRR = IORobot.obtainLength(folderName, readsetFilename + "_Double.fasta")
    
    if not abunHouseKeeper.abunGlobalRRDisable:
        if True:
            alignerSubRoutine(folderName ,referenceFile,  queryFile, mummerLink, header )
    
        dataListRR = alignerRobot.extractMumData(folderName, header + "Out")
        dataListRR = abunHouseKeeper.filterData(dataListRR, lenDicRR)
        for eachitem in dataListRR:
            if eachitem[-2] == "Read164_p" and eachitem[-1] == "Read159_p" :    
                print "debug" , eachitem
            if eachitem[-1] == "Read164_p" and eachitem[-2] == "Read159_p" :    
                print "debug" , eachitem
            
        dataListRR = abunHouseKeeper.filterDataIdentical(dataListRR, lenDicRR)

    else:
        dataListRR = []
    
    header, referenceFile, queryFile = optTypeFileHeader + "CR", contigFilename + "_Double.fasta" , readsetFilename + "_Double.fasta"
    if True:
        alignerSubRoutine(folderName ,referenceFile,  queryFile, mummerLink, header )
        #alignerRobot.useMummerAlign(mummerLink, folderName, header, referenceFile, queryFile)
    
    lenDicCR = dict(lenDicCC.items() + lenDicRR.items())
    dataListCR = alignerRobot.extractMumData(folderName, header + "Out")
    dataListCR = abunHouseKeeper.filterData(dataListCR, lenDicCR)
            
    numberOfNodes = len(lenDicCR) 
    G = graphLib.seqGraph(numberOfNodes)
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
    
    
    
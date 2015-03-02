import json
from itertools import groupby
import abunHouseKeeper
import abunGraphLib

from finisherSCCoreLib import graphLib
from finisherSCCoreLib import IORobot                


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
    G = graphLib.seqGraph(0)
    G.loadFromFile(folderName, contigReadGraph)
    # G.reportEdge()
    lenDicCC = IORobot.obtainLength(folderName, contigFilename+"_Double.fasta")
    
    adjacencyList = [[] for i in range(len(lenDicCC))]
    
    N1 = len(lenDicCC)
    
    
    # # Debug
    # for i in range(14):
    #    debugGraphPath(i, 2, G, N1)
    # # End Debug
    
    for i in range(len(lenDicCC)):
        adjacencyList[i] = abunGraphLib.findAllReachable(i, N1, G) 
        print "i, adjacencyList[i] : ", i , adjacencyList[i]
    
    # ## (b) formation of bipartite graph
    if optionToRun == "tandem" :
        newAdjacencyList = adjacencyList
    elif optionToRun == "xphase": 
        newAdjacencyList = abunGraphLib.filterEdge(adjacencyList, folderName, contigFilename)
    
    G2 = abunGraphLib.seqGraphWt(N1 * 2)
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
                
        
        repeatList.append([abunHouseKeeper.getDistinct(leftList), abunHouseKeeper.getDistinct(rightList)])
           
    with open(folderName + repeatFilename, 'w') as outfile:
        json.dump(repeatList, outfile)

    
    json_data = open(folderName + repeatFilename, 'r')
    loadData = json.load(json_data)
    
    
    assert(loadData == repeatList)
    
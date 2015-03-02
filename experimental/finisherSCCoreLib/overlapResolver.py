import IORobot
from operator import itemgetter
from itertools import groupby
import graphLib
import houseKeeper


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



# ## 1) Fetch best successor and best predecessor for each contigs (I: noEmbed.fasta ;  O:  left_connect , right_connect )
def fetchSuccessor(folderName , mummerLink): 
    
    print "fetchSuccessor"
    left_connect, right_connect = [], [] 
        
    print "Direct greedy"
    numberOfContig, dataSet = IORobot.obtainLinkInfo(folderName, mummerLink, "noEmbed", "greedy")
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
    
    G = graphLib.seqGraph(numberOfNodes)
        
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
    
    G2 = graphLib.seqGraph(0)
    G2.loadFromFile(folderName, "condensedGraph.txt")
    
    houseKeeper.compareGraphUnitTest(G, G2)
    G.reportDummyUsefulNode()
    G.reportEdge()
    
    graphFileName = "condensedGraph.txt"
    contigFile = "noEmbed_Double.fasta"
    outContigFile = "improved.fasta"
    outOpenList = "openZone.txt"
    
    IORobot.readContigOut(folderName, mummerLink, graphFileName, contigFile, outContigFile, outOpenList)
   

# ## 3) X-phased seqGraph (I: startList, graphNodes; O: startList, graphNodes )


# ## 4) EC reduction (I: startList, graphNodes ; O: startList, graphNodes )
def ECReduction(folderName , mummerLink):
    print "ECReduction" 
    

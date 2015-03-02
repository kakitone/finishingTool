import json
import abunGraphLib

from finisherSCCoreLib import IORobot
from finisherSCCoreLib import graphLib


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
    G = abunGraphLib.seqGraphWt(0)
    G.loadFromFile(folderName, contigReadGraph)
    Grev = abunGraphLib.formReverseGraph(G)
    
    json_data = open(folderName + repeatFilename, 'r')
    repeatList = json.load(json_data)
    
    lenDicCC = IORobot.obtainLength(folderName, contigFilename+"_Double.fasta")
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
                
            
            abunGraphLib.markReachableIndices(G, Grev, kkIn, kkOut, N1)
            
            # 2. Marks inside nodes
            singleMissList, allPassList = abunGraphLib.markInsideNodes(G, kkIn, kkOut)
            for i in range(4): 
                print "len(singleMissList[i]), len(allPassList)", len(singleMissList[i]), len(allPassList)

            # 3. Finds start/end of repeat
            myStartIndex, myEndIndex = abunGraphLib.markStartEndNodes(G, rIn, rOut, singleMissList, allPassList)
            print myStartIndex, myEndIndex
            
            # 4. Find repeat interior by shortest path joining S/E
            repeatPathway = abunGraphLib.markInterior(G , myStartIndex, myEndIndex, N1)
            print "repeatPathway", repeatPathway
            #checkPathLength(repeatPathway, G, N1, folderName)
            
            # 5. Find flanking region by shortest path search again
            flankingPathsList = abunGraphLib.markFlankingRegion(G, rIn, rOut, myStartIndex, myEndIndex, N1)
            print flankingPathsList
            
            # 6. Find associated reads by graph node query
            flankingList, repeatList = abunGraphLib.markAssociatedReads(G, singleMissList, allPassList)
            
            # ## Experimental
            repeatList = allPassList
            
            # ## End Experimental
            for eachlist in flankingList:
                print len(eachlist), len(repeatList)
            
            bigDumpList.append([flankingList, repeatList, repeatPathway, flankingPathsList])
        

     


    # 7. Format return and move on to the phasing 
    with open(folderName + repeatSpec, 'w') as outfile:
        json.dump(bigDumpList, outfile)

    
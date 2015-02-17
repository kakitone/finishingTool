import IORobot
import graphLib


def loadEdgeFromBlockedReads(folderName):
    extraEdges = []
    f = open(folderName + 'extraConnect.txt', 'r')
    tmp = f.readline().rstrip()
    while len(tmp) > 0 :
        myInfo = tmp.split(';')
        inNode = myInfo[0].split('_')
        outNode = myInfo[1].split('_')
        inWt = int(myInfo[2])
        outWt = int(myInfo[3])
        
        myInIndex , myOutIndex = 1997, 1997 
        if inNode[1] == 'p':
            myInIndex = 2 * int(inNode[0])
        else:
            myInIndex = 2 * int(inNode[0]) + 1 
        
        if outNode[1] == 'p':
            myOutIndex = 2 * int(outNode[0])
        else:
            myOutIndex = 2 * int(outNode[0]) + 1 
        
        extraEdges.append([myInIndex, myOutIndex, inWt, outWt])
        tmp = f.readline().rstrip()
    
    f.close()
    print extraEdges
    return extraEdges




def xPhased(folderName , mummerLink):
    # ## Repeat resolution  [Proxy for MB]
    # 1. Re-form the contig string graph with ALL connections from contigs only V
    # 2. Log down the reads and associated blocked contigs V 
    # 3. Use reads to connect;
    # 4. Transform graph by identifying 1 successor/predecessor case ; Condense(important);
    # 5. Read out contigs
    
    numberOfContig, dataSet = IORobot.obtainLinkInfo(folderName, mummerLink, "improved2", "mb")
    
    lenDic = IORobot.obtainLength(folderName, "improved2_Double.fasta")
    
    confidenLenThres = 0 
    
    G = graphLib.seqGraph(numberOfContig)
    extraEdges = loadEdgeFromBlockedReads(folderName)
    
    for eachitem in dataSet:
        # print eachitem
        wt, myin, myout = eachitem
        myInData = myin[6:].split('_')
        myOutData = myout[6:].split('_')
        
        if myInData[1] == 'p':
            offsetin = 0
        else:
            offsetin = 1
        
        if myOutData[1] == 'p':
            offsetout = 0
        else:
            offsetout = 1
            
        i = int(myInData[0]) * 2 + offsetin
        j = int(myOutData[0]) * 2 + offsetout
        
        ck = False
        
        for eachedge in extraEdges:
            mystart, myend, len1, len2 = eachedge[0], eachedge[1], eachedge[2] , eachedge[3]
            if [i, j] == [mystart, myend] and min(len1, len2) >= wt and lenDic[myin] >= confidenLenThres and lenDic[myout] >= confidenLenThres:
                ck = True
                
        if ck:
            G.insertEdge(i, j, wt)
    
    
    # G.reportEdge()
    G.MBResolve()
    G.reportEdge()
        
    G.saveToFile(folderName, "condensedGraphMB.txt")
    graphFileName = "condensedGraphMB.txt"
    contigFile = "improved2_Double.fasta"
    outContigFile = "improved3.fasta"
    outOpenList = "openZoneMB.txt"
    
    IORobot.readContigOut(folderName, mummerLink, graphFileName, contigFile, outContigFile, outOpenList)
    
    
    # ## Repeat resolution  [Proxy for phasing step]
    # 6. Find out the repeat region by MSA
    # 7. Find out the location of SNPs and extend across repeat 
    # [short cut : use contig creator : your job here is to get data into the correct formats]
    
    
    
    
    print "xPhased"


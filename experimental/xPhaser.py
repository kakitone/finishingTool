import commonLib
from itertools import groupby
from operator import itemgetter
# ## A new X-phasing step

# ## Helper functions 

# ## Disjoint Union Data Structure
class clusterElem(object):
    def __init__(self, index):
        self.rank = 0
        self.parent = self
        self.id = index
        self.terminatingLoc = -1
        self.childList = []
        self.voteList = []
    def findSuccessor(self):
        self.voteList.sort()
        curMax = -1
        curVal = -1
        for key, items in groupby(self.voteList):
            maxCnt = 0 
            for eachitem in items:
                maxCnt += 1
            if maxCnt > curMax:
                maxCnt = 0
                curVal = key
        return curVal
            
        
def find(x):
    if x.parent == x:
        return x
    else:
        return find(x.parent)
       
def union(x, y):
    
    xRoot = find(x)
    yRoot = find(y)
    
    if xRoot == yRoot:
        return 0
    
    if xRoot.rank < yRoot.rank:
        xRoot.parent = yRoot
        yRoot.childList.append(xRoot)
        
    elif xRoot.rank > yRoot.rank:
        yRoot.parent = xRoot
        xRoot.childList.append(yRoot)
    else:
        yRoot.parent = xRoot
        xRoot.childList.append(yRoot)
        xRoot.rank = xRoot.rank + 1 
        
    return 1

def familyList(x):
    root = find(x)
    stack = []
    familyKmers = [root]
    
    stack.append(root)
    while (len(stack) > 0):
        item = stack.pop(0)
        for eachsubitem in item.childList:
            familyKmers.append(eachsubitem)
            stack.append(eachsubitem)
            
    return familyKmers
                    
def parseContigName(contigName, side):
    myInfo = contigName[6:].split('_')
    index = 4 * int(myInfo[0])
    if myInfo[1] == 'd':
        index = index + 2
    
    if side == 'R':
        index = index + 1 
    
    return index    

def parseIDToName(id):

    index = id / 4
    name = "Contig" + str(index) + "_"
    reminder = id % 4
    if reminder == 0 :
        name = name + "p_L"
    elif reminder == 1:
        name = name + "p_R"
    elif reminder == 2:
        name = name + "d_L"
    elif reminder == 3:
        name = name + "d_R"
    return name 

def checkBasicRequirement(recordFromMummer, percentThres, lenThres):
    start1, end1, start2, end2, match1, match2, percentMatch, name1, name2 = recordFromMummer
    ck = True
    if name1 == name2:
        ck = False
    
    if start2 > end2: 
        ck = False
        
    if percentMatch < percentThres:
        ck = False
        
    if min(match1, match2) < lenThres:
        ck = False
    return ck
    
def checkSameSideRequirement(recordFromMummer, lenDic):
    resultOfCk = 'N'
    ck = True
    
    thresLarge = 50000 
    thresSmall = 7
    percentThres = 50
    lenThres = 100
    
    terminatingLoc = []
    
    # Format : [1, 1525398, 1, 1525398, 1525398, 1525398, 100.0, 'Contig0_d', 'Contig0_d']
    start1, end1, start2, end2, match1, match2, percentMatch, name1, name2 = recordFromMummer
    
    ck = checkBasicRequirement(recordFromMummer, percentThres, lenThres)
    
    if min(start1, start2) < thresSmall and max(start1, start2) < thresLarge:
        resultOfCk = 'L'
        terminatingLoc = [end1, end2]
        
    
    if min(lenDic[name1] - end1, lenDic[name2] - end2) < thresSmall and max(lenDic[name1] - end1, lenDic[name2] - end2) < thresLarge:
        resultOfCk = 'R'
        terminatingLoc = [start1, start2]
    
    if ck:
        return terminatingLoc, resultOfCk
    else:
        return [], 'N'
    
def checkOppositeSideRequirement(recordFromMummer, lenDic):
    
    thresSmall = 7
    percentThres = 50
    lenThres = 100

    # Format : [1, 1525398, 1, 1525398, 1525398, 1525398, 100.0, 'Contig0_d', 'Contig0_d']
    start1, end1, start2, end2, match1, match2, percentMatch, name1, name2 = recordFromMummer
    
    ck = checkBasicRequirement(recordFromMummer, percentThres, lenThres)
    
    if ck and min(start1, start2) < thresSmall and min(lenDic[name1] - end1, lenDic[name2] - end2) + 1 < thresSmall:
        if start1 < thresSmall:
        
            return True, [name2, name1, start2, end1]
        else:
            return True , [name1, name2, start1, end2]
    else:
        return False, []
        
def connectContigs(toPhase, toRemove, toBR, folderName, mummerLink):
    print "\nConnect Contigs"
    tmpList = []
    delThres =20000
    lenDic = commonLib.obtainLength(folderName, "improved3_Double.fasta")
    
    for eachitem in toRemove:
        tmpList.append(eachitem/2)
        
    tmpList.sort()
    removeContigIndexList = []
    for key, items in groupby(tmpList):
        name = "Contig"+ str(key)+"_p"
        if lenDic[name] < delThres:
            removeContigIndexList.append(2*key)
            removeContigIndexList.append(2*key+1)
        
    print "removeContigIndexList", removeContigIndexList
    
    
    
    ### toRemove ===> remove both strand when detected
    G = commonLib.seqGraph(len(lenDic))
    
    ### hack ! make the nodeIndexList to be empty for empty nodes

    for eachnode in G.graphNodesList:
        if eachnode.nodeIndex in removeContigIndexList:
            eachnode.nodeIndexList = []
    
    # form a graph, .condense, then use readContigOut
    ### add edge 
    for eachedge in toBR:
        i = eachedge[0]/2
        j = eachedge[1]/2
        wt = eachedge[3]+1
        print "i, j, wt", i, j, wt
        G.insertEdge(i, j, wt)
    
    tmpFileName = "xphasebonus"
    G.condense()
    G.saveToFile(folderName,tmpFileName )
    
    commonLib.readContigOut(folderName, mummerLink, tmpFileName, "improved3_Double.fasta", "improved4.fasta", tmpFileName+"Open")
    
    
def obtainContigStatistics(folderName):
    print "defineRegionOfInterest"
    # commonLib.writeToFile_Double1(folderName, "improved3.fasta", "improved3_Double.fasta", "contig")
    # commonLib.useMummerAlign(mummerLink, folderName, "phasing", "improved3_Double.fasta", "improved3_Double.fasta")
    
    dataList = commonLib.extractMumData(folderName, "phasing" + "Out")
    lenDic = commonLib.obtainLength(folderName, "improved3_Double.fasta")
    
    print "Record length of contigs"
    for eachitem in lenDic:
        print lenDic[eachitem], eachitem    
        
    return dataList, lenDic

def gatherAssociatedContigs(folderName, lenDic,dataList):
    print "\nPerform alignment and group associated contigs"
    # Convention : 0_p_L, 0_p_R, 0_d_L, 0_d_R 
    N = len(lenDic) * 2
    clusterList = []
    for i in range(N):
        clusterList.append(clusterElem(i))
        if i % 2 == 0:
            clusterList[i].terminatingLoc = 0
        else:
            clusterList[i].terminatingLoc = lenDic[parseIDToName(i)[0:-2]]
    
    oppoPairList = []
    for eachitem in dataList:
        terminatingLoc, resultOfCk = checkSameSideRequirement(eachitem, lenDic)
        isOppMatch, pair = checkOppositeSideRequirement(eachitem, lenDic)
        
        if isOppMatch:
            index1 = parseContigName(pair[0], 'R')
            index2 = parseContigName(pair[1], 'L')
            oppoPairList.append([index1, index2, pair[2], pair[3]])
            
        
        if resultOfCk == 'L' or resultOfCk == 'R':
            index1 = parseContigName(eachitem[-2], resultOfCk)
            index2 = parseContigName(eachitem[-1], resultOfCk)
            
            union(clusterList[index1], clusterList[index2])
            
            if resultOfCk == 'L':
                if clusterList[index1].terminatingLoc < terminatingLoc[0]:
                    clusterList[index1].terminatingLoc = terminatingLoc[0]
                
                if clusterList[index2].terminatingLoc < terminatingLoc[1]:
                    clusterList[index2].terminatingLoc = terminatingLoc[1]
            
            elif resultOfCk == 'R':
                if clusterList[index1].terminatingLoc > terminatingLoc[0]:
                    clusterList[index1].terminatingLoc = terminatingLoc[0]
                
                if clusterList[index2].terminatingLoc > terminatingLoc[1]:
                    clusterList[index2].terminatingLoc = terminatingLoc[1]
            

    headList = []
    for eachitem in clusterList:
        if find(eachitem) == eachitem:
            headList.append(eachitem)
    
    for eachitem in headList:
        
        for eachsub in familyList(eachitem):
            print parseIDToName(eachsub.id), eachsub.terminatingLoc,
        print 
   
    nFamily = len(headList)
    
    # Define the match of inInfo vs outInfo [matchList]
    oppoPairList.sort()
    for key, items in groupby(oppoPairList, itemgetter(0, 1)):
        
        # print parseIDToName(key[0]), parseIDToName(key[1])
        find(clusterList[key[0]]).voteList.append(key[1])
        find(clusterList[key[1]]).voteList.append(key[0])
    
    
    matchList = []
    for eachitem in headList:
        if eachitem.id % 2 == 1:
            successorIndex = eachitem.findSuccessor()
            if successorIndex != -1:
                matchList.append([eachitem.id, successorIndex])
            
    repeatList = []
    for eachitem in  matchList:
        if eachitem[0] != -1 and eachitem[1] != -1:
            
            inList = []
            for eachsubitem in familyList(find(clusterList[eachitem[0]])):
                inList.append([eachsubitem.id, eachsubitem.terminatingLoc])
                
            outList = []
            for eachsubitem in familyList(find(clusterList[eachitem[1]])):
                outList.append([eachsubitem.id, eachsubitem.terminatingLoc])
            
             
            repeatList.append([inList, outList])
            
    # Filter the embedded contigs
    globalRemoveList = []
    
    for eachitem in repeatList:
        inList = eachitem[0]
        outList = eachitem[1]
        
        print "len(inList), len(outList)", len(inList), len(outList)
        print inList 
        print outList
        
        toRemoveList = []
        for aitem in inList:
            for bitem in outList:
                if aitem[0] / 2 == bitem[0] / 2 :
                    toRemoveList.append([aitem, bitem])
                    globalRemoveList.append(aitem[0] / 2)
                
        for eachsub in toRemoveList:
            if eachsub[0] in inList:
                inList.remove(eachsub[0])
            if eachsub[1] in outList:
                outList.remove(eachsub[1])
   
     
    print "\nRepeats and in/out contigs"
    for i in range(len(repeatList)):
        print "(repeatList[i][0]),(repeatList[i][1]): ", (repeatList[i][0]), (repeatList[i][1]) 
    print "globalRemoveList: ", globalRemoveList
    print "oppoPairList", oppoPairList
    
    
    return repeatList, globalRemoveList, oppoPairList


def defineRepeatAndFlankingContigs(folderName, dataList, lenDic, repeatList, globalRemoveList, oppoPairList):
    ### Here is the problem ... should define transitive connectivity to define BR... otherwise, dangerous...
    ### Or one can just keep them 
    
    # Define the repeat 
    contigDic = commonLib.loadContigsFromFile(folderName, "improved3_Double.fasta") 
    newRepeatList = []
    newBRList = []
    
    
    
    print "\nRepeat interior and defining flanking region"
    for eachrepeat in repeatList:
        # ## Get the initial trial
        inReadList = []
        outReadList = []
        for eachitem in eachrepeat[0]:
            inReadList.append(eachitem[0])
        for eachitem in eachrepeat[1]:
            outReadList.append(eachitem[0])
        
        tmpLink = []
        for eachoppoPair in oppoPairList:
            if eachoppoPair[0] in inReadList and eachoppoPair[1] in outReadList:
                tmpLink = eachoppoPair
                break
            
        if len(tmpLink) > 0:
            f1Read, f2Read = tmpLink[0], tmpLink[1]
            f1 , a1, f2, a2 = -1, tmpLink[2], -1 , tmpLink[3]
            
            for eachitem in eachrepeat[0]:
                if eachitem[0] == f1Read:
                    f1 = eachitem[1]
            for eachitem in eachrepeat[1]:
                if eachitem[0] == f2Read:
                    f2 = eachitem[1]
            
            print "f1Read, f2Read, f1, a1, f2, a2:\t ", f1Read, f2Read, f1, a1, f2, a2
            
            f1tilde , f2tilde = f1, f2
            # ## Refine it
            for myrecord in dataList:
                myid = parseContigName(myrecord[-2], 'R')       
                otherid = parseContigName(myrecord[-1], 'R')
                
                if myid == f1Read and otherid != myid and otherid in inReadList:
                    if checkSameSideRequirement(myrecord, lenDic):
                        myStart = myrecord[0] 
                        if myStart > f1tilde:
                            f1tilde = myStart
            
            for myrecord in dataList:
                myid = parseContigName(myrecord[-2], 'L')       
                otherid = parseContigName(myrecord[-1], 'L')
                
                if myid == f2Read and otherid != myid and otherid in outReadList:
                    if checkSameSideRequirement(myrecord, lenDic):
                        myEnd = myrecord[1] 
                        if myEnd < f2tilde:
                            f2tilde = myEnd
                            
            # ## Output the loc indices and read from the real contig to get the repeat out
            print "f1Read, f2Read, f1tilde, a1, f2tilde, a2:  \t", f1Read, f2Read, f1tilde, a1, f2tilde, a2
            
            
            f1Read_parsed = parseIDToName(f1Read)[0:-2]
            f2Read_parsed = parseIDToName(f2Read)[0:-2]
            
            print "f1Read_parsed, f2Read_parsed", f1Read_parsed, f2Read_parsed, lenDic[f1Read_parsed], lenDic[f2Read_parsed]
            
            if a2 < f2tilde:
                repeatSegment = contigDic[f1Read_parsed][f1tilde:] + contigDic[f2Read_parsed][a2:f2tilde]
            else:
                repeatSegment = contigDic[f1Read_parsed][f1tilde:f2tilde - a2]
            
            
            print "len(repeatSegment)", len(repeatSegment)
            # ## Put to repeat, remove from repeat, add toBR 
            
            if a2 >= f2tilde  or a1<= f1tilde:
                
                newBRList.append([f1Read, f2Read, a1, a2])
                tmpSeg = [[], [], repeatSegment]
                for eachin in eachrepeat[0]:
                    if eachin[0] != f1Read:
                        tmpSeg[0].append(eachin)
                
                for eachout in eachrepeat[1]:
                    if eachout[0] != f2Read:
                        tmpSeg[1].append(eachout)
                
                
                if len(tmpSeg[0]) == 1 and len(tmpSeg[1]) == 1:
                    # TODO
                    inIndex = tmpSeg[0][0][0]
                    outIndex = tmpSeg[1][0][0]
                    found = False
                    
                    # if repeat exists , then fill in the blanks
                    # otherwise, fill 0, 0. 
                    a1New, a2New = -1, 0
                    for  secondrecord in dataList:
                        isOppMatch, pair = checkOppositeSideRequirement(secondrecord, lenDic)
                        
                        if isOppMatch:
                            index1 = parseContigName(pair[0], 'R')
                            index2 = parseContigName(pair[1], 'L')
                            
                            if index1 == inIndex and index2 == outIndex:
                                oppoPairList.append([index1, index2, pair[2], pair[3]])
                        
                    if not found:
                        newBRList.append([inIndex, outIndex, -1, 0])
                    else:
                        newBRList.append([inIndex, outIndex, a1New, a2New])
                        
                elif len(tmpSeg[0]) >= 1 and len(tmpSeg[1]) >= 1:
                    newRepeatList.append(tmpSeg) 
            else:
                
                if len(eachrepeat[0]) == 1 and len(eachrepeat[1]) == 1:
                    newBRList.append([f1Read, f2Read, a1, a2])
                elif len(eachrepeat[0]) >= 1 and len(eachrepeat[1]) >= 1:
                    newRepeatList.append([ repeatSegment, eachrepeat[0], eachrepeat[1]])
            
    # Format output 
    # Rmk: if only 1 copy is left, addToBR; if 0 left, remove that repeat    
    toPhase = newRepeatList
    toRemove = globalRemoveList
    toBR = newBRList
    
    print "Items to be returned to next step:"
    print "toPhase", len(toPhase), toPhase
    print "toRemove", len(toRemove), toRemove
    print "toBR", len(toBR) , toBR
    
    return toPhase, toRemove, toBR

# 1 : Find repeat segment and flanking region
    # a ) Use tip info to find related reads 
    # b ) Perform MSA
    # c ) Find flanking region and repeat region
    
# ## Input : improved3_double.fasta
# ## Output : toPhase = [[repeat , inInfo, outInfo]]
# ##          toRemove = [contigNums]
# ##          toBR = [ [inContigName, outContigName, anchorIn, anchorOut]]

# ## repeat = (ACTTACC) if exist. 
# ## inInfo[0] = [contigName(contig1_p), endUsed('R'), terminating loc(1999997)]
# ## outInfo[0] = [contigName(contig2_p), endUsed('L'), terminating loc(50)]
# ## Need to put all ends of all contigs into such a format 






def defineRegionOfInterest(folderName , mummerLink):
    

    dataList, lenDic= obtainContigStatistics(folderName)
    repeatList, globalRemoveList, oppoPairList = gatherAssociatedContigs(folderName, lenDic,dataList)
    toPhase, toRemove, toBR = defineRepeatAndFlankingContigs(folderName, dataList, lenDic, repeatList, globalRemoveList, oppoPairList)
    connectContigs(toPhase, toRemove, toBR, folderName, mummerLink)
  
  

# 2 : Find sites of polymorphism and anchors
    # a) Determine sites of polymorphism [both indel and sub]
    # b) Choose anchors to use 
        # i) Sliding window to find range
        # ii) k-mean to find strong classifiers
        # iii) decision tree to find fine classifiers
def findSitesOfPolymorphism(folderName , mummerLink):
    print "\nfindSitesOfPolymorphism"



# 3 : Do bridging based on anchors and from new contigs 
    # a) Find the associated raw reads
    # b) Establish the vote for each read
    # c) Run Hungarian and do threshold check
    # d) Merge and extend across
def bridgingAcross(folderName , mummerLink):
    print "\nbridgingAcross" 
    

    
    
def mainFlowOfXPhaser(folderName , mummerLink):
    defineRegionOfInterest(folderName , mummerLink)
    findSitesOfPolymorphism(folderName , mummerLink)
    bridgingAcross(folderName , mummerLink)



folderName , mummerLink = "Ph1/", "MUMmer3.23/"
mainFlowOfXPhaser(folderName , mummerLink)

import cleaner
import loggingIndel
import math

 
def findVoteScore(in1List, leftSeg, parameterRobot):
    H0Score, H1Score = 0,0 
    threshold = 10
    perr = parameterRobot.p
    
    print " len(in1List) ",  len(in1List)
    
    for i in range(len(in1List)):
        # Determine F/B

        read = in1List[i].longread
        
        
        score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignment(read, leftSeg[0], parameterRobot)       
        
        reverseread = cleaner.reverseStrand(read)
        scoreRev, returnalignedSeq1Rev, returnalignedSeq2Rev , startiRev, startjRev , endiRev, endjRev = cleaner.SWAlignment(reverseread, leftSeg[0], parameterRobot)
        
        if score > scoreRev:
            forwardCk = True
            scoreLeft0 = score
        else:
            forwardCk =False
            scoreLeft0 = scoreRev
            
        if forwardCk :
            score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignment(read, leftSeg[1], parameterRobot)
            scoreLeft1 = score
            
        else:
            score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignment(reverseread, leftSeg[1], parameterRobot)
            scoreLeft1 = score
            
        print "scoreLeft0, scoreLeft1     " , scoreLeft0, scoreLeft1     
            
        skipCounting = False
        MAPScore= [0,0]
        # Filering reads 

        if startj > threshold or  endj < len(leftSeg[0]) - threshold:
            skipCounting = True

        # MAP decision rule 
        if forwardCk : 
            for i in range(2):
                #tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath= cleaner.SWAlignmentFixRefMP(leftSeg[i],read,  parameterRobot)
                #lengthSeed = len(leftSeg[i])
                #countConfirm, countDel, countIns, countSub = cleaner.countEdits(returnalignedSeq2, returnalignedSeq1) # care for the 1,2 and its fcn
                #confirmNoDelete,confirmNoInsert = countConfirm ,  lengthSeed - countIns
                #MAPScore[i] = math.log(1-perr)*confirmNoDelete + math.log(1-perr)*confirmNoInsert + math.log(perr)*countDel + math.log(perr*perr/3)*countSub + math.log(perr/3)*countIns + math.log(numPath)
                tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath= cleaner.SWAlignmentFixRefMPQS(leftSeg[i],read,  parameterRobot)
                MAPScore[i] =tempScore
        
        else:
            for i in range(2):
                #tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath= cleaner.SWAlignmentFixRefMP(leftSeg[i],reverseread,  parameterRobot)
                #lengthSeed = len(leftSeg[i])
                #countConfirm, countDel, countIns, countSub = cleaner.countEdits(returnalignedSeq2, returnalignedSeq1) # care for the 1,2 and its fcn
                #confirmNoDelete,confirmNoInsert = countConfirm ,  lengthSeed - countIns
                #MAPScore[i] = math.log(1-perr)*confirmNoDelete + math.log(1-perr)*confirmNoInsert + math.log(perr)*countDel + math.log(perr*perr/3)*countSub + math.log(perr/3)*countIns + math.log(numPath)
                tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath= cleaner.SWAlignmentFixRefMPQS(leftSeg[i],reverseread,  parameterRobot)
                MAPScore[i] = tempScore
    

        
        if not skipCounting:
            H0Score += MAPScore[0] 
            H1Score += MAPScore[1]

            
    return H0Score , H1Score

def decideHypothesis(in1List, in2List, leftSeg, parameterRobot):
    leftMap  = []
    H0, H1 = 0 ,0 
    H0Score, H1Score = findVoteScore(in1List, leftSeg, parameterRobot)
    H0 = H0 + H0Score
    H1 = H1 + H1Score
    
    H1Score, H0Score = findVoteScore(in2List, leftSeg, parameterRobot)
    H0 = H0 + H0Score
    H1 = H1 + H1Score
    
    if H0 > H1:
        leftMap= [0,1]
    else:
        leftMap = [1,0]
        
        
    print "H0, H1", H0, H1
    
    return leftMap, H0, H1

def decideMiddlePiece(commonList, leftSeg, rightSeg, parameterRobot):
    middleMap = []
    threshold = 10 
    perr = parameterRobot.p
    
    H0 , H1 = 0 ,0 
    for i in range(len(commonList)):
        read = commonList[i].longread
        
        score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignment(read, leftSeg[0], parameterRobot)
        readRev = cleaner.reverseStrand(read)
        scoreRev, returnalignedSeq1Rev, returnalignedSeq2Rev , startiRev, startjRev , endiRev, endjRev = cleaner.SWAlignment(readRev, leftSeg[0], parameterRobot)
        
        if score > scoreRev :
            forwardCk = True
            scoreL0 = score
        else:
            forwardCk = False
            scoreL0 = scoreRev
        
        if forwardCk :
            tmpread = read
        else:
            tmpread = readRev
            
        scoreL1, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignment(tmpread, leftSeg[1], parameterRobot)
        scoreR0 , returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignment(tmpread, rightSeg[0], parameterRobot)
        scoreR1 , returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignment(tmpread, rightSeg[1], parameterRobot)
                            
        skipCounting = False

        # Filering reads 

        if startj > threshold or  endj < len(leftSeg[0]) - threshold:
            skipCounting = True

        # MAP decision rule 
        if not skipCounting:
            
            MAPScoreL= [0,0]
            for i in range(2):
                #tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath= cleaner.SWAlignmentFixRefMP(leftSeg[i],tmpread,  parameterRobot)
                #lengthSeed = len(leftSeg[i])
                #countConfirm, countDel, countIns, countSub = cleaner.countEdits(returnalignedSeq2, returnalignedSeq1) # care for the 1,2 and its fcn
                #confirmNoDelete,confirmNoInsert = countConfirm ,  lengthSeed - countIns
                #MAPScoreL[i] = math.log(1-perr)*confirmNoDelete + math.log(1-perr)*confirmNoInsert + math.log(perr)*countDel + math.log(perr*perr/3)*countSub + math.log(perr/3)*countIns + math.log(numPath)
                tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath= cleaner.SWAlignmentFixRefMPQS(leftSeg[i],tmpread,  parameterRobot)
                MAPScoreL[i] =tempScore
                
            MAPScoreR= [0,0]
            for i in range(2):
                #tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath= cleaner.SWAlignmentFixRefMP(rightSeg[i],tmpread,  parameterRobot)
                #lengthSeed = len(leftSeg[i])
                #countConfirm, countDel, countIns, countSub = cleaner.countEdits(returnalignedSeq2, returnalignedSeq1) # care for the 1,2 and its fcn
                #confirmNoDelete,confirmNoInsert = countConfirm ,  lengthSeed - countIns
                #MAPScoreR[i] = math.log(1-perr)*confirmNoDelete + math.log(1-perr)*confirmNoInsert + math.log(perr)*countDel + math.log(perr*perr/3)*countSub + math.log(perr/3)*countIns + math.log(numPath)
                tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath= cleaner.SWAlignmentFixRefMPQS(rightSeg[i],tmpread,  parameterRobot)
                MAPScoreR[i] = tempScore
            
            H0Score  =  max( MAPScoreL[0] + MAPScoreR[0] , MAPScoreL[1] + MAPScoreR[1])
            H1Score =  max( MAPScoreL[1] + MAPScoreR[0], MAPScoreL[0] + MAPScoreR[1]) 
            
            H0  = H0 + H0Score
            H1 = H1 + H1Score 
            
            
        
    
    if H0 > H1:
        middleMap = [0,1]
    else:
        middleMap = [1,0]
        
    print "Middle : H0, H1", H0, H1
    
    return middleMap , H0, H1

def linkSegments(inSeg, middleSeg, outSeg, parameterRobot):
    
    print "\nNew Contig"

    middleRead = middleSeg.longread
    revmiddleRead = cleaner.reverseStrand(middleRead)
  
    extendedRead = []

    score, scoreRev = -1, -1 
    inIndex = 0
    startjTmp = -1
    
    while ((max(score, scoreRev)< 30 or startjTmp> 10) and inIndex < len(inSeg)): 
        inRead = inSeg[inIndex].longread
        revinRead = cleaner.reverseStrand(inRead)
        score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj  = cleaner.SWAlignment(inRead, middleRead, parameterRobot)
        scoreRev, returnalignedSeq1Rev, returnalignedSeq2Rev , startiRev, startjRev , endiRev, endjRev  = cleaner.SWAlignment(revinRead, middleRead, parameterRobot)
        if score > scoreRev:
            startjTmp = startj
            print "starti, startj , endi, endj", starti, startj , endi, endj
        else:
            startjTmp = startjRev
            print "startiRev, startjRev , endiRev, endjRev", startiRev, startjRev , endiRev, endjRev
            
        print "startjTmp", startjTmp
            
        print "score, scoreRev,inIndex " ,score, scoreRev  ,inIndex

        inIndex += 1
        
    
    if score  > scoreRev : 
        for eachbase in inRead[0:starti]:
            extendedRead.append(eachbase)
            
        print "len(extendedRead)",  len(extendedRead)
        
        for eachbase in middleRead:
            extendedRead.append(eachbase)
            
        print "len(extendedRead)", len(extendedRead)
    else: 
        
        for eachbase in revinRead[0:startiRev]:
            extendedRead.append(eachbase)
        
        print "len(extendedRead)",  len(extendedRead)
        for eachbase in middleRead:
            extendedRead.append(eachbase)
            
        print "len(extendedRead)",  len(extendedRead)

    score, scoreRev = -1, -1 
    # hack
    outIndex = 0
    endiTmp = 1000
    outRead, revoutRead = [] , [] 
    while ((max(score, scoreRev)< 30  or endiTmp < len(middleRead)-10)and outIndex < len(outSeg)):     
        
        
        outRead = outSeg[outIndex].longread
        revoutRead = cleaner.reverseStrand(outRead)
        score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj  = cleaner.SWAlignment(middleRead,outRead , parameterRobot)
        scoreRev, returnalignedSeq1Rev, returnalignedSeq2Rev , startiRev, startjRev , endiRev, endjRev  = cleaner.SWAlignment(middleRead, revoutRead, parameterRobot)
        
        
        if score  > scoreRev : 
            endiTmp = endi
            print "starti, startj , endi, endj", starti, startj , endi, endj
        else:
            endiTmp = endiRev
            print " startiRev, startjRev , endiRev, endjRev ",  startiRev, startjRev , endiRev, endjRev 
            
        print "score, scoreRev, outIndex" ,score, scoreRev, outIndex
          
        
        outIndex += 1 
    
    
    correctionTerm = len(middleRead) - endiTmp  
    if score  > scoreRev : 
        for eachbase in outRead[endj+correctionTerm:len(outRead)]: 
            extendedRead.append(eachbase)
            
        print "len(extendedRead) forw", len(extendedRead)

    else:
        for eachbase in  revoutRead[endjRev+correctionTerm:len(revoutRead)]:
            extendedRead.append(eachbase)
           
        print "len(extendedRead) rev", len(extendedRead)

    return extendedRead


def shortenSeg(segList, lookRange):
    newSegmentList = []

    indexOfDifference = -1 
    for i in range(len(segList[0])):
        if segList[0][i] != segList[1][i]:
            indexOfDifference = i
    
    startIndex, endIndex = max(0, indexOfDifference -lookRange ), min(indexOfDifference + lookRange +1, len(segList[0]))
    segList = [segList[0][startIndex:endIndex], segList[1][startIndex:endIndex]]
    
    return segList
     

def readExtender(in1List, in2List, out1List, out2List, commonList,parameterRobot,longReadToUse, debug = False):
    # Assume the incoming reads are all in forward strand for simplicity 
    #print "longReadToUse",longReadToUse
    
    if debug == False:
        for eachcommon in commonList:
            #print eachcommon.index
            if eachcommon.index == longReadToUse:
                longSeedRead =eachcommon
        #print "longSeedRead.SNPlocList", longSeedRead.SNPlocList
        
        leftSeg = longSeedRead.leftSegList
        rightSeg = longSeedRead.rightSegList
        
        print "leftSeg, rightSeg"
        print leftSeg, rightSeg
        
    else:
        leftSeg, rightSeg = loggingIndel.loadSegList(parameterRobot.defaultFolder)
        print "leftSeg, rightSeg"
        print leftSeg, rightSeg
    
    
    if len(leftSeg) == 0 or len(rightSeg) == 0 :
        return -1
    
    leftMap = []
    rightMap = []
    
    #lookRangeLbdd, lookRangeUbdd = 15, 16 
    lookRangeLbdd, lookRangeUbdd, stepSize = 9, 10 , 1 
    
    contig1, contig2 = [] , []
    

    ### Align rawLongReads to Seed Reads
    # Use Mapper 
    HTList = []
    
    difference = -1 
    for i in range(lookRangeLbdd, lookRangeUbdd, stepSize):
        lookRange = i
        leftSegRun = shortenSeg(leftSeg, lookRange )
        leftMapTmp, H0Tmp, H1Tmp  = decideHypothesis(in1List, in2List, leftSegRun, parameterRobot)
        if abs(H0Tmp - H1Tmp) > difference:
            leftMap, H0, H1  = leftMapTmp, H0Tmp, H1Tmp  
            difference =  abs(H0Tmp - H1Tmp)
    HTList.append([H0, H1])

    
    
    difference = -1 
    for i in range(lookRangeLbdd, lookRangeUbdd, stepSize):
        lookRange = i
        rightSegRun = shortenSeg(rightSeg, lookRange )
        print out1List, out2List
        
        rightMapTmp, H0Tmp, H1Tmp = decideHypothesis(out1List, out2List, rightSegRun, parameterRobot)    
        if abs(H0Tmp - H1Tmp) > difference:
            rightMap, H0, H1= rightMapTmp, H0Tmp, H1Tmp
            difference =  abs(H0Tmp - H1Tmp)      
    HTList.append([H0, H1])
    
    
    
    difference = -1 
    for i in range(lookRangeLbdd, lookRangeUbdd,stepSize):
        lookRange = i
        leftSegRun = shortenSeg(leftSeg, lookRange )
        rightSegRun = shortenSeg(rightSeg, lookRange )
        middleMapTmp , H0Tmp, H1Tmp= decideMiddlePiece(commonList, leftSegRun, rightSegRun, parameterRobot)
        if abs(H0Tmp - H1Tmp) > difference:
            middleMap , H0, H1= middleMapTmp , H0Tmp, H1Tmp
            difference =  abs(H0Tmp - H1Tmp)        
    HTList.append([H0, H1])
    
    
    
    loggingIndel.logHTVote(HTList, parameterRobot.defaultFolder)
    print "leftMap",  leftMap
    
    print "rightMap" , rightMap
    
    print "middleMap", middleMap 
    
    ### Merge seedreads based on output
    # Use some Merging 
    leftSegIndex = leftMap[0]
    rightSegIndex = middleMap[leftSegIndex]
    
    
    
    if  rightMap[0] == rightSegIndex:
        result= 0
    else:
        result= 1 
        
    if False:
        if rightMap[0] == rightSegIndex:
            if len(in1List) == 0 or len(out1List) == 0:
                contig1 = linkSegments(in2List, commonList[0], out2List, parameterRobot)
                contig2 = []
                
            elif len(in2List) == 0 or len(out2List) == 0:
                contig1 = linkSegments(in1List, commonList[0], out1List, parameterRobot)
                contig2 = []
                
            else:
                contig1 = linkSegments(in1List, commonList[0], out1List, parameterRobot)
                contig2 = linkSegments(in2List, commonList[0], out2List, parameterRobot)
    
            print "00 ,11 "
        elif rightMap[1] == rightSegIndex:
            if len(in1List) == 0 or len(out2List) == 0:
                contig1 = linkSegments(in2List, commonList[0], out1List, parameterRobot)
                contig2 = []
            elif len(in2List) == 0 or len(out1List) == 0:
                contig1 = linkSegments(in1List, commonList[0], out2List, parameterRobot)
                contig2 = []
            else: 
                contig1 = linkSegments(in1List, commonList[0], out2List, parameterRobot)
                contig2 = linkSegments(in2List, commonList[0], out1List, parameterRobot)
            print "01, 10"
            
            
        loggingIndel.logFinishedContigs([contig1, contig2], parameterRobot.defaultFolder)
        
        
    

    return result



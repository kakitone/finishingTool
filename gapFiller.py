from itertools import groupby
from operator import itemgetter
import os
import sys

import alignerRobot
import houseKeeper
import IORobot

# ## 6) Fill gap(I: improved.fasta, openZone,txt ; O: improved2.fasta )
def fillGap(folderName , mummerLink):
    print "fillGap"
    # ## Load an extra field V
    # ## migrate the existing code to here  V 
    # ## add it the functionality of checking filtered list V 
    # ## testing 
    
    formRelatedReadsFile(folderName, mummerLink)
    os.system("cp raw_reads.part* "+ folderName)
    os.system("rm raw_reads.part*")
    
    extractEdgeSet(folderName, mummerLink)
    os.system("cp relatedReads_Double.part* "+ folderName)
    os.system("rm relatedReads_Double.part*")
    
def loadOpenList(folderName):
    
    f = open(folderName + 'openZone.txt', 'r')
    tmp = f.readline().rstrip()
    N = int(tmp)
    usableJunction = [[ 0, 0 ] for i in range(N) ]
    
    tmp = f.readline().rstrip()
    while len(tmp) > 0:
        myInfo = tmp.split(',')
        contigIndex = int(myInfo[0][5:])
        if myInfo[1] == 'noprev':
            usableJunction[contigIndex][0] = 1
        if myInfo[1] == 'nonext' :
            usableJunction[contigIndex][1] = 1
             
        tmp = f.readline().rstrip()
    
    f.close() 
    
    print "usableJunction", usableJunction
    return usableJunction

def formRelatedReadsFile(folderName, mummerLink):    
    # Find associated read and extract into a file associatedReads.fasta
    # Input: contigs.fasta, cleaned_Reads.fasta 
    # Output: relatedReads.fasta

    # ## Extract heads of the contigs
    print ">formRelatedReadsFile"
    
    f = open(folderName + "improved.fasta", 'r')
    f2 = open(folderName + "improvedTrunc.fasta", 'w')
    temp = f.readline()
    tempContig = ""
    thres = 400
    runningIndex = 0
    endThres = 10 
    
    while len(temp) > 0:
        if temp[-1] == '\n':
            temp = temp[0:-1]
        
        
        if temp[0] == '>':

            if len(tempContig) > 0:
                IORobot.writeToFile(f2, runningIndex, tempContig[0:thres])
                runningIndex = runningIndex + 1
                
                IORobot.writeToFile(f2, runningIndex, tempContig[-thres:])
                runningIndex = runningIndex + 1 
                
                                
                IORobot.writeToFile(f2, runningIndex, houseKeeper.reverseComplement(tempContig[0:thres]))
                runningIndex = runningIndex + 1
                
                IORobot.writeToFile(f2, runningIndex, houseKeeper.reverseComplement(tempContig[-thres:]))
                runningIndex = runningIndex + 1
                
                tempContig = ""
        else:
            tempContig = tempContig + temp
        
        temp = f.readline()

    IORobot.writeToFile(f2, runningIndex, tempContig[0:thres])
    runningIndex = runningIndex + 1
    
    IORobot.writeToFile(f2, runningIndex, tempContig[-thres:])
    runningIndex = runningIndex + 1
                  
    IORobot.writeToFile(f2, runningIndex, houseKeeper.reverseComplement(tempContig[0:thres]))
    runningIndex = runningIndex + 1
    
    IORobot.writeToFile(f2, runningIndex, houseKeeper.reverseComplement(tempContig[-thres:]))
    runningIndex = runningIndex + 1
    
    
    f2.close()
    f.close()
    
    # ## Write double stranded reads
    IORobot.writeToFile_Double1(folderName, "improved.fasta", "improved_Double.fasta", "contig")
    # writeToFile_Double1(folderName, "raw_reads.fasta", "raw_reads_Double.fasta","read")
    
    # ## Apply MUMMER on them using cleanedReads against them
    assoiatedReadIndex = []
    nameList = []
    
    numberOfFiles = 20
    
    if True:
        bindir = os.path.abspath(os.path.dirname(sys.argv[0]))
        command = bindir + "/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + folderName + "raw_reads.fasta"
        os.system(command)
    
    
    workerList = []
    
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
            
        outputName, referenceName, queryName, specialName=  "outGapFillRaw"+indexOfMum , "improvedTrunc.fasta", "raw_reads.part-" + indexOfMum + ".fasta", "fromMum" + indexOfMum 
        workerList.append([outputName, referenceName, queryName, specialName])
    
    
    
    if True:
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, workerList, houseKeeper.globalParallel ,True)
        
        # alignerRobot.useMummerAlign(mummerLink, folderName, "out", "improvedTrunc.fasta", "raw_reads.part-" + indexOfMum + ".fasta", True, "fromMum" + indexOfMum )
        
        '''
        command = mummerLink + "nucmer --maxmatch --nosimplify -p " + folderName + "out " + folderName + "improvedTrunc.fasta raw_reads.part-" + indexOfMum + ".fasta"
        os.system(command)

        command = mummerLink + "show-coords -r " + folderName + "out.delta > " + folderName + "fromMum" + indexOfMum
        os.system(command)
        '''
        

    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
        f = open(folderName + "fromMum" + indexOfMum, 'r')
    
        for i in range(6):
            tmp = f.readline()
        
        while len(tmp) > 0:
            infoArr = tmp.split('|')
            myArr = infoArr[-1].split('\t')
            rdGpArr = infoArr[-1].split('\t')
            contigName = rdGpArr[0].rstrip().lstrip()
            readName = rdGpArr[1].rstrip().lstrip()
            
            endSegArr = infoArr[0].split(" ")
            pos = []
            for eachitem in endSegArr:
                if len(eachitem) > 0:
                    pos.append(int(eachitem))
                    
            startPos = pos[0]
            endPos = pos[1]
            if startPos < endThres and endPos > thres - endThres:
                assoiatedReadIndex.append(myArr[1])
                nameList.append([int(contigName.split('_')[1]), readName])
            tmp = f.readline()
        
        f.close()
    
    
    nameList.sort()

    assoiatedReadIndex.sort()
    
    # print "assoiatedReadIndex", assoiatedReadIndex
    
    ckIndex = 0
    f = open(folderName + "associatedNames.txt", 'w')
    oneItem = 0
    keyFound = []
    for key, items in groupby(assoiatedReadIndex):
        
        countItem = 0
        for eachitem in items:
            countItem += 1
            
        if countItem == 1:
            
            oneItem += 1
        else:
            key = key.rstrip()
            if not key in keyFound:
                f.write(key + '\n')
                keyFound.append(key)

        ckIndex += 1
    
    print "ckIndex,oneItem: ", ckIndex, oneItem
    f.close()

    fFilter = open(folderName + "associatedNames.txt", 'r')
    
    fout = open(folderName + "associatedNames2.txt", 'w') 
    
    maxCount = 12000
    mytmpDum = fFilter.readline() 
    i = 0
    while i < maxCount and len(mytmpDum) > 0:
        fout.write(mytmpDum)  
        mytmpDum = fFilter.readline() 
        i = i + 1
        
    fout.close()   
    fFilter.close()

    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' " + folderName + "associatedNames2.txt " + folderName + "raw_reads.fasta > " + folderName + "relatedReads.fasta"
    os.system(command)
    
    IORobot.writeToFile_Double1(folderName, "relatedReads.fasta", "relatedReads_Double.fasta", "read")
    

    
def extractEdgeSet(folderName, mummerLink, option="nopolish"):
    # Tasks: reconstruct the string  graph
    
    # Input : relatedReads_Double.fasta, conig_Double.fasta
    # Intermediate files: fromMum_overlap , fromMum_overlap
    # Output: connectivity of eachNode: InList, OutList [critical]
    #         connectivity of eachNode: arrow representation with size [optional]
    
    
    # ## Perform MUMMER alignment
    print ">Extract Edge set"
    contigOnlyLengthDic = IORobot.obtainLength(folderName, "improved.fasta")
    
    # print lengthDic
    lengthDic = IORobot.findContigLength(folderName, "improved")
    
    numberOfContig = len(contigOnlyLengthDic)*2

    K = 400
    thres = 5
    
    
    # ## Apply MUMMER on them using cleanedReads against them
    IORobot.truncateEndOfContigs(folderName, "improved_Double.fasta", "smaller_improvedContig.fasta", 25000, lengthDic)
    dataSet = []
    
    numberOfFiles = 20
    

    if True:
        bindir = os.path.abspath(os.path.dirname(sys.argv[0]))
        command = bindir + "/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + folderName + "relatedReads_Double.fasta"
        os.system(command)
        
        
    workerList = [] 
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
        
        outputName, referenceName, queryName, specialName=  "outGapFillRefine"+indexOfMum , "smaller_improvedContig.fasta",  "relatedReads_Double.part-" + indexOfMum + ".fasta",  "fromMumRefine" + indexOfMum
        workerList.append([outputName, referenceName, queryName, specialName])
    
        
        
    if True:
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, workerList, houseKeeper.globalParallel ,True)
        
        # alignerRobot.useMummerAlign(mummerLink, folderName, "outRefine", "smaller_improvedContig.fasta", "relatedReads_Double.part-" + indexOfMum + ".fasta", True,  "fromMumRefine" + indexOfMum)
        
    
    for dummyI in range(1, numberOfFiles + 1):
        tmpSet = IORobot.obtainLinkInfoReadContig(dummyI, mummerLink, folderName,thres, lengthDic, K)
        dataSet = dataSet + tmpSet
    
    # ## repeat aware
    usableJunction = loadOpenList(folderName)
    dataSet, blockedSet = filterRepeatEnd(dataSet, usableJunction)
    # ## repeat aware end
    
    dataSet.sort()
    matchPair = formMatchPairFromReadInfo(dataSet)
    
    # Bug fix on repeat detection from reads alone
    matchPair = filterRepeatPair(matchPair)
    # end bug fix
    
    # print matchPair

    bestMatchPair = []
    
    for key, items in groupby(matchPair, itemgetter(0, 1)):
        maxvalue = -1
        maxLenPair = []
        for eachitem in items:
            if eachitem[2] > maxvalue:
                maxvalue = eachitem[2]
                maxLenPair = [eachitem[3], eachitem[4], eachitem[5]]
        bestMatchPair.append([key[0], key[1], maxvalue, maxLenPair[0], maxLenPair[1], maxLenPair[2]])
    
    contigList, leftConnect, rightConnect, rawReadList = formbestpair(bestMatchPair,numberOfContig)
    print "contigList", contigList
    
    writeContigReadCombine(blockedSet, dataSet, folderName, rawReadList, numberOfContig, contigList, leftConnect, option, rightConnect)


def formbestpair(bestMatchPair,numberOfContig):
    rawReadList = []
    
    bestMatchPair.sort(key=lambda tup: tup[2], reverse=True)
    print "bestMatchPair", bestMatchPair

    leftConnect = [[-1, -1, -1] for i in range(numberOfContig)]
    rightConnect = [[-1, -1, -1] for i in range(numberOfContig)]
    
    for eachitem in bestMatchPair:

        prefixContig = eachitem[0]
        suffixContig = eachitem[1]
        
        if leftConnect[suffixContig][0] == -1 and rightConnect[prefixContig][0] == -1:
            leftConnect[suffixContig][0] = prefixContig 
            leftConnect[suffixContig][1] = eachitem[4]
            leftConnect[suffixContig][2] = eachitem[5]
            
            rightConnect[prefixContig][0] = suffixContig
            rightConnect[prefixContig][1] = eachitem[3]
            rightConnect[prefixContig][2] = eachitem[5]
            
            rawReadList.append(eachitem[5])
    
    startList = []
    print "leftConnect", leftConnect
    for i in range(len(leftConnect)):
        if leftConnect[i][0] == -1:
            startList.append(i)
    
    print "startList", startList
    checkLoopList = [False for i in range(len(leftConnect))]
        
    for eachitem in startList:
        tmp = eachitem
        checkLoopList[tmp] = True
        while rightConnect[tmp][0] != -1:
            tmp = rightConnect[tmp][0]
            if tmp != -1:
                checkLoopList[tmp] = True
    
    for dumdumi in range(len(leftConnect)):
        if checkLoopList[dumdumi] == False:
            startList.append(dumdumi)
            tmp = dumdumi
            
            while checkLoopList[tmp] == False :
                if tmp != -1:
                    checkLoopList[tmp] = True
                    
                tmp = rightConnect[tmp][0]


    contigList = []
    print "startList", startList
    for eachitem in startList:
        tmp = eachitem
        myList = [tmp]
        mystart = tmp
        while rightConnect[tmp][0] != -1 and rightConnect[tmp][0] != mystart:
            tmp = rightConnect[tmp][0]
            if tmp != -1:
                myList.append(tmp)
        contigList.append(myList)
    
    return contigList, leftConnect, rightConnect, rawReadList


def writeContigReadCombine(blockedSet, dataSet, folderName, rawReadList, numberOfContig, contigList, leftConnect, option, rightConnect):
    # ## repeat aware logging
    # print "myExtraLinkList", myExtraLinkList
    # ## end repeat aware logging

    myExtraLinkList = loggingReadsToRepeat(blockedSet + dataSet, contigList)    
    i = 0
    fOriginal = open(folderName + "improved.fasta", 'r')
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) > 0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead + tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)  
    fOriginal.close()
    
    # ## Put the needed rawReads into the RAM using Dictionary
    
    fAppendRaw = open(folderName + "appendRaw.txt", 'w')
    for eachraw in rawReadList:
        fAppendRaw.write(eachraw)
        fAppendRaw.write('\n')
    fAppendRaw.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' " + folderName + "appendRaw.txt " + folderName + "relatedReads_Double.fasta > " + folderName + "rawToAppend.fasta"
    os.system(command)
    
    rawRead = {}
    
    fOriginal = open(folderName + "rawToAppend.fasta", 'r')
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    tmpName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            
            if len(tmpRead) > 0:
                rawRead[tmpName] = tmpRead
                tmpRead = ""
                
            tmpName = tmp[1:]
        else:
            tmpRead = tmpRead + tmp
            
        tmp = fOriginal.readline().rstrip()
        
    rawRead[tmpName] = tmpRead
    # ## End
    
    seqToPrint = []
    contigUsed = [False for i in range(numberOfContig / 2)]
    storedStrand = [[-1, 'n'] for i in range(numberOfContig)]
    
    finalList = []
    for eachContig, i in zip(contigList, range(len(contigList))):
        tmpList = []
        for eachitem in eachContig:

            readNum = eachitem / 2
            if contigUsed[readNum] == False:
                seqToPrint.append(eachitem)
                tmpList.append(eachitem)
                contigUsed[readNum] = True
                # ## mark ouput strandinfo
                storedStrand[eachitem] = [len(finalList), 'p']
                
                
        if len(tmpList) > 0:
            finalList.append(tmpList)
    
    
    for kkk in range(len(storedStrand)):
        if storedStrand[kkk][1] == 'n':
            if kkk % 2 == 0:
                storedStrand[kkk][0] = storedStrand[kkk + 1][0]
                storedStrand[kkk][1] = 'd'
            else:
                storedStrand[kkk][0] = storedStrand[kkk - 1][0]
                storedStrand[kkk][1] = 'd'
    
    # ## begin stored output blocked pairs
    blockExtraStored(storedStrand, myExtraLinkList, folderName)
    # ## end output blocked pairs stored
    

    fImproved = open(folderName + "improved2.fasta", 'w')
    
    for eachcontig, dummyIndex in zip(finalList, range(len(finalList))):
        fImproved.write(">Segkk" + str(dummyIndex) + '\n')
        tmpStore = -1997
        tmpStore2 = -1998
        tmpStore3 = -1999
        
        for eachseg, hidum in zip(eachcontig, range(len(eachcontig))):
            readNum = eachseg / 2
            orientation = eachseg % 2
            
            newStart = 0 
            
            x , y , l = tmpStore, leftConnect[eachseg][1], tmpStore2
            extraRead = ""
            if hidum == 0:
                newStart = 0
            else:
                
                if l < x + y:
                    newStart = x + y - l
                else:
                    newStart = 0
                    if option == 'polish':
                        print "Missing polish"
                        extraRead = tmpStore3[x:l - y]
                        # extraRead = performPolishing(leftConnect[eachseg][0], eachseg, tmpStore3[x:l-y],  dataSet, folderName)
                    else:
                        extraRead = tmpStore3[x:l - y]
    
            print extraRead[0:10], len(extraRead)
            
            fImproved.write(extraRead)
            
            if orientation == 0:
                fImproved.write(readSet[readNum][newStart:])   

            else:
                fImproved.write(houseKeeper.reverseComplement(readSet[readNum])[newStart:])
            
            if rightConnect[eachseg][1] != -1:
                tmpStore = rightConnect[eachseg][1]
                tmpStore2 = len(rawRead[rightConnect[eachseg][2]])
                tmpStore3 = rawRead[rightConnect[eachseg][2]]
                
        fImproved.write('\n')
        
    fImproved.close()



def filterRepeatEnd(dataSet, usableJunction):
    newDataSet = []
    blockedSet = []
    # Format : dataSet.append((readName, contigName, 'R',matchLen)) 
    for eachitem in dataSet:
        readName, contigName, seekDir, matchLen = eachitem
        # print contigName
        myInfo = contigName[6:].split('_')
        index, orient = int(myInfo[0]), myInfo[1]
        
        # print "index, orient", index, orient
        check = False
        if orient == 'p':
            if usableJunction[index][0] == 1 and seekDir == 'L':
                check = True
            elif usableJunction[index][1] == 1 and seekDir == 'R':
                check = True
                
        elif orient == 'd':
            if usableJunction[index][1] == 1 and seekDir == 'L':
                check = True
            elif usableJunction[index][0] == 1 and seekDir == 'R':
                check = True

        if check:
            newDataSet.append(eachitem)
        else:
            blockedSet.append(eachitem)
    


    return newDataSet, blockedSet


def filterRepeatPair(matchPair):
    newMatchPair = []
    # Format : [[3, 5, 2486, 2532, 2486, 'Read48_d'] ]
    inNoList = []
    outNoList = []
    
    matchPair.sort(key = itemgetter(0))
    for key, items in groupby(matchPair, itemgetter(0)):
        ct = 0
        anotherSideList = []
        for eachitem in items:
            ct = ct +1 
            anotherSideList.append(eachitem[1])
        
        if len(set(anotherSideList)) > 1 :
            inNoList.append(key)
    
    matchPair.sort(key= itemgetter(1))
    for key, items in groupby(matchPair, itemgetter(1)):
        ct = 0
        anotherSideList = []
        for eachitem in items:
            ct = ct +1
            anotherSideList.append(eachitem[0])
         
        if len(set(anotherSideList)) > 1 :
            outNoList.append(key)
    
    for eachitem in matchPair:
        if not eachitem[0] in inNoList and not eachitem[1] in outNoList:
            newMatchPair.append(eachitem)
    
    return newMatchPair


def loggingReadsToRepeat(blockedSet, contigList):
    
    # for eachitem in blockedSet:
    #    print eachitem
    
    matchPair = formMatchPairFromReadInfo(blockedSet)
    matchPair.sort()
    
    print "loggingReadsToRepeat" 
    
    startList = []
    endList = []
    
    returnList = []
    
    for eachitem in contigList:
        startList.append(eachitem[0])
        endList.append(eachitem[-1])
    
    # print startList, endList
    for eachitem in matchPair: 
        mystart = eachitem[0]
        myend = eachitem[1]
        len1, len2 = eachitem[3], eachitem[4]
        
        
        if myend in startList and mystart in endList :
            returnList.append([mystart, myend, len1, len2])
            # print eachitem
    
    returnList.sort()

    return returnList

        
def blockExtraStored(storedStrand, myExtraLinkList, folderName):
    print "blockExtraStored"
    myExtraLinkList.sort()
    f = open(folderName + 'extraConnect.txt', 'w')
    for key, items in groupby(myExtraLinkList, itemgetter(0, 1)):
        mystart = key[0]
        myend = key[1]
        maxLen = -1
        storedLenPair = [-1, -1]
        count = 0
        for eachsub in items:
            print "\t", eachsub
            count += 1
            len1, len2 = eachsub[2], eachsub[3]
            if min(len1, len2) > maxLen:
                maxLen = min(len1, len2)
                storedLenPair = [len1, len2]
        print "Copy count", count, maxLen, storedLenPair
        # if count > 1:
        f.write(str(storedStrand[mystart][0]) + '_' + storedStrand[mystart][1] + ';' + str(storedStrand[myend][0]) + '_' + storedStrand[myend][1] + ';' + str(storedLenPair[0]) + ';' + str(storedLenPair[1]) + '\n')
            
    f.close()
    # print storedStrand 
    # print myExtraLinkList
        

def formMatchPairFromReadInfo(dataSet):
    dataSet.sort()
    
    matchPair = []
    for key, items in groupby(dataSet, itemgetter(0)):
        print "key", key
        left = []
        right = []
        
        for subitem in items:

            myArr = subitem[1].split('_')
            orientation = myArr[1]

            if orientation == 'p' :
                contigNum = int(myArr[0][6:]) * 2
            else:
                contigNum = int(myArr[0][6:]) * 2 + 1
            
            if subitem[2] == 'L':
                left.append([contigNum, subitem[3]])
            else:
                right.append([contigNum, subitem[3]])

        for eachleft in left:
            for eachright in right:
                leftIndex , rightIndex = eachleft[0], eachright[0]
                leftLen, rightLen = eachleft[1], eachright[1]
                
                if leftIndex != rightIndex and leftIndex / 2 != rightIndex / 2:
                    matchPair.append([rightIndex, leftIndex, min(leftLen, rightLen), rightLen, leftLen, key])

    matchPair.sort()
    return matchPair



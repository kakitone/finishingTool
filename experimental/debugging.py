import os
from operator import itemgetter
from itertools import groupby


def reverseComplement(myStr):
    myNewStr = myStr[::-1]
    myNewStr2 = ""
    for i in range(len(myNewStr)):
        if myNewStr[i] == 'A' or myNewStr[i] == 'a':
            myNewStr2 += 'T'
            
        elif  myNewStr[i] == 'T' or myNewStr[i] == 't':
            myNewStr2 +='A'
            
        elif  myNewStr[i] == 'C' or myNewStr[i] == 'c':
            myNewStr2 += 'G'
            
        elif  myNewStr[i] == 'G' or myNewStr[i] == 'g':
            myNewStr2 += 'C'
            
    return myNewStr2



def findReadsForSeekers(leftSeekerList, rightSeekerList):
    print "finding associated reads "

def writeToSeekerFile(fout, name, tmp, truncate):
    
    reverseTmp = reverseComplement(tmp)
    
    segmentName = ">"+name + "_p" + "_r"
    fout.write(segmentName)
    fout.write('\n')
    start, end = max(0, len(tmp) - truncate), len(tmp)
    fout.write(tmp[start:end])
    fout.write('\n')
    
    segmentName = ">"+name + "_p" + "_l"
    fout.write(segmentName)
    fout.write('\n')
    start, end = 0, min(truncate,len(tmp))
    fout.write(tmp[start:end])
    fout.write('\n')

    segmentName = ">"+name + "_d" + "_r"
    fout.write(segmentName)
    fout.write('\n')
    start, end = max(0, len(tmp) - truncate), len(tmp)
    fout.write(reverseTmp[start:end])
    fout.write('\n')

    segmentName = ">"+name + "_d" + "_l"
    fout.write(segmentName)
    fout.write('\n')
    start, end = 0, min(truncate,len(tmp))
    fout.write(reverseTmp[start:end])
    fout.write('\n')
    
def test2():
    print "Hello world !"
    ### TODO: successive reads extension 
    # Input: improved2.fasta
    
    # Steps to be done: 
    # a) read in the data V 
    # b) Iteratively find seekers
    # c) check set intersection and extend 
    # d) polish and output to files
    
    # Output: improved3.fasta

    # Set parameters if any
    truncate = 50000
    folderName = "illuminaContig/"
    numberOfScan = 3
    numberOfFiles = 20
    mummerLink = "MUMmer3.23/"
    thres = 400
    endThres = 10 
    
    
    # a)
    fout = open(folderName + "seekers.fasta", 'w')
    
    f = open(folderName+ "improved2.fasta", 'r')
    tmp= f.readline().rstrip()
    
    connectivityDict = {}
    name = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            name = tmp[1:]
        else:
            writeToSeekerFile(fout, name, tmp,truncate)  
            
            connectivityDict[name+'_p'+'_r' ] = []
            connectivityDict[name+'_p'+'_l' ] = []
            connectivityDict[name+'_d'+'_r' ] = []
            connectivityDict[name+'_d'+'_l' ] = []
                        
        tmp = f.readline().rstrip()
        
    f.close()
    fout.close()
    
    # b)
    
    for run in range(numberOfScan): 
        overlapList = []
        
        for dummyI in range(1, numberOfFiles+1):
            indexOfMum = ""
            if dummyI < 10:
                indexOfMum = "0"+str(dummyI)
            else:
                indexOfMum = str(dummyI)
            
            seekerLenDic = {}
            seekerLenDic = findSeekerLen(folderName)
            
            if True:
                command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"seekerout "+ folderName+ "seekers.fasta "+folderName+"raw_reads.part-"+indexOfMum+".fasta"
                os.system(command)
            
                command  = mummerLink +"show-coords -r "+folderName+"seekerout.delta > "+folderName+"fromMumseekerout"+indexOfMum 
                os.system(command)
                
            
            f = open(folderName + "fromMumseekerout"+indexOfMum, 'r')
        
            for i in range(6):
                tmp = f.readline()
            
            while len(tmp) > 0:
                infoArr = tmp.split('|')
                myArr = infoArr[-1].split('\t')
                rdGpArr = infoArr[-1].split('\t')
                
                contigName = rdGpArr[0].rstrip().lstrip()
                readName = rdGpArr[1].rstrip().lstrip()
                
                contigPosArr = infoArr[0].split(" ")
                readPosArr = infoArr[1].split(" ")
                
                overlapLenArr = infoArr[2].split()
                
                pos = []
                for eachitem in contigPosArr:
                    if len(eachitem) > 0:
                        pos.append(int(eachitem))
                        
                startPosContig = pos[0]
                endPosContig = pos[1]
                
                pos = []
                for eachitem in readPosArr:
                    if len(eachitem) > 0:
                        pos.append(int(eachitem))
                        
                startPosRead = pos[0]
                endPosRead = pos[1]
                
                pos = []
                for eachitem in overlapLenArr:
                    if len(eachitem) > 0:
                        pos.append(int(eachitem))
                        
                overlapLenContig = pos[0]
                overlapLenRead = pos[1]
                
                
                contigNameDetail, contigStrand, contigSeeker = contigName.split('_')
                
                rdInfo = readName.split('/')
                rdInfo2 = rdInfo[-1].rstrip().split('_')
                extensionLen = int(rdInfo2[1]) - int(rdInfo2[0]) - overlapLenRead
                
                
                if readName == "m121024_074656_42178_c100389662550000001523034410251204_s1_p0/1562/0_6817":
                    
                    print "contigName, readName, startPosContig, endPosContig,startPosRead, endPosRead: ", contigName, readName
                    print "startPosContig, endPosContig,startPosRead, endPosRead",  startPosContig, endPosContig,startPosRead, endPosRead
                    print "contigNameDetail, contigStrand, contigSeeker", contigNameDetail, contigStrand, contigSeeker
                    
                    print "seekerLenDic[contigName],extensionLen +overlapLenRead", seekerLenDic[contigName],extensionLen +overlapLenRead
                    print "extensionLen", extensionLen

                 
                checked = meetMatchRequirement(contigNameDetail, contigStrand, contigSeeker, startPosContig, endPosContig,startPosRead, endPosRead,seekerLenDic,extensionLen +overlapLenRead)
                

                if readName == "m121024_074656_42178_c100389662550000001523034410251204_s1_p0/1562/0_6817":
                    print "checked", checked

                
                if checked:
                    if startPosRead < endPosRead:
                        overlapList.append([contigNameDetail,contigStrand,contigSeeker,readName,overlapLenContig,extensionLen,'f'] )
                    else :
                        overlapList.append([contigNameDetail,contigStrand,contigSeeker,readName,overlapLenContig,extensionLen,'r'] )
                        
                tmp  = f.readline()
            
            f.close()
        
        overlapList.sort()
        nextRoundList = []
        
        # Bug : will write to non-assigned one... 
        for key, items in groupby(overlapList, itemgetter(slice(0,3))):
            print "key", key
            
            maxExtension = -10
            for eachitem in items:
                contigNameDetail,contigStrand,contigSeeker,readName,overlapLenContig,extensionLen,orient  = eachitem

                lookUpIndex = contigNameDetail+"_"+contigStrand+"_"+contigSeeker
                 
                readContent = [readName, orient]
                print readContent
                if not readContent in  connectivityDict[lookUpIndex] :
                    connectivityDict[lookUpIndex].append(readContent)
                    
                if extensionLen >  maxExtension:
                    maxExtension, maxRead, maxOrient =extensionLen,  readName, orient
            
            
            # important: last entry is a duplicate... to indicate the extension one.
            if maxExtension > -10:
                readContent = [maxRead, maxOrient]
                connectivityDict[lookUpIndex].append(readContent)
            
        
        writeSeeker(folderName, connectivityDict)
        
    
    # c) 
    matchingList = extendThrough(folderName, connectivityDict)
    
    # d) 
    polishAndFill(folderName, matchingList, connectivityDict)


def findSeekerLen(folderName):
    seekerLen = {}
    
    f = open(folderName + "seekers.fasta", 'r')
    tmp = f.readline().rstrip()
    name = ""
    
    while len(tmp) > 0:
        if tmp[0] == '>':
            name = tmp[1:]
        else:
            seekerLen[name] = len(tmp)
            
        tmp = f.readline().rstrip()
        
        
    f.close()
    return seekerLen


def meetMatchRequirement(contigNameDetail, contigStrand, contigSeeker, startPosContig, endPosContig,startPosRead, endPosRead,seekerLenDic,extensionLen):
    thres = 8
    checked =  False
    
    contigLen = seekerLenDic[contigNameDetail+"_"+contigStrand+"_"+contigSeeker]
    readLen = extensionLen

    if contigSeeker == 'r':
        if endPosContig > contigLen - thres and  startPosRead < thres and startPosRead < endPosRead:
            checked = True
        
        elif endPosContig > contigLen - thres and  startPosRead > readLen - thres and startPosRead >= endPosRead:
            checked = True
        
    elif contigSeeker == 'l':
        if startPosContig < thres and endPosRead > readLen - thres and startPosRead < endPosRead:
            checked = True
        elif startPosContig < thres and endPosRead < thres and startPosRead >= endPosRead:
            checked = True
                
    return checked

def writeSeeker(folderName, connectivityDict):
    #print "writeSeeker"
    print "connectivityDict: ", connectivityDict
    
    f = open(folderName + "seekers.txt",'w')
    for eachitem in connectivityDict:
        if len(connectivityDict[eachitem]) > 0:
            f.write(connectivityDict[eachitem][-1][0])
            f.write('\n')
    f.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "+folderName+"seekers.txt "+folderName+"raw_reads.fasta > "+folderName+"seekerReads.fasta"
    os.system(command)
    
    
    readDic ={}
    f = open(folderName + "seekerReads.fasta", 'r')
    tmp = f.readline().rstrip()
    name = ""
    myread = ""
    
    while len(tmp) >0:
        if tmp[0] == '>':
            readDic[name] = myread
            myread = "" 
            name = tmp[1:]
        else:
            myread = myread + tmp
            
            
        tmp = f.readline().rstrip()
    
    readDic[name] = myread
    f.close()
    
    fout = open(folderName + "seekers.fasta" , 'w')
    
    for eachitem in connectivityDict:
        if len(connectivityDict[eachitem]) > 0:
            readName = connectivityDict[eachitem][-1][0]
            orient = connectivityDict[eachitem][-1][1]
            name = eachitem
            
            fout.write(">"+name+"\n")
            if orient == 'f':
                fout.write(readDic[readName])
            else:
                fout.write(reverseComplement(readDic[readName]))
                
            fout.write('\n')
         
    fout.close()
    
    
def extendThrough(folderName, connectivityDict):
    matchingList = []
    for eachitem1 in connectivityDict:
        contigDetail1, primalDual1, leftRight1 =  eachitem1.split('_')
        for eachitem2 in connectivityDict:
            contigDetail2, primalDual2, leftRight2 =  eachitem2.split('_')
            if contigDetail1 != contigDetail2 and leftRight1 != leftRight2 :
                list1 = []
                list2 = []
                
                for eachsubitem in connectivityDict[eachitem1]:
                    list1.append(eachsubitem[0]+eachsubitem[1])
                
                for eachsubitem in connectivityDict[eachitem2]:
                    list2.append(eachsubitem[0]+eachsubitem[1])
                
                
                if len(set(list1).intersection(set(list2))) > 0:
                    matchingList.append([eachitem1, eachitem2])
                    
    print "matchingList",  matchingList 

    return matchingList 

def polishAndFill(folderName, matchingList, connectivityDict):
    print "polishAndFill"
    
    
def test1():
    folderName = "fungi/"
    
    
    f = open(folderName+ "index1.txt", 'w')
    for i in range(16, 17):
        f.write("scaffold_"+str(i)+'\n')
    f.close()
    
    f = open(folderName+ "index2.txt", 'w')
    for i in range(15, 16):
        f.write("scaffold_"+str(i)+'\n')
    f.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "+folderName+"index1.txt "+folderName+"contigs.fasta > "+folderName+"smallContigs1.fasta"
    os.system(command)
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "+folderName+"index2.txt "+folderName+"contigs.fasta > "+folderName+"smallContigs2.fasta"
    os.system(command)
    
    genomeName = "reference.fasta"
    
    os.chdir("gepard-1.30")
    os.system("./gepardcmd.sh -seq1 ../"+folderName+"smallContigs1.fasta"+" -seq2 ../"+folderName+"smallContigs2.fasta -matrix matrices/edna.mat -outfile "+"../"+folderName+"shortOverlap3.png")
    os.chdir("../")  
    

def testCases():
    testFolders = ["EcoliTestRun/", "Ph1/", "MRub/","S_cerivisea/"]
    #testFolders = [ "S_cerivisea/"]
    
    for eachfile in testFolders:
        command = "python finisherSC.py " + eachfile + " MUMmer3.23/"
        os.system(command)
        
        
test2()
#testCases()

import os
from itertools import groupby
from collections import deque
from operator import itemgetter
import sys
import time
from scipy import weave
import numpy as np


### Error correction
import numpy as np
from scipy import weave
from itertools import groupby
from operator import itemgetter

def updateVote(newAlignedSeq1, newAlignedSeq2,tempVoteTable, istart, jstart, iend, jend ):

    temp = 1
        
    Lseq = min(len(newAlignedSeq1), len(newAlignedSeq2))
    cumCount = 0 
    insCumCount = 0
    
    for j in range(Lseq): 
        if newAlignedSeq2[j] == 0 :
            # format : [number of insert, base , count ] 
            subbase = newAlignedSeq1[j] -1 
            currentIns = [insCumCount, subbase, 1 ]
            existInIns = False
            for eachins, dummyindex in zip(tempVoteTable.insCount[cumCount + jstart-1], range(len(tempVoteTable.insCount[cumCount + jstart-1]))):
                #print "eachins" ,eachins
                if eachins[0:2] == currentIns[0:2]:
                    tempVoteTable.insCount[cumCount + jstart-1][dummyindex][2] = eachins[2]+1 
                    existInIns = True
                    #print "debug",tempVoteTable.insCount[cumCount + jstart-1][dummyindex][2] 
                    
            if not existInIns :
                tempVoteTable.insCount[cumCount -1  +jstart].append(currentIns)
            
            insCumCount = insCumCount + 1 
        elif newAlignedSeq1[j] == 0:
            tempVoteTable.delCount[cumCount+ jstart] += 1 
            
        elif newAlignedSeq2[j] == newAlignedSeq1[j]:
            tempVoteTable.confirmCount[cumCount + jstart] += 1 
            
        elif newAlignedSeq2[j]  != newAlignedSeq1[j]  :
            subbase = newAlignedSeq1[j] -1 
            tempVoteTable.subCount[cumCount+jstart][subbase] += 1 
            
            
        if newAlignedSeq2[j] != 0:
            if (j< len(newAlignedSeq2) -1 and  newAlignedSeq2[j+1] != 0) or j==len(newAlignedSeq2) -1:
                tempVoteTable.confirmNoIns[cumCount + jstart] += 1 
                
            cumCount = cumCount + 1 
            insCumCount = 0 
                
class voteTable(object):
    def __init__(self, index, eachlongread):
        self.index = index    
        self.longread = eachlongread

        
        
    def initVoteTable(self):
        L = len(self.longread)
        
        self.delCount = [0 for i in range(L)]
        self.confirmCount =[0 for i in range(L)]
        # 1, 2,3 ,4 - > 0, 1, 2, 3 
        self.subCount = [[0 , 0 ,0 ,0 ] for i in range(L)]
        
        self.insCount = [ [] for i in range(L) ]
        self.confirmNoIns= [0 for i in range(L)]
        
    def polished(self):
        
        print "Before", len(self.longread)
        returnString = ""
        insTotal, delTotal = 0,0
        for i in range(len(self.longread)):
            if (self.confirmCount[i]) >= self.delCount[i]:  
                if self.longread[i] == 1:
                    returnString += 'A'
                elif self.longread[i] == 2:
                    returnString += 'C'
                elif self.longread[i] == 3:
                    returnString += 'G'
                elif self.longread[i] == 4:
                    returnString += 'T'    
            else:
                delTotal += 1
                    
            maxScore = -1
            base = -1
            self.insCount[i].sort()
            
            for key, items in groupby(self.insCount[i], itemgetter(0)):
                maxScore = -1
                base = -1
                if key == 0:
                    for eachitem in items:
                        if eachitem[2] > maxScore:
                            maxScore = eachitem[2]
                            base = eachitem[1]
            
            if len(self.insCount[i]) > 0 and (self.confirmNoIns[i]) < maxScore:
                #print "self.confirmNoIns[i], maxScore, self.insCount[i]", self.confirmNoIns[i], maxScore, self.insCount[i]
                if base == 0:
                    returnString += 'A'
                elif base == 1:
                    returnString += 'C'
                elif base == 2:
                    returnString += 'G'
                elif base == 3:
                    returnString += 'T'

                insTotal += 1
        print "After: length, insTotal, delTotal",  len(returnString), insTotal, delTotal
        return returnString
        
 


def SWAlignment(seq1 , seq2):
    score  = 0 

    wts = -10
    wti = -1 
    wtd = -1
    wtm=  1
    
    m = len(seq1) + 1
    n = len(seq2) +1
    
    H = np.zeros([m,n], dtype = np.float64)
    B = np.zeros([m,n], dtype = np.float64)

    # Assign weights 
    for i in range(m):
        H[i][0] = 0
        B[i][0] = 4
    for j in range(n):
        H[0][j]  =0
        B[0][j] = 4 
        
        
        
    seq1NP = np.zeros(m-1, dtype = np.float64)
    seq2NP = np.zeros(n-1, dtype = np.float64)
    
    for i in range(m-1):
        seq1NP[i] = seq1[i]
    for j in range(n-1) :
        seq2NP[j] = seq2[j]
        
    
    code =\
        """
        int i; 
        int j ;
        double w;

        for (i =1 ;i <m ; i++){
            for (j=1; j<n; j++){
                if (seq1NP[i-1] == seq2NP[j-1]){
                    w = wtm ;
                }
                else{
                    w= wts ;
                }
                
                    H2(i,j) = 0 ;
                    
                    if (H2(i,j) < H2(i-1,j-1) + w) {
                        H2(i,j) = H2(i-1,j-1) + w ; 
                    }

                    if (H2(i,j) < H2(i-1,j)+wtd ) {
                        H2(i,j) = H2(i-1,j)+wtd ;
                    }

                    if  (H2(i,j) < H2(i,j-1) + wti){
                        H2(i,j) = H2(i,j-1) + wti;
                    }

                    if (H2(i-1,j-1) + w == H2(i,j)){
                        B2(i,j) = 1 ; 
                    }
                    else if (H2(i-1,j)+wtd == H2(i,j)){
                        B2(i,j) = 2;
                    }
                    
                    else if  (H2(i,j-1) + wti == H2(i,j)){
                        B2(i,j) = 3;
                    }
                    else if (0 == H2(i,j)) {
                        B2(i,j) = 4 ;
                    }
            }
        }

        """
    
    weave.inline(code, ['H','B','m', 'n', 'wtm','wts','wti','wtd', 'seq1NP', 'seq2NP'])
    

    # Backtrack 
    alignedSeq1 = []
    alignedSeq2 = []
    
    bestindex = np.argmax(H)

 
    
    endi = bestindex / n
    endj = bestindex%n
    
    score = H[endi][endj]
    
    tempi, tempj = endi , endj
    while (B[tempi][tempj] != 4):
        if B[tempi][tempj] == 1:
            alignedSeq1.append(seq1[tempi -1 ])
            alignedSeq2.append(seq2[tempj-1 ]) 
            tempi = tempi -1 
            tempj = tempj -1           

        elif B[tempi][tempj] == 2 :
            alignedSeq1.append(seq1[tempi-1 ])
            alignedSeq2.append(0) 
            tempi = tempi - 1

        elif B[tempi][tempj] == 3 :
            alignedSeq1.append(0)
            alignedSeq2.append(seq2[tempj-1 ])
            tempj = tempj -1

    
    starti, startj = tempi , tempj 
    
    returnalignedSeq1 = reverseString(alignedSeq1)
    returnalignedSeq2 = reverseString(alignedSeq2)
    
    return score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj


def reverseString(str):
    newstr = []
    for index in range(len(str)):
        newstr.append(str[len(str)-1-index])
    
    return newstr

def transformToDesiredForm(alignedSeq1, alignedSeq2):
    newAlignedSeq1, newAlignedSeq2 = [] ,[]

    Lseq = len(alignedSeq1)
    
    for i in range(Lseq):
        newAlignedSeq1.append(alignedSeq1[i])
        newAlignedSeq2.append(alignedSeq2[i])
    
    i =0 
    while (i <= len(newAlignedSeq1) -2 ):
        if newAlignedSeq1[i] == 0 and newAlignedSeq2[i] == 0 :
            newAlignedSeq1.pop(i)
            newAlignedSeq2.pop(i)
                    
        elif newAlignedSeq1[i] == newAlignedSeq2[i+1] and  newAlignedSeq1[i] == newAlignedSeq2[i+1] and  newAlignedSeq2[i] == 0:
            newAlignedSeq2[i] = newAlignedSeq1[i] 
            newAlignedSeq2[i+1] = 0
            i = i +1 
        elif newAlignedSeq2[i] == newAlignedSeq1[i+1] and  newAlignedSeq2[i] == newAlignedSeq1[i+1] and  newAlignedSeq1[i] == 0:
            newAlignedSeq1[i] = newAlignedSeq2[i] 
            newAlignedSeq1[i+1] = 0
            i = i +1 

        elif  newAlignedSeq1[i] == 0 and newAlignedSeq2[i+1] == 0:
            newAlignedSeq1[i] = newAlignedSeq1[i+1]
            newAlignedSeq1[i+1] =0 
            i = i +1 
        elif  newAlignedSeq2[i] == 0 and newAlignedSeq1[i+1] == 0:
            newAlignedSeq2[i] = newAlignedSeq2[i+1] 
            newAlignedSeq2[i+1] = 0   
            i = i +1          
        else: 
            i = i +1 

    return newAlignedSeq1, newAlignedSeq2



def polishing(eachlongread, helperList):
    
    tempVoteTable = voteTable(0, eachlongread) 
    tempVoteTable.initVoteTable()            
    
    for eachshortreaad in helperList:
        score, alignedSeq1, alignedSeq2, istart, jstart, iend, jend = SWAlignment(eachshortreaad, eachlongread )
        newAlignedSeq1, newAlignedSeq2 = transformToDesiredForm(alignedSeq1, alignedSeq2)
        updateVote(newAlignedSeq1, newAlignedSeq2,tempVoteTable,istart, jstart, iend, jend )
            
    cleanedLongRead = tempVoteTable.polished()
    
    return cleanedLongRead







##### House keeping files
def writeToFile_Double1(folderName, fileName1, fileName2, option = "contig"):

    f2 = open(folderName + fileName2, 'w')
    fOriginal = open(folderName + fileName1, 'r')
    
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) >0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead+ tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)
    
    print "len(readSet)", len(readSet)
    
    fOriginal.close()
    
    if option == "contig":
        header = ">Contig"
    else:
        header = ">Read"
    for eachcontig, dum in zip(readSet, range(len(readSet))):
        f2.write(header+ str(dum)+"_p\n")
        f2.write(eachcontig+'\n') 
        f2.write(header+ str(dum)+"_d\n")
        f2.write(reverseComplement(eachcontig)+'\n')
        
    f2.close()
                   
   
def writeToFile_Double2(folderName, fileName1, fileName2):
    # for reads
    
    f = open(folderName + fileName1, 'r')
    f2 = open(folderName + fileName2, 'w')
    
    tmp1 = f.readline().rstrip()
    if len(tmp1) > 0:
        tmp2 = f.readline().rstrip()
        
    while len(tmp1) > 0 and len(tmp2) >0:
        f2.write(tmp1+ "_p")
        f2.write('\n')
        f2.write(tmp2)
        f2.write('\n')
        f2.write(tmp1 + "_d")
        f2.write('\n')
        f2.write(reverseComplement(tmp2))
        f2.write('\n')
        
        tmp1 = f.readline().rstrip()
        
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
    
    f2.close()
    f.close()

def writeToFile(f2, runningIndex,seq):
    f2.write(">Seg_"+ str(runningIndex))
    f2.write('\n')
    f2.write(seq)
    f2.write('\n')
   

def findContigLength(folderName, option):
    
    if option == "contigs":
        print "\n\nfindContigLength(folderName)\n"
        contigLength = {}
        f = open(folderName + "smaller_contigs_Double.fasta", 'r')
        
        
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) >0:
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()
        
        
    else:
        print "\n\nfindContigLength(folderName)\n"
        contigLength = {}
        f = open(folderName + "improved_Double.fasta", 'r')
        
        
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) >0:
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()
        
    
    
        f = open(folderName + "relatedReads_Double.fasta", 'r')
        
    
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) >0:
    
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()

    return contigLength

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

##### New method of using string graph to approximate DB graph

    
def formRelatedReadsFile(folderName,mummerLink):    
    # Find associated read and extract into a file associatedReads.fasta
    # Input: contigs.fasta, cleaned_Reads.fasta 
    # Output: relatedReads.fasta

    ### Extract heads of the contigs
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
                writeToFile(f2, runningIndex,tempContig[0:thres])
                runningIndex = runningIndex +1
                
                writeToFile(f2, runningIndex,tempContig[-thres:])
                runningIndex = runningIndex +1 
                
                                
                writeToFile(f2, runningIndex,reverseComplement(tempContig[0:thres]))
                runningIndex = runningIndex +1
                
                writeToFile(f2, runningIndex,reverseComplement(tempContig[-thres:]))
                runningIndex = runningIndex +1
                
                tempContig = ""
        else:
            tempContig = tempContig + temp
        
        temp = f.readline()

    writeToFile(f2, runningIndex,tempContig[0:thres])
    runningIndex = runningIndex +1
    
    writeToFile(f2, runningIndex,tempContig[-thres:])
    runningIndex = runningIndex +1
                  
    writeToFile(f2, runningIndex,reverseComplement(tempContig[0:thres]))
    runningIndex = runningIndex +1
    
    writeToFile(f2, runningIndex,reverseComplement(tempContig[-thres:]))
    runningIndex = runningIndex +1
    
    
    f2.close()
    f.close()
    
    ### Write double stranded reads
    writeToFile_Double1(folderName, "improved.fasta", "improved_Double.fasta", "contig")
    #writeToFile_Double1(folderName, "raw_reads.fasta", "raw_reads_Double.fasta","read")
    
    ### Apply MUMMER on them using cleanedReads against them
    assoiatedReadIndex = []
    nameList = []
    numberOfFiles = 20
    bindir = os.path.abspath(os.path.dirname(sys.argv[0]))
    command = bindir+"/fasta-splitter.pl --n-parts "+str(numberOfFiles)+" "+ folderName+"raw_reads.fasta"
    os.system(command)
    
    for dummyI in range(1, numberOfFiles+1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0"+str(dummyI)
        else:
            indexOfMum = str(dummyI)
        
        if True:
            command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"out "+ folderName+ "improvedTrunc.fasta raw_reads.part-"+indexOfMum+".fasta"
            os.system(command)
    
            command  = mummerLink +"show-coords -r "+folderName+"out.delta > "+folderName+"fromMum"+indexOfMum
            os.system(command)
        
        f = open(folderName + "fromMum"+indexOfMum, 'r')
    
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
            if startPos < endThres and endPos > thres-endThres:
                assoiatedReadIndex.append(myArr[1])
                nameList.append([int(contigName.split('_')[1]), readName])
            tmp  = f.readline()
        
        f.close()
    
    
    nameList.sort()

    assoiatedReadIndex.sort()
    
    print "assoiatedReadIndex", assoiatedReadIndex
    
    ckIndex = 0
    f = open(folderName+"associatedNames.txt", 'w')
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
                f.write(key+'\n')
                keyFound.append(key)

        ckIndex += 1
    
    print "ckIndex,oneItem: ",ckIndex, oneItem
    f.close()

    fFilter = open(folderName+"associatedNames.txt", 'r')
    
    fout = open(folderName+"associatedNames2.txt",'w') 
    
    maxCount = 12000
    mytmpDum = fFilter.readline() 
    i=0
    while i < maxCount and len(mytmpDum)> 0:
        fout.write(mytmpDum)  
        mytmpDum = fFilter.readline() 
        i = i+1
        
    fout.close()   
    fFilter.close()

    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "+folderName+"associatedNames2.txt "+folderName+"raw_reads.fasta > "+folderName+"relatedReads.fasta"
    os.system(command)
    
    writeToFile_Double1(folderName, "relatedReads.fasta", "relatedReads_Double.fasta","read")
    
    numberOfFiles = 20
    bindir = os.path.abspath(os.path.dirname(sys.argv[0]))
    command = bindir+"/fasta-splitter.pl --n-parts "+str(numberOfFiles)+" "+ folderName+"relatedReads_Double.fasta"
    os.system(command)
    

def extractEdgeSet(folderName, mummerLink, option= "nopolish"):
    # Tasks: reconstruct the string  graph
    
    # Input : relatedReads_Double.fasta, conig_Double.fasta
    # Intermediate files: fromMum_overlap , fromMum_overlap
    # Output: connectivity of eachNode: InList, OutList [critical]
    #         connectivity of eachNode: arrow representatiowith size [optional]
    
    
    ### Perform MUMMER alignment
    print ">Extract Edge set"
    lengthDic = findContigLength(folderName, "improved")
    numberOfContig = 0
    f = open(folderName+ "improved_Double.fasta",'r')
    tmp = f.readline()
    tmp = tmp.rstrip()
    while (len(tmp)>0):
        numberOfContig += 1
        tmp = f.readline()
        tmp = tmp.rstrip()
    
    numberOfContig = numberOfContig /2
    print "numberOfContig", numberOfContig

    f.close()
    K = 400
    
    thres = 5
    # Nodes are contigs, 
    nodes = [i for i in range(numberOfContig)]
    dataSet = []
    
    ### Apply MUMMER on them using cleanedReads against them
    fmyFile = open(folderName+ "improved_Double.fasta", 'r')
    fSmaller = open(folderName+ "smaller_improvedContig.fasta", 'w')

    tmp = fmyFile.readline().rstrip()
    maxSize = 25000
    
    dummySeq = ""
    for i in range(0):
        dummySeq = dummySeq + "A"
        
    
    myName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            fSmaller.write(tmp+'\n')
            myName = tmp[1:]
        else:
            component = tmp[0:min(len(tmp), maxSize )] 
            countComp = len(component)
            fSmaller.write(component)
            
            fSmaller.write(dummySeq)
            countComp = countComp + len(dummySeq)
            
            component =tmp[max(0, len(tmp)-maxSize):len(tmp)]
            fSmaller.write(component)
            countComp = countComp + len(component)
            
            lengthDic[myName] = countComp 
            print "DebugName", myName, countComp
            fSmaller.write('\n')

        tmp = fmyFile.readline().rstrip()

    fSmaller.close()
    fmyFile.close()
    
    numberOfFiles = 20
    for dummyI in range(1, numberOfFiles+1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0"+str(dummyI)
        else:
            indexOfMum = str(dummyI)

        if True:
            command =mummerLink +"nucmer --maxmatch --simplify -p "+folderName+"outRefine "+ folderName+ "smaller_improvedContig.fasta "+ "relatedReads_Double.part-"+indexOfMum+".fasta"
            os.system(command)
            
            command  = mummerLink +"show-coords -r "+folderName+"outRefine.delta > "+folderName+"fromMumRefine"+indexOfMum
            os.system(command)
            
        f = open(folderName + "fromMumRefine"+indexOfMum, 'r')
        for i in range(6):
            tmp = f.readline()
    
        while len(tmp) > 0:
            info = tmp.split('|')
            filterArr =  info[1].split()
            rdGpArr = info[-1].split('\t')
            firstArr = info[0].split()
            matchLenArr = info[2].split()
           
        
            matchLen = int(matchLenArr[1])    
            contigStart, contigEnd =  int( firstArr[0]), int( firstArr[1])
            readStart, readEnd =  int(filterArr[0]) , int(filterArr[1])
            
            contigName = rdGpArr[0].rstrip().lstrip()
            readName = rdGpArr[1].rstrip().lstrip()
                
        
            if readStart < readEnd and matchLen> K and min(contigStart,readStart)  < thres and min(lengthDic[contigName]- contigEnd,  lengthDic[readName] - readEnd) < thres:
                conditionForMatch = True
            else:
                conditionForMatch = False
    
            if conditionForMatch :
                if contigStart < thres:
                    dataSet.append((readName, contigName, 'L',matchLen))
                
                if lengthDic[contigName]- contigEnd < thres :
                    dataSet.append((readName, contigName, 'R',matchLen))
            
            tmp = f.readline()
            
        f.close()  
        
        
    dataSet.sort()
    
    print "dataSet:", dataSet[0]

    matchPair = []
    for key, items in groupby(dataSet, itemgetter(0)):
        left = []
        right = []
        
        for subitem in items:

            myArr =subitem[1].split('_')
            orientation = myArr[1]

            if orientation == 'p' :
                contigNum = int(myArr[0][6:])*2
            else:
                contigNum = int(myArr[0][6:])*2 +1
            
            if subitem[2] == 'L':
                left.append([contigNum, subitem[3]])
            else:
                right.append([contigNum, subitem[3]])

        for eachleft in left:
            for eachright in right:
                leftIndex , rightIndex = eachleft[0], eachright[0]
                leftLen, rightLen = eachleft[1], eachright[1]
                
                if leftIndex != rightIndex:
                    matchPair.append([rightIndex, leftIndex, min(leftLen,rightLen),  rightLen,leftLen, key])

    matchPair.sort()
        
    keyFound = []
    bestMatchPair = []
    rawReadList = []
    
    for key, items in groupby(matchPair, itemgetter(0,1)):
        maxvalue = -1
        maxLenPair = []
        for eachitem in items:
            if eachitem[2] > maxvalue:
                maxvalue = eachitem[2]
                maxLenPair = [eachitem[3], eachitem[4], eachitem[5]]
        bestMatchPair.append([key[0], key[1], maxvalue, maxLenPair[0], maxLenPair[1], maxLenPair[2]])
    
    bestMatchPair.sort( key=lambda tup: tup[2], reverse=True)
    
    print "bestMatchPair", bestMatchPair


    leftConnect = [[-1,-1,-1] for i in range(numberOfContig)]
    rightConnect = [[-1,-1,-1] for i in range(numberOfContig)]
    
    for eachitem in bestMatchPair:

        prefixContig = eachitem[0]
        suffixContig = eachitem[1]
        
        if leftConnect[suffixContig][0] == -1 and rightConnect[prefixContig][0] == -1:
            leftConnect[suffixContig][0] = prefixContig 
            leftConnect[suffixContig][1] =  eachitem[4]
            leftConnect[suffixContig][2] =  eachitem[5]
            
            rightConnect[prefixContig][0] = suffixContig
            rightConnect[prefixContig][1] = eachitem[3]
            rightConnect[prefixContig][2] = eachitem[5]
            
            rawReadList.append(eachitem[5])
    
    startList= []
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
        while rightConnect[tmp][0] != -1 and rightConnect[tmp][0]!=mystart:
            tmp = rightConnect[tmp][0]
            if tmp != -1:
                myList.append(tmp)
        contigList.append(myList)
    
    print "contigList", contigList

            
    i = 0
    fOriginal = open(folderName + "improved.fasta", 'r')
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) >0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead+ tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)  
    fOriginal.close()
    
    ### Put the needed rawReads into the RAM using Dictionary
    
    fAppendRaw = open(folderName+ "appendRaw.txt", 'w')
    for eachraw in rawReadList:
        fAppendRaw.write(eachraw)
        fAppendRaw.write('\n')
    fAppendRaw.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "+folderName+"appendRaw.txt "+folderName+"relatedReads_Double.fasta > "+folderName+"rawToAppend.fasta"
    os.system(command)
    
    rawRead = {}
    
    fOriginal = open(folderName + "rawToAppend.fasta", 'r')
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    tmpName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            
            if len(tmpRead) >0:
                rawRead[tmpName]=tmpRead
                tmpRead = ""
                
            tmpName = tmp[1:]
        else:
            tmpRead = tmpRead+ tmp
            
        tmp = fOriginal.readline().rstrip()
        
    rawRead[tmpName]=tmpRead
    ### End
    
    seqToPrint = []
    contigUsed = [False for i in range(numberOfContig/2)]
    fAllInOne = open(folderName+"allInOne2.fasta",'w' )
    finalList = []
    for eachContig, i in zip(contigList, range(len(contigList))):
        tmpList = []
        for eachitem in eachContig:

            readNum = eachitem/2
            if contigUsed[readNum] ==False:
                seqToPrint.append(eachitem)
                tmpList.append(eachitem)
                contigUsed[readNum] = True
        if len(tmpList) > 0:
            finalList.append(tmpList)
        
    
    fImproved = open(folderName +"improved2.fasta", 'w')
    
    for eachcontig, dummyIndex in zip(finalList, range(len(finalList))):
        fImproved.write(">Segkk"+str(dummyIndex)+'\n')
        tmpStore = -1997
        tmpStore2 = -1998
        tmpStore3 = -1999
        
        for eachseg, hidum in zip(eachcontig, range(len(eachcontig))):
            readNum = eachseg/2
            orientation = eachseg%2
            
            newStart = 0 
            
            x , y , l= tmpStore,leftConnect[eachseg][1],tmpStore2
            extraRead = ""
            if hidum ==0:
                newStart = 0
            else:
                
                if l <x+y:
                    newStart = x+y-l
                else:
                    newStart = 0
                    if option == 'polish':
                        extraRead = performPolishing(leftConnect[eachseg][0], eachseg, tmpStore3[x:l-y],  dataSet, folderName)
                    else:
                        extraRead = tmpStore3[x:l-y]
                    
    
            print extraRead[0:10], len(extraRead)
            
            
            fImproved.write(extraRead)
            
            if orientation == 0:
                fImproved.write(readSet[readNum][newStart:])   

            else:
                fImproved.write(reverseComplement(readSet[readNum])[newStart:])

            
            if rightConnect[eachseg][1] != -1:
                tmpStore = rightConnect[eachseg][1]
                tmpStore2 = len(rawRead[rightConnect[eachseg][2]])
                tmpStore3 = rawRead[rightConnect[eachseg][2]]
                
        fImproved.write('\n')
        
    fImproved.close()
            


    for eachContigIndex,dum in zip(seqToPrint,range(len(seqToPrint))):

        readNum = eachContigIndex/2
        orientation = eachContigIndex%2
        
        fAllInOne.write(">Seg_"+str(dum)+'\n')
        if orientation == 0:
            fAllInOne.write(readSet[readNum])
            fAllInOne.write('\n')
        else:
            fAllInOne.write(reverseComplement(readSet[readNum]))
            fAllInOne.write('\n')
                    
    fAllInOne.close()
    


    for eachcontig, i in zip(finalList, range(len(finalList))):
        fout = open(folderName + "improvedContig2_"+str(i)+".fasta", 'w')
        seqToPrint = []

        for eachitem in eachcontig:

            readNum = eachitem/2
            seqToPrint.append(eachitem)

        print "ImprovedContig ",i 
        for eachhaha in seqToPrint:
            print len(readSet[eachhaha/2])
            
        for eachContigIndex,dum in zip(seqToPrint,range(len(seqToPrint))):

            readNum = eachContigIndex/2
            orientation = eachContigIndex%2
            
            fout.write(">Seg_"+str(dum)+'\n')
            if orientation == 0:
                fout.write(readSet[readNum])
                fout.write('\n')
            else:
                fout.write(reverseComplement(readSet[readNum]))
                fout.write('\n')


        i += 1
        fout.close()




def performPolishing(contigid_R, contigid_L, mySegment,  dataSet, folderName):
    ## R,L refer to R/L seekers
    # dataSet: ('Read100_d', 'Contig0_p', 'R', 684)
    
    # Find reads covering one flanking region
    helperReadsList = []
    readNum = contigid_R/2
    orientation = contigid_R%2
    contigName = ""
    if orientation == 0:
        contigName = "Contig"+ str(readNum) + "_p"
    else:
        contigName = "Contig"+ str(readNum) + "_d"
    
    for key, items in groupby(dataSet, itemgetter(1,2)):
        if key[0] == contigName and key[1] == 'R':
            for eachsubitem in items:
                helperReadsList.append(eachsubitem[0])
    
    readNum = contigid_L/2
    orientation = contigid_L%2
    contigName = ""
    if orientation == 0:
        contigName = "Contig"+ str(readNum) + "_p"
    else:
        contigName = "Contig"+ str(readNum) + "_d"
    
    for key, items in groupby(dataSet, itemgetter(1,2)):
        if key[0] == contigName and key[1] == 'L':
            for eachsubitem in items:
                helperReadsList.append(eachsubitem[0])    

    f = open(folderName + "helperList.txt", 'w')
    for eachitem in helperReadsList:
        f.write(eachitem)
        f.write('\n')
    f.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "+folderName+"helperList.txt "+folderName+"relatedReads_Double.fasta > "+folderName+"SR.fasta"
    os.system(command)

    f = open(folderName +"LR.fasta", 'w')
    f.write(">Seg\n")
    f.write(mySegment)
    f.close()
    
    # Perform error correction using outside tool 
    helperList = []
    command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"helperOut "+ folderName+ "SR.fasta "+ folderName+ "LR.fasta"
    os.system(command)
    
    command  = mummerLink +"show-coords -r "+folderName+"helperOut.delta > "+folderName+"fromhelperOut"
    os.system(command)
    
    f = open(folderName+"fromhelperOut",'r')
    helperNameList = []
    
    for i in range(6):
        tmp = f.readline()

    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr =  info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        matchLenArr = info[2].split()
    
        matchLen = int(matchLenArr[1])    
        helperStart, helperEnd =  int( firstArr[0]), int( firstArr[1])
        readStart, readEnd =  int(filterArr[0]) , int(filterArr[1])
        
        helperName = rdGpArr[0].rstrip().lstrip()
        readName = rdGpArr[1].rstrip().lstrip()
        
        helperNameList.append(helperName)
        
        tmp = f.readline().rstrip()
                
    f.close()

    
    f= open(folderName + "SR.fasta", 'r')
    temp = f.readline().rstrip()
    name = ""
    
    while(len(temp) > 0):
        if temp[0] == '>':
            name = temp[1:]
        else:
            if name in helperNameList:
                helperList.append(temp)
                
        temp = f.readline().rstrip()
    
    
    
    f.close()
    
    mySegmentInt, helperIntList= [],[]

    mySegmentInt = covertAlphToInt(mySegment)
    for eachitem in helperList:
        helperIntList.append(covertAlphToInt(eachitem))
    
        
    # Return the associated read
    if len(helperIntList) > 0:
        polishedSegment = polishing(mySegmentInt,helperIntList )
    else:
        polishedSegment = mySegment
    
    
    return polishedSegment

def covertAlphToInt(mySegmentIn):
    mySegment = []
    for eachitem in mySegmentIn:
        if eachitem == 'A':
            mySegment.append(1)
        elif eachitem == 'C':
            mySegment.append(2)
        elif eachitem == 'G':
            mySegment.append(3)
        elif eachitem == 'T':
            mySegment.append(4)
    return mySegment

def newGraphPipeLine(folderName, mummerLink,option):
    print "newGraphPipeLine"
    formRelatedReadsFile(folderName,mummerLink)
    extractEdgeSet(folderName, mummerLink, option)
    
    
def greedyAlg(mummerLink, folderName):
    
    print "Direct greedy"
    
    thres = 7
    minLen = 400
    

    writeToFile_Double1(folderName, "contigs.fasta", "contigs_Double.fasta", "contig")
    
    fmyFile = open(folderName+ "contigs_Double.fasta", 'r')
    fSmaller = open(folderName+ "smaller_contigs_Double.fasta", 'w')

    tmp = fmyFile.readline().rstrip()
    maxSize = 50000

    myName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            fSmaller.write(tmp+'\n')
            myName = tmp[1:]
        else:
            component = tmp[0:min(len(tmp), maxSize )] 
            countComp = len(component)
            fSmaller.write(component)
            
            component =tmp[max(0, len(tmp)-maxSize):len(tmp)]
            fSmaller.write(component)
            countComp = countComp + len(component)
            

            print "DebugName", myName, countComp
            fSmaller.write('\n')

        tmp = fmyFile.readline().rstrip()

    fSmaller.close()
    fmyFile.close()
    
    if True:
        command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"greedy "+ folderName+ "smaller_contigs_Double.fasta "+ folderName+ "smaller_contigs_Double.fasta"
        os.system(command)
        
        command  = mummerLink +"show-coords -r "+folderName+"greedy.delta > "+folderName+"fromMumGreedy"
        os.system(command)
        
    lengthDic = findContigLength(folderName, "contigs")
    

    f = open(folderName + "fromMumGreedy", 'r')
    for i in range(6):
        tmp = f.readline()
        
    dataSet = []
    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr =  info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        matchLenArr = info[2].split()
       
    
        matchLen = int(matchLenArr[0])    
        contigStart, contigEnd =  int( firstArr[0]), int( firstArr[1])
        contig2Start, contig2End =  int(filterArr[0]) , int(filterArr[1])
        
        contigName = rdGpArr[0].rstrip().lstrip()
        contig2Name = rdGpArr[1].rstrip().lstrip()
            
    
        if contigName != contig2Name and matchLen > minLen and contig2Start < contig2End  and min(contigStart,contig2Start) < thres and min(lengthDic[contigName]- contigEnd,  lengthDic[contig2Name] - contig2End) < thres:
            conditionForMatch = True
        else:
            conditionForMatch = False

        if conditionForMatch :
            if contigStart < thres:
                dataSet.append((matchLen, contig2Name, contigName))
            
        
        tmp = f.readline()
        
    f.close()  
    dataSet.sort(reverse=True)
    
    print "dataSet:", dataSet[0]
    numberOfContig = len(lengthDic)
    
    # [next_item, overlap_length]
    
    leftConnect = [[-1,-1] for i in range(numberOfContig)]
    rightConnect = [[-1,-1] for i in range(numberOfContig)]
    
    for eachitem in dataSet:
        prefix = eachitem[1].split('_')
        suffix = eachitem[2].split('_')
        lengthOfOverlap = int(eachitem[0])
        
        if prefix[1] == 'p':
            prefixContig = int(prefix[0][6:])*2 
        else:
            prefixContig = int(prefix[0][6:])*2 +1
        
        if suffix[1] == 'p':
            suffixContig = int(suffix[0][6:])*2 
        else:
            suffixContig = int(suffix[0][6:])*2 +1
            
        
        if leftConnect[suffixContig][0] == -1 and rightConnect[prefixContig][0] == -1:
            leftConnect[suffixContig][0] = prefixContig 
            leftConnect[suffixContig][1] = lengthOfOverlap
            rightConnect[prefixContig][0] = suffixContig
            rightConnect[prefixContig][1] = lengthOfOverlap
            
    startList= []
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
        while rightConnect[tmp][0] != -1 and rightConnect[tmp][0]!=mystart:
            tmp = rightConnect[tmp][0]
            if tmp != -1:
                myList.append(tmp)
        contigList.append(myList)
            
    i = 0
    fOriginal = open(folderName + "contigs.fasta", 'r')
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) >0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead+ tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)  
    
    seqToPrint = []
    contigUsed = [False for i in range(numberOfContig/2)]
    fAllInOne = open(folderName+"allInOne.fasta",'w' )
    finalList = []
    for eachContig, i in zip(contigList, range(len(contigList))):
        tmpList = []
        for eachitem in eachContig:

            readNum = eachitem/2
            if contigUsed[readNum] ==False:
                seqToPrint.append(eachitem)
                tmpList.append(eachitem)
                contigUsed[readNum] = True
        if len(tmpList) > 0:
            finalList.append(tmpList)
    
    fImproved = open(folderName +"improved.fasta", 'w')
    for eachcontig, dummyIndex in zip(finalList, range(len(finalList))):
        fImproved.write(">Segkk"+str(dummyIndex)+'\n')
        for eachseg in eachcontig:
            readNum = eachseg/2
            orientation = eachseg%2
            overlapLength =  rightConnect[eachseg][1]
            if overlapLength == -1:
                endPt = len(readSet[readNum])
            else:
                endPt = len(readSet[readNum]) - overlapLength
                
            if orientation == 0:
                fImproved.write(readSet[readNum][0:endPt])
            else:
                fImproved.write(reverseComplement(readSet[readNum])[0:endPt])
        fImproved.write('\n')
        
    fImproved.close()
    
def mainFlow(folderName,mummerLink,option):
    greedyAlg(mummerLink, folderName)
    newGraphPipeLine(folderName, mummerLink, option)
    
    
    #os.system("cp raw_reads.part* "+ folderName)
    os.system("rm raw_reads.part*")
    
    #os.system("cp relatedReads_Double.part* "+ folderName)
    os.system("rm relatedReads_Double.part*")
    
    
t0 = time.time()
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

folderName = sys.argv[1]
mummerLink = sys.argv[2]

option = "nopolish"
if len(sys.argv) > 3:
    if sys.argv[3] == "polish" :
        option =  "polish"

mainFlow(folderName,mummerLink, option)
print  "Time",  time.time() - t0

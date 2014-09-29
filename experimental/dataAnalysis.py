import os
from itertools import groupby
from operator import itemgetter
import matplotlib.pyplot as plt

### TODO
### Function to anaylze the data : 
# 1) Mummer the data and get [contigNum, startIndex,endIndex] by choosing the best
# 2) Remove and report any completely embedded segments
# 3) Compute successor overlaps and report NGaps and overlap histogram

def extractMumData(folderName, fileName):
    #"Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    f = open(folderName+fileName,'r')
    dataList = []
    
    for i in range(6):
        tmp = f.readline()

    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr =  info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        
        matchLenArr = info[2].split()
    
        matchLen1 = int(matchLenArr[0])
        matchLen2 = int(matchLenArr[1])    
        percentMatch = float(info[3])
        
        
        helperStart, helperEnd =  int( firstArr[0]), int( firstArr[1])
        readStart, readEnd =  int(filterArr[0]) , int(filterArr[1])
        
        helperName = rdGpArr[0].rstrip().lstrip()
        readName = rdGpArr[1].rstrip().lstrip()
        
        dataList.append([helperStart, helperEnd , readStart, readEnd,matchLen1,matchLen2,percentMatch,helperName,readName ])
    
        tmp = f.readline().rstrip()
                
    f.close()
    
    return dataList


def obtainLength(folderName, fileName):
    f = open(folderName+fileName, 'r')
    tmp = f.readline().rstrip()
    lenDic = {}
    tmplen = 0 
    tmpName = ""
    
    while len(tmp) > 0:
        
        if tmp[0] == '>':
            if tmplen != 0:
                lenDic[tmpName] =tmplen
                tmplen = 0
            tmpName = tmp[1:]
        else:
            tmplen += len(tmp)
        tmp  = f.readline().rstrip()
        
    lenDic[tmpName] =tmplen
    

    f.close()
    
    return lenDic
    
    
def putListToFileO(folderName, sourceFileName, targetFileName, myList):
    f = open(folderName + targetFileName+".txt", 'w')
    for eachitem in myList:
        f.write(eachitem+'\n')
        
    f.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "+folderName+targetFileName+".txt "+folderName+sourceFileName+" > "+folderName+targetFileName+".fasta"
    os.system(command)
    
def transformCoor(dataList):
    #"Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
 
    newList = []
    
    for eachitem in dataList:
        if eachitem[2] < eachitem[3]:
            newList.append(eachitem)
        else:
            tmpitem = eachitem
            tmp = tmpitem[2]
            tmpitem[2] = tmpitem[3]
            tmpitem[3] = tmp
            newList.append(tmpitem)
    
    return newList
def analyseOverlap(folderName = "S_cerivisea/"):
    # "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
   
    print "Hello world"
    
    mummerLink = "MUMmer3.23/"
    
    # 0) Mummer self align
    if False:
        command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"dout "+ folderName+ "reference.fasta "+ folderName+ "reference.fasta"
        os.system(command)
        
        command  = mummerLink +"show-coords -r "+folderName+"dout.delta > "+folderName+"doutOut"
        os.system(command)
    
    if False:
        command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"ralignment "+ folderName+ "repeat.fasta "+ folderName+ "repeat.fasta"
        os.system(command)
        
        command  = mummerLink +"show-coords -r "+folderName+"ralignment.delta > "+folderName+"ralignmentOut"
        os.system(command)        
    
    
    # 1) Mummer the data and get [referenceName, contigName ,  startIndex,endIndex, lengthMatch] by choosing the best
    if True:
        command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"alignment "+ folderName+ "reference.fasta "+ folderName+ "repeat.fasta"
        os.system(command)
        
        command  = mummerLink +"show-coords -r "+folderName+"alignment.delta > "+folderName+"alignmentOut"
        os.system(command)
        
    dataList = extractMumData(folderName, "alignmentOut")
    #print dataList[0]
    
    dataList.sort(key = itemgetter(-1))
    
    alignmentList = []
    for key, items in groupby(dataList, itemgetter(-1)):
        maxMatch = -1
        tmpitem = []
        for eachitem in items:
            if eachitem[4] > maxMatch:
                tmpitem = [ eachitem[-2], eachitem[-1], eachitem[0], eachitem[1] , eachitem[4]]
                maxMatch = eachitem[4]
                
        alignmentList.append(tmpitem)
            

    
    # 2) Remove and report any completely embedded segments : sorted alignmentList into filteredList, embeddedList
    
    alignmentList.sort(key=lambda tup: [tup[0],tup[2]])

    filteredList = []
    embeddedList = []
    for key, items in groupby(alignmentList, itemgetter(0)):
        prevEnd = -1
        for eachitem in items:
            if eachitem[3] <= prevEnd:
                embeddedList.append(eachitem)
            else:
                filteredList.append(eachitem)
                prevEnd = eachitem[3]
                
    print " len(filteredList) , len(embeddedList)",  len(filteredList) , len(embeddedList)

    #for eachitem in alignmentList:
    #    print eachitem  
                  
    #print "------------------------------------------"
    #for eachitem in embeddedList:
    #    print eachitem  
                  

    print "filteredList",filteredList
    lenDic = obtainLength(folderName, 'repeat.fasta')
    print "lenDic", lenDic
    lenDic = obtainLength(folderName, 'reference.fasta')
    print "lenDic", lenDic
    
    
    # 3) Compute successor overlaps and report NGaps and overlap histogram 
    overLapList = []
    overlapLenList = []
    for key, items in groupby(filteredList, itemgetter(0)):
        prevEnd = -1
        
        for eachitem in items:
            if prevEnd != -1:
                overLapList.append([key, -eachitem[2]+prevEnd + 1, eachitem[1]])
                overlapLenList.append(-eachitem[2]+prevEnd + 1)
                
            prevEnd = eachitem[3]

        
    #for eachitem in overlapLenList:
    #    print eachitem
    #print overlapLenList
    plt.plot(overlapLenList,color='r')
    plt.plot([0 for i in range(len(overlapLenList)) ], color='b')
    plt.show()
    ### Problem observed  : mis-assembly @ EC, repeat separated out, low percentage overlap
    
    
    
def removeEmbedded():
    print "Hello world"
    folderName = "S_cerivisea/"
    mummerLink = "MUMmer3.23/"
    
    thres = 3
    
    if True:
        command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"self "+ folderName+ "contigs.fasta "+ folderName+ "contigs.fasta"
        os.system(command)
        
        command  = mummerLink +"show-coords -r "+folderName+"self.delta > "+folderName+"selfOut"
        os.system(command)
    
    dataList = extractMumData(folderName, "selfOut")
    dataList= transformCoor(dataList)
    
    lenDic = obtainLength(folderName, 'contigs.fasta')
    
    removeList = []
    for eachitem in dataList:
        match1, match2, name1, name2 = eachitem[4],eachitem[5], eachitem[7], eachitem[8]
        
        if name1 != name2:
            l1, l2 = lenDic[name1], lenDic[name2]
            
            if abs(l1- match1) < thres and abs(l2-match2) > thres:
                removeList.append(name1)
            elif abs(l1- match1) > thres and abs(l2-match2) < thres:
                removeList.append(name2)
            elif abs(l1- match1) < thres and abs(l2- match2) < thres:
                print "Both shortembedd", eachitem
                
    
    
    nameList = []
    for eachitem in lenDic:
        nameList.append(eachitem)

    print len(nameList)
    
    for eachitem in removeList:
        if eachitem in nameList:
            nameList.remove(eachitem)
    print len(nameList)
    
    putListToFileO(folderName, "contigs.fasta", "remainedContigs", nameList)

    
    
def computeConsensus():
    folderName = "Kurt_new/"
    mummerLink = "MUMmer3.23/"
    
    percentMatchList = []
    itemList =[]

    itemList.append( ["Pedobacter/Pedobacter_heparinus_DSM_2366.reference.fa", "Pedobacter/CA_all_all.fa"] )
    itemList.append(["Pedobacter/Pedobacter_heparinus_DSM_2366.reference.fa", "Pedobacter/all_finSCall.fa"] )
    itemList.append( ["Pedobacter/Pedobacter_heparinus_DSM_2366.reference.fa", "Pedobacter/pbj_all_all.fa"] )
        
    itemList.append( ["Pedobacter/Pedobacter_heparinus_DSM_2366.reference.fa", "Pedobacter/CA_50_all.fa"] )
    itemList.append(["Pedobacter/Pedobacter_heparinus_DSM_2366.reference.fa", "Pedobacter/50_finSCall.fa"] )
    itemList.append( ["Pedobacter/Pedobacter_heparinus_DSM_2366.reference.fa", "Pedobacter/pbj_50_all.fa"] )


    itemList.append(["Pedobacter/Pedobacter_heparinus_DSM_2366.reference.fa", "Pedobacter/CA_30_all.fa"] )
    itemList.append( ["Pedobacter/Pedobacter_heparinus_DSM_2366.reference.fa", "Pedobacter/30_finSCall.fa"] )
    itemList.append(["Pedobacter/Pedobacter_heparinus_DSM_2366.reference.fa", "Pedobacter/pbj_30_all.fa"] )
    
    ### hijacking 
    folderName = "Ph1/"
    itemList  = [["reference.fasta", "improved3.fasta"]]
    
    ### End hijacking
    for eachname in itemList:    
        referenceName = eachname[0]
        improvedName = eachname[1]
        if True:
            command = mummerLink +  "/nucmer --prefix=consensus " + folderName+ referenceName+" " + folderName+improvedName
            os.system(command)
        
            command =  mummerLink + "/show-tiling -v 50 -g 50000 -c  consensus.delta > result.txt"
            os.system(command)
        
        f = open("result.txt", 'r')
        tmp = f.readline().rstrip()
        
        G = 0   
        matchLen = 0 
        while len(tmp) > 0:
            if tmp[0] != '>':
                # -426835    814180    687    1241016    99.94    100.00    +    ctg7180000000036quiver
                myInfo  = tmp.split()
                print myInfo
                matchLen = matchLen + min(0, int(myInfo[2])) +int(myInfo[3]) * float(myInfo[4])*float(myInfo[5]) /10000
            else:
                # >ref_NC_001224_ 85779 bases
                myInfo  = tmp.split()
                G = G + int(myInfo[1])
                print G
                
            tmp = f.readline().rstrip() 
        
        print " Percent match : ", matchLen/G
        f.close() 
        percentMatchList.append([referenceName, improvedName, matchLen/G])
        
    for eachresult in percentMatchList:
        print eachresult 
 
 
 
def checkGapLocation():
    folderName = "Ph1/"    
    mummerLink = "MUMmer3.23/"
    gap1 = [1362691, 1363262]
    gap2 = [4282994, 4287906]
    print "Gap 1 length: " , gap1[1] - gap1[0]
    print "Gap 2 length: " , gap2[1] - gap2[0]
    
    if False:
        command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"gapAnalysis "+ folderName+ "reference.fasta "+ folderName+ "phasingSeedName.fasta"
        os.system(command)
        
        command  = mummerLink +"show-coords -r "+folderName+"gapAnalysis.delta > "+folderName+"gapAnalysisOut"
        os.system(command)        

    
    dataList = extractMumData(folderName, "gapAnalysisOut")
    possibleList = []
    
    for eachitem in dataList:
        startPt = eachitem[0]
        endPt = eachitem[1]
        
        '''
        if startPt < gap1[0] and endPt > gap1[1]:
            print "Bridged Gap 1"
            print eachitem
        
        if startPt < gap2[0] and endPt > gap2[1]:
            print "Bridged Gap 2"
            print eachitem
       '''
        
        k1, k2 = 0.333, 0.667
        mid1 = gap2[0]*(1-k1)+ gap2[1]*k1
        mid2 = gap2[0]*(1-k2)+ gap2[1]*k2
        
        assert(mid1 <= mid2)
        if startPt<  gap2[0]:
            possibleList.append(startPt)
            
        if startPt < gap2[0] and endPt >= mid1:
            print "inside" 
        
        if startPt < mid1 and endPt > mid2:
            print "middle" 
        
        if  startPt < mid2 and endPt > gap2[1]:
            print "outside" 
            
    possibleList.sort(reverse = True)
    
    print possibleList[0:10]


def calculateSeg(lengthToCapture, gap1):
    amountToExtend = (lengthToCapture - (gap1[1] - gap1[0]) )/2
    assert(amountToExtend >0)
    seg1 = [gap1[0] - amountToExtend, gap1[1]+amountToExtend]
    return seg1
    

def detectRepeat():
    gap1 = [1362691, 1363262]
    gap2 = [4282994, 4287906]       
    lengthToCapture = 60000
    seg1 = calculateSeg(lengthToCapture, gap1)
    seg2 = calculateSeg(lengthToCapture, gap2)
    
    folderName = "Ph1/"    
    mummerLink = "MUMmer3.23/"
    
    f = open(folderName +"reference.fasta" , 'r')
    G = ""
    
    tmp =f.readline().rstrip()
    print tmp
    tmp =f.readline().rstrip()
    
    while len(tmp) > 0:
        G = G+ tmp
        tmp =f.readline().rstrip()
        
    print len(G)
        
    fout = open(folderName + "repeatCheck.fasta", 'w')
    fout.write(">Seg0\n")
    fout.write(G[seg1[0]:seg1[1]]+"\n")
    fout.write(">Seg1\n")
    fout.write(G[seg2[0]:seg2[1]]+"\n")
    
    
    
    fout.close()
    f.close()
    
    
    
     
#removeEmbedded()    


#checkGapLocation()

def generateASyntheticCase():
    '''
    Output: reference.fasta, contigs.fasta, improved3.fasta
    Algorithm: 
    1) Call the existing subroutines
    '''

#computeConsensus()
#putListToFileO("S_cerivisea/", "improved2.fasta", "repeat", ["Segkk27", "Segkk76"])
#putListToFileO("Ph1/", "improved3.fasta", "repeatAnalysis", ["Segkk0", "Segkk1", "Segkk4"])

#analyseOverlap()

#detectRepeat()




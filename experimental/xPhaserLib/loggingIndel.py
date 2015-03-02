import csv 
import numpy as np 
import common 

def rawDataSave(folderName, motherGen,noisyReads, cleanReads, longRepeatStat, readsMapping ):
    #print "Save"
    #print "motherGen", motherGen
    #print  "noisyReads", noisyReads[0][0]
    #print "cleanReads", cleanReads
    #print  "longRepeatStat", longRepeatStat
    #print  "readsMapping", readsMapping 
    
    # Save  motherGen
    f = open(folderName+ "motherGenome.txt", 'w')
    for eachbase in motherGen:
        f.write(str(eachbase))
    f.close()
    
    # Save cleanReads
    f = open(folderName+ "cleanShortReads.txt", 'w')
    for eachread, index in  zip(cleanReads[0], range(len(cleanReads[0]))):
        for eachbase in eachread:
            f.write(str(eachbase))
        if index != len(cleanReads[0])-1:
            f.write("\n")
    f.close()
    
    f = open(folderName+ "cleanLongReads.txt", 'w')
    for eachread, index in zip(cleanReads[1], range(len(cleanReads[1]))):
        for eachbase in eachread:
            f.write(str(eachbase))
        if index != len(cleanReads[1])-1:
            f.write("\n")
    f.close()
    
    # Save noisyReads
    f = open(folderName+ "noisyShortReads.txt", 'w')
    for eachread, index in  zip(noisyReads[0], range(len(noisyReads[0]))):
        for eachbase in eachread:
            f.write(str(eachbase))
        if index != len(noisyReads[0])-1:
            f.write("\n")
    f.close()
    
    f = open(folderName+ "noisyLongReads.txt", 'w')
    for eachread, index in zip(noisyReads[1], range(len(noisyReads[1]))):
        for eachbase in eachread:
            f.write(str(eachbase))
        if index != len(noisyReads[1])-1:
            f.write("\n")
    f.close()
    
    # Save longRepeatStat
    f = open(folderName + "longRepeatStat.csv","wb")
    mywriter = csv.writer(f)
    
    mywriter.writerow(["startLoc1", "startLoc2", "length"])
    for eachitem in longRepeatStat:
        mywriter.writerow(eachitem)
    
    f.close()
    
    # Save readsMapping
    f = open(folderName + "readsMapping.csv","wb")
    mywriter = csv.writer(f)
    
    mywriter.writerow(["readNum", "startLoc", "forward/backward", "longshort"])
    for eachitem in readsMapping[0]:
        mywriter.writerow(eachitem + [0])
        
    for eachitem in readsMapping[1]:
        mywriter.writerow(eachitem + [1])
    
    f.close()    
def transformFromFixedListToNumpy(listOfReads):
    L = len(listOfReads[0])
    N = len(listOfReads)
    
    
    returnReads = np.zeros(N*L, dtype = np.int8).reshape(N,L)
    
    for i in range(N):
        for j in range(L):
            returnReads[i][j] = listOfReads[i][j]
    
    return returnReads
    
def transformFromVaryListToNumpy(listOfReads):
    returnList = []
    for i in range(len(listOfReads)):
        L = len(listOfReads[i])
        tempread = np.zeros(L, dtype = np.int8)
        for j in range(L):
            tempread[j] = listOfReads[i][j]
        returnList.append(tempread)
        
    return returnList
            
        
def rawDataLoad(folderName ):
    #print "Save"
    #print "motherGen", motherGen
    #print  "noisyReads", noisyReads[0][0]
    #print "cleanReads", cleanReads
    #print  "longRepeatStat", longRepeatStat
    #print  "readsMapping", readsMapping 
    motherGen,noisyReads, cleanReads, longRepeatStat, readsMapping = [] , [] , [] , [] , [] 
    
    # Load  motherGen
    f = open(folderName+ "motherGenome.txt", 'r')
    temp = f.read()
    motherGen = np.zeros(len(temp), dtype = np.int8) 
    
    for index in range(len(temp)):
        motherGen[index] = int(temp[index])
    f.close()
    
    # Load cleanReads
    f = open(folderName+ "cleanShortReads.txt", 'r')
    shortReadsList = [] 
    temp = f.readline()
    while ( len(temp) > 0 ):
        if temp[-1] == '\n':
            shortReadsList.append(temp[0:-1])
        else: 
            shortReadsList.append(temp)
        temp = f.readline()
        
    cleanShortReads = transformFromFixedListToNumpy(shortReadsList)
    f.close()
    
    f = open(folderName+ "cleanLongReads.txt", 'r')
    longReadsList = [] 
    temp = f.readline()
    while ( len(temp) > 0 ):
        if temp[-1] == '\n':
            longReadsList.append(temp[0:-1])
        else: 
            longReadsList.append(temp)
        temp = f.readline()
        
    cleanLongReads = transformFromFixedListToNumpy(longReadsList)
    f.close()
    
    cleanReads = [cleanShortReads, cleanLongReads]
    
    # Load noisyReads
    f = open(folderName+ "noisyShortReads.txt", 'r')
    shortReadsList = [] 
    temp = f.readline()
    while ( len(temp) > 0 ):
        if temp[-1] == '\n':
            shortReadsList.append(temp[0:-1])
        else: 
            shortReadsList.append(temp)
        temp = f.readline()
        
    noisyShortReads = transformFromVaryListToNumpy(shortReadsList)
    f.close()
    
    f = open(folderName+ "noisyLongReads.txt", 'r')
    longReadsList = [] 
    temp = f.readline()
    while ( len(temp) > 0 ):
        if temp[-1] == '\n':
            longReadsList.append(temp[0:-1])
        else: 
            longReadsList.append(temp)
        temp = f.readline()
        
    noisyLongReads = transformFromVaryListToNumpy(longReadsList)
    f.close()
    
    noisyReads = [noisyShortReads, noisyLongReads]
    
    # Save longRepeatStat
    f = open(folderName + "longRepeatStat.csv","r")
    myreader = csv.reader(f)
    
    counter = 0
    for eachrow in myreader :
        if counter == 0:
            counter =1
        else:
            longRepeatStat.append([int(eachrow[0]), int(eachrow[1]) , int(eachrow[2])])

    f.close()
    
    # Save readsMapping
    f = open(folderName + "readsMapping.csv","r")
    myreader = csv.reader(f)

    shortReadMapping = []
    longReadMapping = [] 
    
    counter = 0 
    for eachrow in myreader : 
        if counter  == 0 :
            counter = 1 
        else: 
            if eachrow[3] == '0' :
                shortReadMapping.append([int(eachrow[0]), int(eachrow[1]), int(eachrow[2])])
            elif eachrow[3] == '1':
                longReadMapping.append([int(eachrow[0]), int(eachrow[1]), int(eachrow[2])])
            else:
                print "Errrrrror "
    
    readsMapping = [shortReadMapping,longReadMapping ]       
    f.close()
    
    return motherGen,noisyReads, cleanReads, longRepeatStat, readsMapping






def logVoteTable(folderName,voteTableList):
    f=open(folderName + "voteTableLog.txt", "w")
    for eachvoteTable in voteTableList:
        # Write the long reads 
        f.write(str(eachvoteTable.index)+ ":")
        for eachbase in eachvoteTable.longread:
            f.write(str(eachbase))
        f.write(":")
        
        
        # Write the delCount
        for eachcount, i in zip(eachvoteTable.delCount, range(len(eachvoteTable.delCount))):
            if i != 0:
                f.write(",")
            f.write(str(eachcount))
            
        f.write(":")
        # Write the confirmCount
        for eachcount, i in zip(eachvoteTable.confirmCount, range(len(eachvoteTable.confirmCount))):
            if i != 0:
                f.write(",")
            f.write(str(eachcount))
        
        f.write(":")
        
        # Write subCount 
        for eachcount, i in zip(eachvoteTable.subCount, range(len(eachvoteTable.subCount))):
            if i != 0:
                f.write(",")
            f.write(str(eachcount))
        
        f.write(":")
        
        # Write the insCount
        for eachcount, i in zip(eachvoteTable.insCount, range(len(eachvoteTable.insCount))):
            if i != 0:
                f.write(",")
            f.write(str(eachcount))
        
        
        f.write(":")
        
        for eachcount, i in zip(eachvoteTable.confirmNoIns, range(len(eachvoteTable.confirmNoIns))):
            if i != 0 :
                f.write(",")
            f.write(str(eachcount))
        
            
        f.write("\n")
        

    
    f.close()
    
def loadVoteTable(folderName):
    f = open(folderName+"voteTableLog.txt", 'r')
    
    temp = f.readline()
    voteTableList = []
    while ( len(temp) > 0 ):
        dataList = temp.split(":")
        
        index = int(dataList[0])
        longread = np.zeros(len(dataList[1]), np.int8)
        for j in range(len(dataList[1])):
            longread[j] = dataList[1][j]
        
        subdataList = dataList[2].split(",")
        delCount = []
        for j in range(len(subdataList)):
            delCount.append(int(subdataList[j]))
        
        subdataList = dataList[3].split(",")
        confirmCount = []
        for j in range(len(subdataList)):
            confirmCount.append(int(subdataList[j]))
            
            
        subdataList = dataList[4].split("],")
        #print subdataList
        subCount = []
        for j in range(len(subdataList)):
            subsubdataList = subdataList[j].split(",")
            tempsub = []
            for k in range(len(subsubdataList)):
                testitem = subsubdataList[k]
                if testitem[0] == "[":
                    testitem = testitem[1:len(testitem)]
                if testitem[-1] == "]":
                    testitem = testitem[0: -1]
                tempsub.append(int(testitem))
                
            subCount.append(tempsub)
        
        #print dataList[5]
        subdataList = dataList[5].split("],[")
        insCount = []
        
        #print subdataList
        
        #print "... insCount " 

        for j in range(len(subdataList)):
            
            subsubList = subdataList[j].split("], [")
            tempentryList = []
             
            # Transform 
            
            for k in range(len(subsubList)):
                if len(subsubList[k]) > 0 and  subsubList[k][0] == '[':
                    subsubList[k] = subsubList[k][1:len(subsubList[k])]
                if len(subsubList[k]) > 0 and subsubList[k][-1] == ']':
                    subsubList[k] = subsubList[k][0:-1]
            
            # Fill in 
            for k in range(len(subsubList)):
                if len(subsubList[k]) > 0 :
                    tempitem = subsubList[k].split(", ")
                    if tempitem[0][0] == '[':
                        tempitem[0] = tempitem[0][1:len(tempitem)]
                    if tempitem[2][-1] == ']':
                        tempitem[2] = tempitem[2][0:-1]

                    tempentryList.append([int(tempitem[0]),int(tempitem[1]),int(tempitem[2])])
                  
            
            insCount.append(tempentryList)
                        
        confirmNoIns = []
        subdataList = dataList[6][0:-1].split(",") 
        
        for j in range(len(subdataList)):
            confirmNoIns.append( int(subdataList[j]))
                        
             
        #print "index", index 
        #print "longread", len(longread)
        #print "delCount", len(delCount)
        #print "confirmCount", len(confirmCount)
        #print "subCount", len(subCount)
        #print "insCount", len(insCount)
        #print "confirmNoIns", len(confirmNoIns)
        
        tempVoteTable = common.voteTable(index, longread) 
        tempVoteTable.initVoteTable()
        tempVoteTable.delCount = delCount
        tempVoteTable.confirmCount = confirmCount
        tempVoteTable.subCount = subCount
        tempVoteTable.insCount = insCount 
        tempVoteTable.confirmNoIns = confirmNoIns
        
        voteTableList.append(tempVoteTable)
        temp = f.readline()

    f.close()
    
    return voteTableList         


def logrMapping(rMapping, folderName):
    f= open(folderName+"rMapping.csv", "wb")
    
    mywriter = csv.writer(f)
    
    mywriter.writerow(["Long reads", "Short reads"])
    
    for eachitem in rMapping:
        mywriter.writerow([eachitem[0], eachitem[1]])
        
    f.close()
    
    
    
def logshortToLongMap(folderName,shortToLongMap):
    print folderName
    f = open(folderName+"shortToLong.csv","wb")
    mywriter = csv.writer(f)
    mywriter.writerow( [ "indexlong", "indexshort",  "jstart", "jend", "istart", "iend"] )
    for eachitem in shortToLongMap:
        mywriter.writerow(eachitem)
         #[ indexlong, indexshort,  jstart, jend, istart, iend] 
    f.close()

    
                    
def loadShortToLongMap(folderName):
    shortToLongMap = []
    f = open(folderName+"shortToLong.csv","r")
    myreader = csv.reader(f)
    
    count = 0
    for eachrow in myreader:
        if count == 0:
            count = 1
        else:
            tempitem = []
            for eachitem in eachrow:
                tempitem.append(int(eachitem))
                
            shortToLongMap.append(tempitem)

    f.close()
    return shortToLongMap


def loadSegList(folderName):
    f = open(folderName + "segmentFile.txt", 'r')
    temp = f.readline()
    listOfSegments = []
    while ( len(temp) > 0 ):
        listOfSegments.append(temp)
        temp = f.readline()
    
    
    leftSegList = []
    rightSegList = []
    
    if len(listOfSegments) == 4:
        for segment in listOfSegments[0:2]:
            tmpSeg = []
            for eachbase in segment:
                if eachbase != '\n':
                    tmpSeg.append(int(eachbase))
            leftSegList.append(tmpSeg)

        for segment in listOfSegments[2:4]:
            tmpSeg = []
            for eachbase in segment:
                if eachbase != '\n':
                    tmpSeg.append(int(eachbase))
            rightSegList.append(tmpSeg)
            
    f.close()
    
    return leftSegList, rightSegList
    
    
    

def logSegList(folderName , leftSegList, rightSegList):
    f = open(folderName + "segmentFile.txt", 'w')
    for eachitem in leftSegList:
        for eachbase in eachitem:
            f.write(str(eachbase))
    
        f.write("\n")
    for eachitem in rightSegList:
        for eachbase in eachitem : 
            f.write(str(eachbase))
        f.write("\n")
        
    f.close()
    
    
    
def logFinishedContigs(contigList, folderName):
    f = open(folderName+ "finishedContigs.txt", 'w')    
    for eachcontig in contigList:
        for eachbase in eachcontig:
            f.write(str(eachbase))
        f.write("\n")
    f.close()
            
    

def logHTVote(HTList, folderName):
    f = open(folderName+ "HTList.txt", 'w')
    for eachvote in HTList:
        f.write(str(eachvote))
    f.write("\n")
    f.close
    
    

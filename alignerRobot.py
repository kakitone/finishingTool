import os
import sys
import houseKeeper
from multiprocessing import Pool
import time

def extractMumData(folderName, fileName):
    # "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    f = open(folderName + fileName, 'r')
    dataList = []
    
    for i in range(6):
        tmp = f.readline()

    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr = info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        
        matchLenArr = info[2].split()
    
        matchLen1 = int(matchLenArr[0])
        matchLen2 = int(matchLenArr[1])    
        percentMatch = float(info[3])
        
        
        helperStart, helperEnd = int(firstArr[0]), int(firstArr[1])
        readStart, readEnd = int(filterArr[0]) , int(filterArr[1])
        
        helperName = rdGpArr[0].rstrip().lstrip()
        readName = rdGpArr[1].rstrip().lstrip()
        
        dataList.append([helperStart, helperEnd , readStart, readEnd, matchLen1, matchLen2, percentMatch, helperName, readName ])
    
        tmp = f.readline().rstrip()
                
    f.close()
    
    return dataList



def useMummerAlign(mummerLink, folderName, outputName, referenceName, queryName, specialForRaw = False, specialName = "", refinedVersion= False):
    nucmerMummer(specialForRaw, mummerLink, folderName, outputName, referenceName, queryName, refinedVersion)
    showCoorMummer(specialForRaw, mummerLink, folderName, outputName, specialName)
    
    

def nucmerMummer(specialForRaw, mummerLink, folderName, outputName, referenceName, queryName,refinedVersion):
    if not refinedVersion:
        if not specialForRaw:
            if houseKeeper.globalFast:
                command = mummerLink + "nucmer -b 50  --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " + folderName + queryName
            else:
                command = mummerLink + "nucmer --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " + folderName + queryName
        else:
            if houseKeeper.globalFast:
                command = mummerLink + "nucmer -b 50  --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " +  queryName
            else:
                command = mummerLink + "nucmer --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " +  queryName
    else:
        command = mummerLink + "nucmer -l 10 --maxmatch -p " + folderName + outputName +" " + folderName +referenceName +" " + queryName
        
    os.system(command)
    

def showCoorMummer(specialForRaw, mummerLink, folderName, outputName, specialName):
    if not specialForRaw:
        command = mummerLink + "show-coords -r " + folderName + outputName + ".delta > " + folderName + outputName + "Out"
    else:
        command = mummerLink + "show-coords -r " + folderName + outputName + ".delta > " + folderName + specialName 
        
    os.system(command)


def combineMultipleCoorMum(specialForRaw, mummerLink, folderName, outputName, specialName, numberOfFiles):
    print ""
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
        showCoorMummer(specialForRaw, mummerLink, folderName, outputName+indexOfMum, specialName+indexOfMum)
    
    
    command =  "head -5 "+ folderName + specialName +"01" + "> " + folderName +specialName
    os.system(command)

    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
    
        command = " tail -n+6 "+ folderName + specialName +indexOfMum +">> " + folderName +specialName
        os.system(command)

def zeropadding(i):
    tmpi = ""

    if i < 10:
        tmpi = "0" + str(i)
    else:
        tmpi = str(i)
    return tmpi
   
   
def useMummerAlignBatch(mummerLink, folderName, workerList, nProc ,specialForRaw = False, refinedVersion = False):
    # Format for workerList : [[outputName, referenceName, queryName, specialName]... ]
    # nProc : a parameter on how many threads should be created each time
    # Goal : parallelize this part  
    if not houseKeeper.globalLarge:
        p = Pool(processes=nProc)
        results = []
        
        for eachitem in workerList:
            outputName, referenceName, queryName, specialName = eachitem
            results.append(p.apply_async(useMummerAlign, args=(mummerLink, folderName, outputName, referenceName, queryName, specialForRaw , specialName, refinedVersion)))
        
        outputlist = [itemkk.get() for itemkk in results]
        print  len(outputlist)
        p.close()
    else:
        '''
        a) Split
        b) align several times
        c) join the query
        raw_reads.part-01
        '''
        p = Pool(nProc)
        results = []
        numberOfFiles = 10
        
        for eachitem in workerList:
            print eachitem
            outputName, referenceName, queryName, specialName = eachitem[0], eachitem[1], eachitem[2] , eachitem[3]
        
            
            bindir =  os.path.abspath(os.path.dirname(sys.argv[0]))   
            command = bindir + "/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + folderName + referenceName
            os.system(command)
            

            if specialForRaw : 
                queryNameMod = queryName
            else:
                queryNameMod = folderName + queryName

            command = bindir + "/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + queryNameMod
            os.system(command)
        
            
        for eachitem in workerList:   
            outputName, referenceName, queryName, specialName = eachitem[0], eachitem[1], eachitem[2] , eachitem[3]
            for i in range(1, numberOfFiles+1):
                for j in range(1, numberOfFiles+1):
                    if specialForRaw : 
                        tmpRefName , tmpQryName = referenceName[0:-6] + ".part-" + zeropadding(i) +".fasta",  queryName[0:-6] + "-" + zeropadding(j) + ".fasta"
                    else:
                        tmpRefName , tmpQryName = referenceName[0:-6] + ".part-" + zeropadding(i) +".fasta",  queryName[0:-6] + ".part-" + zeropadding(j) + ".fasta"
                    
                    
                    results.append(p.apply_async(nucmerMummer, args =(specialForRaw, mummerLink, "", folderName + outputName +zeropadding(i)+zeropadding(j), tmpRefName, tmpQryName)))
      
        
        outputlist = [itemkk.get() for itemkk in results]
        print len(outputlist)
        p.close()
        
        for eachitem in workerList:
            outputName, referenceName, queryName, specialName = eachitem           
            if not specialForRaw:
                outNameMod =  folderName + outputName + "Out"
            else:
                outNameMod = folderName + specialName 
        
        
            tmpName = folderName + outputName +zeropadding(1)+zeropadding(1) + ".delta"
        
            command = mummerLink + "show-coords -r " + tmpName + "| head -5 > " + outNameMod
            os.system(command)
            
            for i in range(1, numberOfFiles+1):
                for j in range(1, numberOfFiles+1):
                    
                    tmpName = folderName + outputName +zeropadding(i)+zeropadding(j) + ".delta"
                    command = mummerLink + "show-coords -r " + tmpName + "| tail -n+6 >> " + outNameMod
                    os.system(command)
        

def transformCoor(dataList):
    # "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
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

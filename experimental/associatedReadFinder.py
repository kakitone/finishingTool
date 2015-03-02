

import os

from operator import itemgetter
from itertools import groupby

from finisherSCCoreLib import IORobot
from finisherSCCoreLib import alignerRobot
from finisherSCCoreLib import houseKeeper

import abunHouseKeeper


def checkSatisfy(eachitem, lenDic):
    # "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    threshold = 40
     
    seedName = eachitem[-1]
    detectedName = eachitem[-2]
    
    # seedIndices = seedName.split("/")[-1].split("_")
    # detectIndices = detectedName.split("/")[-1].split("_")
    # lenSeed = int(seedIndices[1]) - int(seedIndices[0])
    # lenDetect = int(detectIndices[1]) - int(detectIndices[0])
    
    lenSeed = lenDic[seedName]
    lenDetect = lenDic[detectedName]
    
    checkSeed = False
    checkDetect = False 
    
    if eachitem[0] < threshold or eachitem[1] > lenSeed - threshold: 
        checkSeed = True
    
    if min(eachitem[2], eachitem[3]) < threshold or max(eachitem[2], eachitem[3]) > lenDetect - threshold:
        checkDetect = True
    
    if checkSeed and checkDetect: 
        return True
    else:
        return False


def getAllAssociatedReads(folderName, mummerLink,forFastaName):
    '''
    Input : relatedReads.fasta, raw_reads.fasta 
    Output : all_associated_reads.fasta
    
     Algorithm : 
        a) Get all the associated reads
        b) Loop for N=1 times : ==> this correspond 4 reads to link between the bridge in total
            i) Align the raws and tmp_seedReads
            ii) Put the new reads into the SeedReads
    '''
    
    header, referenceFile, queryFile = "seedReads", forFastaName + ".fasta" , "raw_reads.fasta"
    command = "cp " + folderName + "relatedReads.fasta " + folderName + referenceFile
    os.system(command)
    
    N = abunHouseKeeper.abunGlobalReadSearchDepth
    
    print "N: ", N
    if N >0 :
        for trial in range(N):
            print "trial", trial
            numberOfFiles = 20
            
            if True: 
                workerList = []
                
                for dummyI in range(1, numberOfFiles + 1):
                    indexOfMum = ""
                    if dummyI < 10:
                        indexOfMum = "0" + str(dummyI)
                    else:
                        indexOfMum = str(dummyI)
                    
                    outputName, referenceName, queryName, specialName= header+indexOfMum, referenceFile, "raw_reads.part-"+ indexOfMum + ".fasta",  header + indexOfMum
                    workerList.append([outputName, referenceName, queryName, specialName])
    
                alignerRobot.useMummerAlignBatch(mummerLink, folderName, workerList, houseKeeper.globalParallel ,False)
            
            dataList = []
            
            for i in range(1, 1+numberOfFiles): 
                if i < 10:
                    indexOfMum = "0" + str(i)
                else:
                    indexOfMum = str(i)
                dataList = dataList+ alignerRobot.extractMumData(folderName, header+ str(indexOfMum)+"Out")
            
            
            filterList = []
            
            lenDicRR = IORobot.obtainLength(folderName, queryFile)
            
            print "len(dataList)", len(dataList)
            for eachitem in dataList:
                if checkSatisfy(eachitem, lenDicRR):
                    filterList.append(eachitem)
                
            filterList.sort(key=itemgetter(-1))
            newReads = []
            
            for key, items in groupby(filterList, itemgetter(-1)):
                newReads.append(key)
                                        
            
            f = open(folderName + forFastaName + ".txt", 'w')
            
            for eachitem in newReads:
                f.write(eachitem + "\n")
            f.close()
                
            command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' " + folderName + forFastaName + ".txt " + folderName + "raw_reads.fasta > " + folderName + forFastaName + ".fasta"
            os.system(command)
    else:
        os.system("cp " + folderName + "relatedReads.fasta " + folderName + forFastaName + ".fasta")




    
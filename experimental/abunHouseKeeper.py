from itertools import groupby
import os 

abunGlobalAvoidrefine = False
abunGlobalReadSearchDepth = 1
abunGlobalRRDisable = False

def replaceFiles( folderName, replacedName) :
    commandList = []
    commandList.append("cp " + folderName + "improved3.fasta " + folderName + "improved3_backup.fasta")
    commandList.append("cp " + folderName + "improved3_Double.fasta " + folderName + "improved3_backup.fasta")
    
    commandList.append("cp " + folderName + replacedName + " "+folderName + "improved3.fasta")
    commandList.append("cp " + folderName + replacedName[0:-6]+"_Double.fasta " + folderName + "improved3_Double.fasta")
    
    for eachcommand in commandList:
        print eachcommand
        os.system(eachcommand)

def parseEdgeNameToID(name, mytype):

    if mytype == 'C':
        dataInfo = name[6:].split('_')
    elif mytype == 'R':
        dataInfo = name[4:].split('_')
        
    contigNum = int(dataInfo[0])
    
    if dataInfo[1] == 'p':
        id = contigNum * 2
    elif dataInfo[1] == 'd':
        id = contigNum * 2 + 1
        
    return id 


def getDistinct(myList):
    newList = [] 
    myList.sort()
    
    for key, items in groupby(myList):
        newList.append(key)
    
    return newList


def filterData(dataList, lenDic):
    newDataList = []
    for eachitem in dataList:
        if headTailMatch(eachitem, lenDic):
            newDataList.append(eachitem)

    return newDataList

def filterDataIdentical(dataList, lenDic):
    newDataList = []
    for eachitem in dataList:
        if not identicalItem(eachitem, lenDic):
            newDataList.append(eachitem)

    return newDataList

def identicalItem(eachitem, lenDic):
    start1, end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]
    thres = 20 
    
    endPt1 = lenDic[eachitem[-2]] - end1
    endPt2 = lenDic[eachitem[-1]] - end2
    
    if abs(start1-start2) < thres and abs(endPt1 - endPt2) < thres:
        return True
    else:
        return False

def headTailMatch(eachitem, lenDic):
    start1, end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]
    l1, l2 = lenDic[eachitem[-2]], lenDic[eachitem[-1]]
    thres = 10 
    
    
    diffContig , forwardStrand, headtailoverlap = False, False , False
    
    if eachitem[-2] == eachitem[-1]:
        diffContig = False
    else:
        diffContig = True

    if start2 > end2 : 
        forwardStrand = False
    else:
        forwardStrand = True
    
    
    if (start1 <= thres and end2 >= l2 - thres) or (end1 >= l1 - thres and start2 <= thres): 
        headtailoverlap = True
    else:
        headtailoverlap = False
        
    
    if diffContig and forwardStrand and headtailoverlap:
        return True
    else:
        return False
    
    
    return True


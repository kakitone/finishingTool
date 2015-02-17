import os 
import sys

globalFast = False
globalParallel = 1
globalLarge = False


def quastEvaluate(folderName, quastLink, originalName, improvedNameList, referenceName):
    
    # ./quast.py ~/git/myFinisher/finishingTool/S_cerivisea/contigs.fasta  ~/git/myFinisher/finishingTool/S_cerivisea/improved.fasta -R ~/git/myFinisher/finishingTool/S_cerivisea/reference.fasta
    bindir = os.path.abspath(os.path.dirname(sys.argv[0]))   
    header = bindir + "/" + quastLink + "quast.py" + " "
    
    originalContigPath = folderName + originalName + " "
    improvedContigPath = "" 
    for eachname in improvedNameList:
        improvedContigPath = improvedContigPath + folderName + eachname + " "
        
    
    referencePath = "-R " + folderName + referenceName + " "
    
    command = header + originalContigPath + improvedContigPath + referencePath 
    
    os.system(command)
    
    bindir = os.path.abspath(os.path.dirname(sys.argv[0]))    
    command = "cp " + bindir + "/quast_results/latest/report.txt " + folderName + "assemblyAssessment.txt"
    os.system(command)

def compareGraphUnitTest(G, G2):
    assert(len(G.graphNodesList) == len(G2.graphNodesList))
    for index in range(len(G.graphNodesList)):
        assert(G.graphNodesList[index].nodeIndex == G2.graphNodesList[index].nodeIndex)
        assert(G.graphNodesList[index].nodeIndexList == G2.graphNodesList[index].nodeIndexList)
        assert(G.graphNodesList[index].overlapList == G2.graphNodesList[index].overlapList)
        assert(G.graphNodesList[index].listOfPrevNodes == G2.graphNodesList[index].listOfPrevNodes)
        assert(G.graphNodesList[index].listOfNextNodes == G2.graphNodesList[index].listOfNextNodes)
        assert(G.graphNodesList[index].visited == G2.graphNodesList[index].visited)


def performMapping(folderName, mummerLink, mapcontigsname):
    print "performMapping"
    info = mapcontigsname.split("_")
    
    oldRef, newRef = info[0], info[1]
    print oldRef, newRef
    command = mummerLink + "nucmer --maxmatch --nosimplify  -p " + folderName + "outMapper " + folderName + newRef + " " + folderName + oldRef
    os.system(command)
    
    command = mummerLink + "show-tiling -v 50 -g 50000 -c " + folderName + "outMapper.delta > " + folderName + "mappingResults.txt"
    os.system(command)
    
    command = "more " + folderName + "mappingResults.txt"
    os.system(command)

    

# ## 7) Compare with reference (I: improved.fasta, improved2.fasta, reference.fasta ; O : assembly assessment report )

def compareWithReference(folderName , mummerLink):
    print "compareWithReference"
    quastEvaluate(folderName, "quast-2.3/", originalName="contigs.fasta", improvedNameList=["noEmbed.fasta", "improved.fasta", "improved2.fasta", "improved3.fasta"] , referenceName="reference.fasta")
    # quastEvaluate(folderName, "quast-2.3/", originalName = "contigs.fasta", improvedNameList= ["noEmbed.fasta", "improved.fasta"] , referenceName= "reference.fasta" )
    

def checkingPath(folderName, mummerLink):
    
    pathExists = True
    newFolderName, newMummerLink = "", ""
    
    if folderName[-1] == "/":
        newFolderName = folderName
    else:
        newFolderName = folderName + "/" 
    
    if mummerLink[-1] == "/":
        newMummerLink = mummerLink
    else:
        newMummerLink = mummerLink + "/" 
    
    if not os.path.isdir(newFolderName):
        pathExists = False
        print "Not exists : " + newFolderName
    
    if not os.path.isdir(newMummerLink):
        pathExists = False
        print "Not exists : " + newMummerLink
    
    if not os.path.exists(newFolderName + "contigs.fasta"):
        pathExists = False
        print "Not exists : " + newFolderName + "contigs.fasta"
    
    if not os.path.exists(newFolderName + "raw_reads.fasta"):
        pathExists = False
        print "Not exists : " + newFolderName + "raw_reads.fasta"
    
    
    return pathExists, newFolderName, newMummerLink


    
def reverseComplement(myStr):
    myNewStr = myStr[::-1]
    myNewStr2 = ""
    for i in range(len(myNewStr)):
        if myNewStr[i] == 'A' or myNewStr[i] == 'a':
            myNewStr2 += 'T'
            
        elif  myNewStr[i] == 'T' or myNewStr[i] == 't':
            myNewStr2 += 'A'
            
        elif  myNewStr[i] == 'C' or myNewStr[i] == 'c':
            myNewStr2 += 'G'
            
        elif  myNewStr[i] == 'G' or myNewStr[i] == 'g':
            myNewStr2 += 'C'
        elif myNewStr[i] == 'N' or myNewStr[i] == 'n':
            myNewStr2 += 'N'
        else:
            print myNewStr[i]
            
    return myNewStr2



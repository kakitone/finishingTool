import os


def mainFlow():
    folderName = "Ph1/"
    genomeName = "reference.fasta"
    
    os.chdir("gepard-1.30")
    os.system("./gepardcmd.sh -seq1 ../"+folderName+genomeName+" -seq2 ../"+folderName+"contigs.fasta -matrix matrices/edna.mat -outfile "+"../"+folderName+"original.png")
    os.chdir("../")  

    os.chdir("gepard-1.30")
    os.system("./gepardcmd.sh -seq1 ../"+folderName+genomeName+" -seq2 ../"+folderName+"allInOne2.fasta -matrix matrices/edna.mat -outfile "+"../"+folderName+"finalAllInOne.png")
    os.chdir("../")    

    for i in range(2):
        os.chdir("gepard-1.30")
        os.system("./gepardcmd.sh -seq1 ../"+folderName+genomeName+" -seq2 ../"+folderName+"improvedContig2_"+str(i)+".fasta -matrix matrices/edna.mat -outfile "+"../"+folderName+"Contig"+ str(i)+ ".png")
        os.chdir("../")  

def mainFlowDor():
    print "Hello"

    folderName = "S_cerivisea/"
    mummerLink  = "MUMmer3.23/"
    
    command =mummerLink + "nucmer -mumreference -c 1000 -l 100 -banded -d 10 -p "+folderName+"out "+folderName+"improved3.fasta "+folderName+"reference.fasta"
    os.system(command)
    
    command  = mummerLink +"delta-filter -1 "+folderName+"out.delta > "+folderName+"out.1delta"
    os.system(command)

    command  = mummerLink +"mummerplot -large -fat -t png "+folderName+"out.1delta"
    os.system(command)    
    
mainFlowDor()
#mainFlow()
import os


def mainFlow():
    folderName = "EcoliTestRun/"
    genomeName = "reference.fasta"
    
    os.chdir("gepard-1.30")
    os.system("./gepardcmd.sh -seq1 ../"+folderName+genomeName+" -seq2 ../"+folderName+"contigs.fasta -matrix matrices/edna.mat -outfile "+"../"+folderName+"original.png")
    os.chdir("../")  

    os.chdir("gepard-1.30")
    os.system("./gepardcmd.sh -seq1 ../"+folderName+genomeName+" -seq2 ../"+folderName+"allInOne2.fasta -matrix matrices/edna.mat -outfile "+"../"+folderName+"finalAllInOne.png")
    os.chdir("../")    

    for i in range(5):
        os.chdir("gepard-1.30")
        os.system("./gepardcmd.sh -seq1 ../"+folderName+genomeName+" -seq2 ../"+folderName+"improvedContig2_"+str(i)+".fasta -matrix matrices/edna.mat -outfile "+"../"+folderName+"Contig"+ str(i)+ ".png")
        os.chdir("../")  


mainFlow()
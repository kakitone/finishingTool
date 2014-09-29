### This page describe the interface and pipeline for the xpmbfinisher
### Input : contigs.fasta, raw_reads.fasta
### Output : improved2.fasta 

import breakerLib
import embedLib
import dbLib
import gapLib
import branchLib
import MBALib
import xPhaseLib
import eulerLib
import cleanLib



def xpmbPipeline():
    folderName  = "" 
    mummerLink = ""
    
    ### break mis-assembled junctions (I: contigs.fasta , O: contigsBroken.fasta )
    breakerLib.breakWrongMerge(folderName,mummerLink)
    
    ### remove completely embedded contigs (I: contigsBroken.fasta, O: contigsNoEmbed.fasta )
    embedLib.removeEmbedded(folderName,mummerLink)
    
    ### mummer Form memory efficient graph (I:contigsNoEmbed.fasta, O:initialGraph.g , endContigs.fasta, fmapping.csv )
    dbLib.formDBGraph(folderName,mummerLink)

    ### streaming to fill gaps (I:raw_reads.fasta, initialGraph.g, fmapping.csv,  O:gapGraph.g, fmapping2.csv)
    gapLib.gapFill(folderName,mummerLink)
    
    ### branchClearing (I:gapGraph.g, O:noBranchGraph.g )
    branchLib.branchClearing(folderName, mummerLink)
    
    ### repeat resolution by MBA (I:noBranchGraph.g, fmapping2.csv , O:bridgedGraph.g)
    MBALib.MBA(folderName,mummerLink)
    
    ### repeat resolution by xPhasing (I:bridgedGraph.g, fmapping2.csv, raw_reads.fasta , O: phasedGraph.g)
    xPhaseLib.xphasing(folderName,mummerLink)
    
    ### euler cycle fining (I: phasedGraph.g , O: sequence.txt )
    eulerLib.ECFinding(folderName,mummerLink)
    
    ### fine-tuning and polishing (I: sequence.txt, contigsNoEmbed.fasta, fmapping2.csv, O: improved2.fasta)
    cleanLib.polishing(folderName,mummerLink)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
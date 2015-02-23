FinisherSC
=============

## Introduction and contact ##
If you encounter any issue, please feel free to directly contact me at kklam@eecs.berkeley.edu
My name is  Ka-Kit Lam. 

Our paper is located at http://biorxiv.org/content/early/2014/11/03/010215

## Command to run the tool ##

Here is the command to run the tool:

	python finisherSC.py destinedFolder mummerPath
	
If you are running on server computer and would like to use multiple threads, then the following commands can generate 20 threads to run FinisherSC. 

	python finisherSC.py -par 20 destinedFolder mummerPath

## Input and output files ##

Note that it is assumed that you have "contigs.fasta" , "raw_reads.fasta" in the destinedFolder and mummerPath is your path to Mummer. Therefore, please rename the files in the destinedFolder if needed. For "finisherSC.py", the key output file is **improved3.fasta** in destinedFolder. 


## A sample test run ##
Here is an example run with the pre-installed data files and softwares(e.g. Mummer, Gepard) [You may go to http://www.eecs.berkeley.edu/~kakitone/package.zip  to download the the whole package to have a quick test. But we note here that all data and softwares besides FinisherSC have their copyright belonging to their original authors. We provide the sandbox here for users' convenience only. Moreover, all the updates and versioning will be at GitHub page. Users are encouraged to git clone our project for later usage after they successfully go through the example run].

1. You may run the test cases as python finisherSC.py EcoliTestRun/ MUMmer3.23/
2. After the processing, take the file whose name is **improved3.fasta** in destinedFolder and this is your output file.
3. If you want to visualize the output against the reference, you can use software like Gepard.

## Options ##
a)Fast mode \[ -f True \] (e.g. for large and complex genomes , you may want to use this mode to get results faster; However, there will be a bit of computation vs quality trade-off if you want to use this mode; We suggest that you use this mode only when you got stuck in the standard mode): 

	python finisherSC.py destinedFolder mummerPath -f True 
[It means that you want to use the fast mode ]

b)Parallel mode. \[-par numberOfThreads\] You can now use -par 20 to run finisherSC on 20 threads. The command is 
        
	python finisherSC.py -par 20 destinedFolder mummerPath


c)Break down large contig file \[-l True\]. If your contigs.fasta is too big, then you want to use this option to break it down for alignment. 

	python finisherSC.py -l True destinedFolder mummerPath


d)Pick from previous jobs \[-p pickupFilename\] (e.g. you have got improved2.fasta but the cluster node timeout and you want to start from improved2.fasta instead of from scratch)

	python finisherSC.py destinedFolder mummerPath  -p improved2.fasta

e)Mapping between old contigs and new contigs \[-o referenceName_QueryName \](e.g. you want to know how the contigs from contigs.fasta are mapped to improved3.fasta )

	python finisherSC.py destinedFolder mummerPath -o contigs.fasta_improved3.fasta
[It will then output the alignment of improved3.fasta against contigs.fasta. The output will be shown in the terminal and in mappingResults.txt at the destinedFolder]

f)Help \[-h\]. You can always use -h to get the usage suggestion. 

	python finisherSC.py -h


g)There is an experimental improvement based on repeat phasing and the further improved file is improved4.fasta. An illustration of that part is given in Fig. 7 of http://arxiv.org/abs/1402.6971 . To experiment with that, first run finisherSC.py as before. After that, you can issue the following command:

	python experimental/xPhaser.py destinedFolder mummerPath

h) There is also an experimental improvement resolving tandem repeat. The further improved file is tademResolved.fasta. To experiment with that, first run finisherSC.py as before. After that, you can issue the following command:

	python experimental/tSolver.py destinedFolder mummerPath


## Data files for experiments ##
We are using the benchmark datasets available online in PacBio DevNet. And we can run FinisherSC to process those data to completion. For (1, 2, 3), we use -par 20 option and run it on a server computer and we get them finished in a couple of hours to a day. For (4), we run without any options on a laptop and get it finished in a couple of minutes. 
  
1. Drosophila : https://github.com/PacificBiosciences/DevNet/wiki/Drosophila-sequence-and-assembly 
2. Saccharomyces cerevisiae W303 : https://github.com/PacificBiosciences/DevNet/wiki/Saccharomyces-cerevisiae-W303-Assembly-Contigs
3. C. elegans : https://github.com/PacificBiosciences/DevNet/wiki/C.-elegans-data-set
4. E. coli, M. ruber and P. heparinus : http://files.pacb.com/software/hgap/index.html

The data are supported by the original authors. But we provide an example on how to download and transform data types so it is easier for users to validate FinisherSC on benchmark data . 

1. Download data with links specified in file_list

        for f in `cat file_list`; do wget --force-directories $f; done 
	
2. Download the DEXTRACTOR

        git clone https://github.com/thegenemyers/DEXTRACTOR.git
	
3. Transform .bax.h5 files to .fasta file
	
        find . -name '*.bax.h5' | xargs DEXTRACTOR/dextract  > raw_reads.fasta



## Remarks ##
1. After passing through the finishingTool, it is advised to use a polishing tool(e.g. Quiver) to polish the improved contigs( i.e. improved3.fasta )
2. If you have problem running scipy.weave, try to remove ~/.python27_compiled 

## Annotation of FinisherSC ##
Main functional components of FinisherSC : 

> finisherSC.py : Starting point of FinisherSC
	
> nonRedundantResolver.py : Filter completely embedded contigs (aka Subroutine 1 in the paper's supplementary material)
	
> overlapResolver.py : Utilize existing overlaps of contigs to resolve repeats and merge contigs (aka Subroutine 2,3 in the paper's supplementary material)
	
> gapFiller.py : Utilize raw reads data to fill gaps(aka Subroutine 4, 5 in the paper's supplementary material)
	
> twoRepeatOneBridgeSolver.py : Resolve repeats with two copies where only one copy is spanned by some reads (aka Subroutine 6,7 in the paper's supplementary material)

Helper Functions :

> IORobot.py : various functions that support reading and writing of data

> alignerRobot.py : library containing ways to parse MUMmer results

> houseKeeper.py : various functions to support house keeping of data

> graphLib.py : library containing various graph string graph operations 

> debugging.py : a debugging point for hacking into various functions 

> unittester.py : integration test package

> viewer.py : create dot plots 


An older version of the tool (no longer supported) : 

> finisher.py : an independent and functioning(but no longer supported) version based on a greedy algorithm

Relevant third parties software: 

> fasta-splitter.pl : a tool to break FASTA file into smaller trucks in commandline

> MUMmer : an alignment tool to perform the mapping of reads and contigs

> gepard-1.30 : a good tool to create dotplots

Experimental folder: 

> There are various experimental functions there. The two important ones are tandem repeat resolver(tSolver.py) and approximate repeat phaser(xPhaser.py). 

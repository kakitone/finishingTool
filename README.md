finishingTool
=============
If you encounter any issue, please feel free to directly contact me at kklam@eecs.berkeley.edu
My name is  Ka-Kit Lam. 

Our paper is located at http://biorxiv.org/content/early/2014/10/14/010215 

## Command to run the tool ##

Here is the command to run the tool:

	python finisherSC.py destinedFolder mummerPath
	
## Input and output files ##

Note that it is assumed that you have "contigs.fasta" , "raw_reads.fasta" in the destinedFolder and mummerPath is your path to Mummer. Therefore, please rename the files in the destinedFolder if needed. For "finisherSC.py", the key output file is **improved3.fasta** in destinedFolder. 


## A sample test run ##
Here is an example run with the pre-installed data files and softwares(e.g. Mummer, Gepard) [Please go to https://www.dropbox.com/sh/xjpt8xf5g1xf0ek/bGmZvt9Zfd to Download the EColi-Data set as test cases coz it is too big for GitHub].

1. You may run the test cases as python finisherSC.py EcoliTestRun/ MUMmer3.23/
2. After the processing, take the file whose name is **improved3.fasta** in destinedFolder and this is your output file.
3. If you want to visualize the output against the reference, you can use software like Gepard.

## Options[Supplementary usage] ##
a)Fast mode (e.g. for large and complex genomes , you may want to use this mode to get results faster; However, there will be a bit of computation vs quality trade-off if you want to use this mode; We suggest that you use this mode only when you got stuck in the standard mode): 

	python finisherSC.py destinedFolder MUMmer3.23 -f True 
[It means that you want to use the fast mode ]

b)Pick from previous jobs (e.g. you have got improved2.fasta but the cluster node timeout and you want to start from improved2.fasta instead of from scratch)

	python finisherSC.py destinedFolder MUMmer3.23  -p improved2.fasta

c)Mapping between old contigs and new contigs(e.g. you want to know how the contigs from contigs.fasta are mapped to improved3.fasta )

	python finisherSC.py destinedFolder MUMmer3.23 -o contigs.fasta_improved3.fasta
[It will then output the alignment of improved3.fasta against contigs.fasta. The output will be shown in the terminal and in mappingResults.txt at the destinedFolder]

d)Help. You can always use -h to get the usage suggestion. 

	python finisherSC.py -h


## Remarks ##
1. After passing through the finishingTool, it is advised to use a polishing tool(e.g. Quiver) to polish the improved contigs( i.e. improved3.fasta )
2. If you have problem running scipy.weave, try to remove ~/.python27_compiled 
3. There is an experimental improvement based on repeat phasing and the further improved file is improved4.fasta. An illustration of that part is given in Fig. 7 of http://arxiv.org/abs/1402.6971 . To experiment with that, first run finisherSC.py as before. After that, you can issue the following command:

	python experimental/newPhasing.py destinedFolder mummerPath




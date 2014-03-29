Gepard  (GEnome PAir - Rapid Dotter)

Version: 1.30
Freeware 

Author:   Jan Krumsiek
Email:    krumsiek@gmx.de
Web:      http://mips.gsf.de/services/analysis/gepard


Reference
=========
Krumsiek J, Arnold R, Rattei T
Gepard: A rapid and sensitive tool for creating dotplots on genome scale.
Bioinformatics 2007; 23(8): 1026-8. PMID: 17309896


Contents
========
1. System requirements
2. Tutorial
3. Startup scripts
4. Command line tool


1. System requirements
======================
Gepard requires the Java Runtime Environment Version 5.0 or higher which is 
available at http://www.java.com/download/


2. Tutorial
===========
This download archive contains a file called "tutorial.html" which is a quick-
start guide for your first steps with the program.


3. Startup scripts
==================
There are three startup scripts for both Windows and Un*x systems:

1. Normal memory (256MB)      Win: gepard.bat         Un*x: gepard.sh
2. High memory (512MB)        Win: gepard_himem.bat   Un*x: gepard_himem.sh
3. Very high memory (1024MB)  Win: gepard_vhimem.bat  Un*x: gepard_vhimem.sh


The following table shows the approximate sequence sizes (self-plot) you can 
handle with each memory setting (including Gepard suffix array calculation by 
Gepard):

 256 MB   ~10 million base pairs
 512 MB   ~20 million base pairs
1024 MB   ~40 million base pairs

If none of the startup scripts works for you, you need to call Java manually
from the command line of your system. Simply copy-paste the "java ..." line
from the respective script file into your console.


4. Command line tool
====================
Since version 1.20 Gepard contains a comand line functionality. In its current 
version it only supports OFFLINE dotplots which means that you have to online
access to the genomes in the PEDANT database from the command line.

To start the comand line tool of Gepard use "gepardcmd.sh" on Un*x systems and
"gepardcmd.bat" on Windows systems. If the command line startup script does 
not work on your system you need to start the program manually as outlined
above in section 3.


import time
import argparse

import nonRedundantResolver
import overlapResolver
import gapFiller
import twoRepeatOneBridgeSolver
import houseKeeper

###################################################### Starting point
def mainFlow(folderName , mummerLink, pickupname, mapcontigsname):      
    print "Go Bears! ! !" 
    
    print "pickupname, mapcontigsname", pickupname, mapcontigsname
    
    
    if not pickupname in ["noEmbed.fasta", "improved.fasta", "improved2.fasta"]:
        nonRedundantResolver.removeEmbedded(folderName , mummerLink)
     
    if not pickupname in ["improved.fasta", "improved2.fasta"]:
        overlapResolver.fetchSuccessor(folderName , mummerLink)
        overlapResolver.formSeqGraph(folderName , mummerLink)
    
    if not pickupname in ["improved2.fasta"]:
        gapFiller.fillGap(folderName , mummerLink)
    
    twoRepeatOneBridgeSolver.xPhased(folderName , mummerLink)
    
    # ECReduction(folderName , mummerLink )
    # compareWithReference(folderName , mummerLink)
    
    if mapcontigsname != None:
        houseKeeper.performMapping(folderName, mummerLink, mapcontigsname)
        
    print "<3 Do cool things that matter <3"

# folderName = "S_cerivisea/"
# mummerLink = "MUMmer3.23/"

t0 = time.time()

parser = argparse.ArgumentParser(description='FinisherSC : a repeat-aware tool to upgrade de-novo assembly with long reads')
parser.add_argument('folderName')
parser.add_argument('mummerLink')
parser.add_argument('-p', '--pickup', help='Picks up existing work (input is noEmbed.fasta, improved.fasta or improved2.fasta)', required=False)
parser.add_argument('-o', '--mapcontigs', help='Maps new contigs to old contigs(input is of the format of contigs.fasta_improved3.fasta which means improved3.fasta will be mapped back to contigs.fasta; Output can be found in mappingResults.txt in the destinedFolder;)', required=False)
parser.add_argument('-f', '--fast', help= 'Fast aligns contigs (input is True)', required=False)
parser.add_argument('-par', '--parallel', help= 'Fast aligns contigs (input is maximum number of threads)', required=False)
parser.add_argument('-l', '--large', help= 'Large number of contigs/large size of contigs (input is True)', required=False)


args = vars(parser.parse_args())

print "args", args
pathExists, newFolderName, newMummerLink = houseKeeper.checkingPath(args['folderName'], args['mummerLink'])

if args['fast'] == "True":
    houseKeeper.globalFast = True
else:
    houseKeeper.globalFast = False

if args['parallel'] != None:
    houseKeeper.globalParallel = int(args['parallel'])
else:
    houseKeeper.globalParallel = 1


if args['large'] == "True":
    houseKeeper.globalLarge = True
else:
    houseKeeper.globalLarge = False


if pathExists:
    mainFlow(newFolderName, newMummerLink, args['pickup'], args['mapcontigs'])
else:
    print "Sorry. The above folders or files are missing. If you continue to have problems, please contact me(Ka-Kit Lam) at kklam@eecs.berkeley.edu"

print  "Time", time.time() - t0

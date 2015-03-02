import abunSplitter

from finisherSCCoreLib import houseKeeper
import time
import argparse
import abunHouseKeeper

t0 = time.time()

parser = argparse.ArgumentParser(description='aSplitter')
parser.add_argument('folderName')
parser.add_argument('mummerLink')

parser.add_argument('-f', '--fast', help= 'Fast aligns contigs (input is True)', required=False)
parser.add_argument('-par', '--parallel', help= 'Fast aligns contigs (input is maximum number of threads)', required=False)
parser.add_argument('-l', '--large', help= 'Large number of contigs/large size of contigs (input is True)', required=False)


parser.add_argument('-rp', '--replace', help= 'Input files to aSplitter(e.g. noEmbed.fasta, improved.fasta, improved2.fasta or improved3.fasta)', required=False)
parser.add_argument('-ar', '--avoidrefine', help= 'Avoid refined abundance estimation (input is True)', required=False)
parser.add_argument('-rs', '--readsearch', help= 'Number of linking reads across a gap  (input is number of such linking reads/2)', required=False)
parser.add_argument('-rd', '--RRDisable', help= 'Whether one should disable Read to Read overlap check (input is True)', required=False)


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


if args['avoidrefine'] == "True":
    abunHouseKeeper.abunGlobalAvoidrefine = True
else:
    abunHouseKeeper.abunGlobalAvoidrefine = False


if args['readsearch'] != None:
    abunHouseKeeper.abunGlobalReadSearchDepth = int(args['readsearch']) 
else:
    abunHouseKeeper.abunGlobalReadSearchDepth = 1


if args['replace'] != None : 
    replacedName = args['replace']
    abunHouseKeeper.replaceFiles( newFolderName, replacedName) 

if args['RRDisable'] == "True":
    abunHouseKeeper.abunGlobalRRDisable = True
else:
    abunHouseKeeper.abunGlobalRRDisable = False



if pathExists:
    abunSplitter.mainFlow(newFolderName, newMummerLink)
else:
    print "Sorry. The above folders or files are missing. If you continue to have problems, please contact me(Ka-Kit Lam) at kklam@eecs.berkeley.edu"

print  "Time", time.time() - t0

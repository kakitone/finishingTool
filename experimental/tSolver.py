import tandemRepeatSolver

import time
import argparse

t0 = time.time()

parser = argparse.ArgumentParser(description='tSolver')
parser.add_argument('folderName')
parser.add_argument('mummerLink')

args = vars(parser.parse_args())
print "args", args
tandemRepeatSolver.mainFlowForTandemResolve(args['folderName'] , args['mummerLink'])
print  "Time", time.time() - t0


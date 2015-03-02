import matplotlib.pyplot as plt
from finisherSCCoreLib import IORobot

lenDic = {}
coverageDic = {}

lenDic = IORobot.obtainLength("/Users/kakitlam/", "abun.fasta")

f = open("/Users/kakitlam/Documents/abundata", 'r')
tmp = f.readline()

while len(tmp) > 0:
	if len(tmp) > 10:
		myitem = tmp[0:-1].split()
		coverageDic[myitem[0]] = float(myitem[1])
	tmp = f.readline()

f.close()

myList = []
baseCt = {}

for eachitem in lenDic:
	myList.append(lenDic[eachitem]*coverageDic[eachitem])
	baseCt[eachitem] = lenDic[eachitem]*coverageDic[eachitem]


for eachitem in lenDic :
	print eachitem,  baseCt[eachitem]



for eachitem in lenDic :
	print eachitem, lenDic[eachitem]


for eachitem in lenDic :
        print eachitem, coverageDic[eachitem]

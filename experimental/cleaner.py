import loggingIndel
#import shortToLongMapper
import common
import numpy as np
import math

from scipy import weave

def reverseStrand(originalGen):
    #print originalGen
    reverseGenome = np.zeros(len(originalGen), dtype = np.int8)
    round = len(originalGen)
    # A, C, G, T = 1, 2, 3, 4 
    code =\
    """
    
    int i; 

    for (i =1 ;i <round ; i++){
        if (originalGen[i] == 1){
            reverseGenome[i] = 4; 
        }
        else if (originalGen[i] == 2){
            reverseGenome[i] = 3; 
        }else if (originalGen[i] == 3){
            reverseGenome[i] = 2; 
        }else if (originalGen[i] == 4){
            reverseGenome[i] = 1; 
        }
    }

    """
    
    weave.inline(code, ['reverseGenome', 'originalGen', 'round'])
    return reverseGenome    



def reverseString(str):
    newstr = []
    for index in range(len(str)):
        newstr.append(str[len(str)-1-index])
    
    return newstr

def printSeq(str2):
    printStr = ""

    for eachchar in str2:
        if eachchar != 0:
            printStr = printStr + str(eachchar)
        elif eachchar == 0 :
            printStr = printStr + " "
    print printStr + "\n"

def transformToReverseStrand(readList):
    reverseStrandList = []
    # ACGT = 1234
    for eachread in readList:
        tempread= reverseStrand(eachread)                
        reverseStrandList.append(tempread)
                
    return reverseStrandList 

def meetRequirement(score, alignedSeq1, alignedSeq2,eachshortread, eachlongread,istart, jstart, iend, jend, threshold, liid, overhang ):
    check = True
    #print "score , ", score 
    #print "istart, jstart, iend, jend", istart, jstart, iend, jend
    #print "alignedSeq1", alignedSeq1
    #print "alignedSeq2", alignedSeq2
    if ( iend - istart ) < liid :

        check = False 
    
    if score < (iend - istart)* threshold:
       # print iend, istart, jend, jstart,  threshold, score
        check = False

    if min(len(eachshortread) - iend, len(eachlongread) - jend)  > overhang or min(istart, jstart) > overhang :
        check  =False 

    #print check 
    
    return check 

def SWAlignment(seq1 , seq2, parameterRobot):
    score  = 0 

    wts = parameterRobot.editsub
    wti = parameterRobot.editins 
    wtd = parameterRobot.editdel
    wtm=  parameterRobot.editmatch
    
    #print "wts, wti, wtd, wtm",wts, wti, wtd, wtm
    
    m = len(seq1) + 1
    n = len(seq2) +1
    
    #H = np.zeros(m*n, dtype = np.float64).reshape(m,n)
    #B = np.zeros(m*n, dtype = np.float64).reshape(m,n)
    H = np.zeros([m,n], dtype = np.float64)
    B = np.zeros([m,n], dtype = np.float64)

    # Assign weights 
    for i in range(m):
        H[i][0] = 0
        B[i][0] = 4
    for j in range(n):
        H[0][j]  =0
        B[0][j] = 4 
        
        
        
    seq1NP = np.zeros(m-1, dtype = np.float64)
    seq2NP = np.zeros(n-1, dtype = np.float64)
    
    for i in range(m-1):
        seq1NP[i] = seq1[i]
    for j in range(n-1) :
        seq2NP[j] = seq2[j]
        
 #   print "seq1NP, seq2NP"
 #   print seq1NP, seq2NP
 #   print "seq1, seq2"
 #   print seq1, seq2
    
    code =\
        """
        int i; 
        int j ;
        double w;

        for (i =1 ;i <m ; i++){
            for (j=1; j<n; j++){
                if (seq1NP[i-1] == seq2NP[j-1]){
                    w = wtm ;
                }
                else{
                    w= wts ;
                }
                
                    H2(i,j) = 0 ;
                    
                    if (H2(i,j) < H2(i-1,j-1) + w) {
                        H2(i,j) = H2(i-1,j-1) + w ; 
                    }

                    if (H2(i,j) < H2(i-1,j)+wtd ) {
                        H2(i,j) = H2(i-1,j)+wtd ;
                    }

                    if  (H2(i,j) < H2(i,j-1) + wti){
                        H2(i,j) = H2(i,j-1) + wti;
                    }

                    if (H2(i-1,j-1) + w == H2(i,j)){
                        B2(i,j) = 1 ; 
                    }
                    else if (H2(i-1,j)+wtd == H2(i,j)){
                        B2(i,j) = 2;
                    }
                    
                    else if  (H2(i,j-1) + wti == H2(i,j)){
                        B2(i,j) = 3;
                    }
                    else if (0 == H2(i,j)) {
                        B2(i,j) = 4 ;
                    }
            }
        }

        """
    
    weave.inline(code, ['H','B','m', 'n', 'wtm','wts','wti','wtd', 'seq1NP', 'seq2NP'])
    
    
#    for i in range(1, m):
#        for j in range(1, n):
#            if seq1[i-1] == seq2[j-1]:
#                w = wtm 
#            else:
#                w= wts
            
#            H[i][j] = max(H[i-1][j-1] + w, H[i-1][j]+wtd, H[i][j-1] + wti, 0)
#            if H[i-1][j-1] + w == H[i][j] :
#                B[i][j] = 1 
#            elif H[i-1][j]+wtd == H[i][j] :
#                B[i][j] = 2
#            elif  H[i][j-1] + wti == H[i][j] :
#                B[i][j] = 3
#            elif 0 == H[i][j]:
#                B[i][j] = 4 

    # Backtrack 
    alignedSeq1 = []
    alignedSeq2 = []
    
    bestindex = np.argmax(H)

 
    
    endi = bestindex / n
    endj = bestindex%n
    
    #print endi , endj
    score = H[endi][endj]
    
   # print "score", score 
   # print "B", B[1]
   # print "H", H
    
    tempi, tempj = endi , endj
    while (B[tempi][tempj] != 4):
        if B[tempi][tempj] == 1:
            alignedSeq1.append(seq1[tempi -1 ])
            alignedSeq2.append(seq2[tempj-1 ]) 
            tempi = tempi -1 
            tempj = tempj -1           

        elif B[tempi][tempj] == 2 :
            alignedSeq1.append(seq1[tempi-1 ])
            alignedSeq2.append(0) 
            tempi = tempi - 1

        elif B[tempi][tempj] == 3 :
            alignedSeq1.append(0)
            alignedSeq2.append(seq2[tempj-1 ])
            tempj = tempj -1

    
    starti, startj = tempi , tempj 
    
    returnalignedSeq1 = reverseString(alignedSeq1)
    returnalignedSeq2 = reverseString(alignedSeq2)
    
    #printSeq(returnalignedSeq1)
    #printSeq(returnalignedSeq2)
    return score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj


def SWAlignmentFixRef(seq1, seq2, parameterRobot):
    score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = 0,0,0,0,0,0,0


    wts = parameterRobot.editsub
    wti = parameterRobot.editins 
    wtd = parameterRobot.editdel
    wtm=  parameterRobot.editmatch
    
    
    
    m = len(seq1) + 1
    n = len(seq2) +1
    
    H = np.zeros(m*n, dtype = np.int32).reshape(m,n)
    B = np.zeros(m*n, dtype = np.int32).reshape(m,n)
    # Assign weights 
    for i in range(m):
        H[i][0] = wti*i
        B[i][0] = 4
    for j in range(n):
        H[0][j]  =0
        B[0][j] = 4 
        
    for i in range(1, m):
        for j in range(1, n):
            if seq1[i-1] == seq2[j-1]:
                w = wtm 
            else:
                w= wts
            
            H[i][j] = max (H[i-1][j-1] + w, H[i-1][j]+wtd, H[i][j-1] + wti, H[i][0])
            if H[i-1][j-1] + w == H[i][j] :
                B[i][j] = 1 
            elif H[i-1][j]+wtd == H[i][j] :
                B[i][j] = 2
            elif  H[i][j-1] + wti == H[i][j] :
                B[i][j] = 3
            elif H[i][0] == H[i][j]:
                B[i][j] = 4 

    # Backtrack 
    alignedSeq1 = []
    alignedSeq2 = []
    
    bestindex = np.argmax(H[:][m-1])

 
    
    endi = m -1 
    endj = bestindex
    
    #print endi , endj
    score = H[endi][endj]

    tempi, tempj = endi , endj
    while (B[tempi][tempj] != 4):
        if B[tempi][tempj] == 1:
            alignedSeq1.append(seq1[tempi -1 ])
            alignedSeq2.append(seq2[tempj-1 ]) 
            tempi = tempi -1 
            tempj = tempj -1           

        elif B[tempi][tempj] == 2 :
            alignedSeq1.append(seq1[tempi-1 ])
            alignedSeq2.append(0) 
            tempi = tempi - 1

        elif B[tempi][tempj] == 3 :
            alignedSeq1.append(0)
            alignedSeq2.append(seq2[tempj-1 ])
            tempj = tempj -1

    
    
    for dummy in range(tempi-1,-1,-1):
        alignedSeq1.append(seq1[dummy])
        alignedSeq2.append(0)
        
    starti, startj = tempi , tempj 
    
    returnalignedSeq1 = reverseString(alignedSeq1)
    returnalignedSeq2 = reverseString(alignedSeq2)
    
    #printSeq(returnalignedSeq1)
    #printSeq(returnalignedSeq2)
    
    return score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj


def SWAlignmentFixRefMP(seq1,seq2,  parameterRobot):
    score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj ,numPath= 0,0,0,0,0,0,0, 0

    wts = parameterRobot.editsub
    wti = parameterRobot.editins 
    wtd = parameterRobot.editdel
    wtm=  parameterRobot.editmatch

    m = len(seq1) + 1
    n = len(seq2) + 1
    
    H = np.zeros(m*n, dtype = np.int32).reshape(m,n)
    B = np.zeros(m*n*4, dtype = np.int32).reshape(m,n,4)

    seq1NP = np.zeros(m-1, dtype = np.float64)
    seq2NP = np.zeros(n-1, dtype = np.float64)
    
    for i in range(m-1):
        seq1NP[i] = seq1[i]
    for j in range(n-1) :
        seq2NP[j] = seq2[j]
        
    # Assign weights 
    for i in range(m):
        H[i][0] = wti*i
        B[i][0][3]= 1
        
    for j in range(n):
        H[0][j]  =0
        B[0][j][3] = 1
    
    
    code =\
    """
    #include <math.h>   
    
    int i; 
    int j ;
    double w;
    int tmp ;

    for (i =1 ;i <m ; i++){
        for (j=1; j<n; j++){
            if (seq1NP[i-1] == seq2NP[j-1]){
                w = wtm ;
            }
            else{
                w= wts ;
            }
            
                H2(i,j) = H2(i,0) ;
                
                if (H2(i,j) < H2(i-1,j-1) + w) {
                    H2(i,j) = H2(i-1,j-1) + w ; 
                }

                if (H2(i,j) < H2(i-1,j)+wtd ) {
                    H2(i,j) = H2(i-1,j)+wtd ;
                }

                if  (H2(i,j) < H2(i,j-1) + wti){
                    H2(i,j) = H2(i,j-1) + wti;
                }

                if (H2(i-1,j-1) + w == H2(i,j)){
                    B3(i,j,0) =1 ;
                }
                
                if (H2(i-1,j)+wtd == H2(i,j)){
                    B3(i,j,1) =1 ;
                }
                
                if  (H2(i,j-1) + wti == H2(i,j)){
                    B3(i,j,2) =1 ;
                }
                
                if (H2(i,0) == H2(i,j)) {
                    B3(i,j,3) =1 ;
                }
        }
    }

    """
    
    weave.inline(code, ['H','B','m', 'n', 'wtm','wts','wti','wtd', 'seq1NP', 'seq2NP'])


#    for i in range(1, m):
#        for j in range(1, n):
#            if seq1[i-1] == seq2[j-1]:
#               w = wtm 
#            else:
#                w= wts
            
#            H[i][j] = max (H[i-1][j-1] + w, H[i-1][j]+wtd, H[i][j-1] + wti, H[i][0])
#            
#            if H[i-1][j-1] + w == H[i][j] :
#                B[i][j].append(1) 
#                
#            if H[i-1][j]+wtd == H[i][j] :
#                B[i][j].append(2)
#                
#            if  H[i][j-1] + wti == H[i][j] :
#                B[i][j].append(3)
#                
#            if H[i][0] == H[i][j]:
#                B[i][j].append(4) 

    # Backtrack 
    alignedSeq1 = []
    alignedSeq2 = []
    
    bestindex = np.argmax(H[:][m-1])

 
    
    endi = m -1 
    endj = bestindex
    
    #print endi , endj
    score = H[endi][endj]
    #print B
    
    tempi, tempj = endi , endj
    
    
    currentIndex = np.argmax(B[tempi][tempj][:]) +1 
    #print currentIndex
    #print currentIndex ,B[tempi][tempj]
    
    while (currentIndex != 4):
        if currentIndex == 1:
            alignedSeq1.append(seq1[tempi -1 ])
            alignedSeq2.append(seq2[tempj-1 ]) 
            tempi = tempi -1 
            tempj = tempj -1           

        elif currentIndex == 2 :
            #print tempi, seq1
            alignedSeq1.append(seq1[tempi-1 ])
            alignedSeq2.append(0) 
            tempi = tempi - 1

        elif currentIndex == 3 :
            alignedSeq1.append(0)
            alignedSeq2.append(seq2[tempj-1 ])
            tempj = tempj -1
            
        currentIndex = np.argmax(B[tempi][tempj][:]) +1 
        #print currentIndex
    
    
    for dummy in range(tempi-1,-1,-1):
        alignedSeq1.append(seq1[dummy])
        alignedSeq2.append(0)
        
    starti, startj = tempi , tempj 
    
    returnalignedSeq1 = reverseString(alignedSeq1)
    returnalignedSeq2 = reverseString(alignedSeq2)
    
    #printSeq(returnalignedSeq1)
    #printSeq(returnalignedSeq2)
    
    # Count Path
    numPath = 0
    tempi, tempj = endi , endj
    stack = [[tempi, tempj]]
    
    for j in range(n):
        if B[1][j][3] == 1 and np.sum(B[1][j][:]) > 1 :
            B[1][j][:] = 0
            B[1][j][3] = 1
    
    while (len(stack) > 0):
        targetIndex = stack.pop()
        tempi, tempj = targetIndex[0], targetIndex[1]
        tmpList = []
        for i in range(4):
            if B[tempi][tempj][i] > 0:
                tmpList.append(i+1)
        #print "len(B[tempi][tempj])",len(B[tempi][tempj]),B[tempi][tempj]
        for ptrDir in tmpList:
            #print tempi, tempj ,ptrDir
            if ptrDir == 1:
                stack.append([tempi-1, tempj-1])

            elif ptrDir== 2 :
                stack.append([tempi-1, tempj])
    
            elif ptrDir== 3 :
                stack.append([tempi, tempj-1])
                
            elif ptrDir == 4:
                numPath += 1
        
    
    #print numPath
    return score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath



def SWAlignmentFixRefMPQS(seq1,seq2,  parameterRobot):
    score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj ,numPath= 0,0,0,0,0,0,0, 0

    p  = parameterRobot.p 
    wts = math.log(p/3*p)
    wti = math.log(p/3) 
    wtd = math.log(p)
    wtm=  2*math.log(1-p)



    m = len(seq1) + 1
    n = len(seq2) + 1
    
    H = np.zeros(m*n, dtype = np.float64).reshape(m,n)
    B = np.zeros(m*n*4, dtype = np.float64).reshape(m,n,4)

    seq1NP = np.zeros(m-1, dtype = np.float64)
    seq2NP = np.zeros(n-1, dtype = np.float64)
    
    for i in range(m-1):
        seq1NP[i] = seq1[i]
    for j in range(n-1) :
        seq2NP[j] = seq2[j]
        
    # Assign weights 
    for i in range(m):
        H[i][0] = wtd*i
        B[i][0][3]= 1
        
    for j in range(n):
        H[0][j]  =0
        B[0][j][3] = 1
    
    
    code =\
    """
    #include <math.h>   
    
    int i; 
    int j ;
    double w;
    int tmp ;

    for (i =1 ;i <m ; i++){
        for (j=1; j<n; j++){
            if (seq1NP[i-1] == seq2NP[j-1]){
                w = wtm ;
            }
            else{
                w= wts ;
            }
            
                H2(i,j) = H2(i,0) ;
                
                if (H2(i,j) < H2(i-1,j-1) + w) {
                    H2(i,j) = H2(i-1,j-1) + w ; 
                }

                if (H2(i,j) < H2(i-1,j)+wtd ) {
                    H2(i,j) = H2(i-1,j)+wtd ;
                }

                if  (H2(i,j) < H2(i,j-1) + wti){
                    H2(i,j) = H2(i,j-1) + wti;
                }

                if (H2(i-1,j-1) + w == H2(i,j)){
                    B3(i,j,0) =1 ;
                }
                
                if (H2(i-1,j)+wtd == H2(i,j)){
                    B3(i,j,1) =1 ;
                }
                
                if  (H2(i,j-1) + wti == H2(i,j)){
                    B3(i,j,2) =1 ;
                }
                
                if (H2(i,0) == H2(i,j)) {
                    B3(i,j,3) =1 ;
                }
        }
    }

    """
    
    weave.inline(code, ['H','B','m', 'n', 'wtm','wts','wti','wtd', 'seq1NP', 'seq2NP'])


    # Backtrack 
    alignedSeq1 = []
    alignedSeq2 = []
    
    bestindex = np.argmax(H[:][m-1])

    endi = m -1 
    endj = bestindex
    
    #print endi , endj
    score = H[endi][endj]
    #print "score", score
    
    tempi, tempj = endi , endj
    
    
    currentIndex = np.argmax(B[tempi][tempj][:]) +1 
    #print currentIndex
    #print currentIndex ,B[tempi][tempj]
    
    while (currentIndex != 4):
        if currentIndex == 1:
            alignedSeq1.append(seq1[tempi -1 ])
            alignedSeq2.append(seq2[tempj-1 ]) 
            tempi = tempi -1 
            tempj = tempj -1           

        elif currentIndex == 2 :
            #print tempi, seq1
            alignedSeq1.append(seq1[tempi-1 ])
            alignedSeq2.append(0) 
            tempi = tempi - 1

        elif currentIndex == 3 :
            alignedSeq1.append(0)
            alignedSeq2.append(seq2[tempj-1 ])
            tempj = tempj -1
            
        currentIndex = np.argmax(B[tempi][tempj][:]) +1 
        #print currentIndex
    
    
    for dummy in range(tempi-1,-1,-1):
        alignedSeq1.append(seq1[dummy])
        alignedSeq2.append(0)
        
    starti, startj = tempi , tempj 
    
    returnalignedSeq1 = reverseString(alignedSeq1)
    returnalignedSeq2 = reverseString(alignedSeq2)
    
    #printSeq(returnalignedSeq1)
    #printSeq(returnalignedSeq2)
    
    # Count Path
    numPath = 0
    tempi, tempj = endi , endj
    stack = [[tempi, tempj]]
    
    for j in range(n):
        if B[1][j][3] == 1 and np.sum(B[1][j][:]) > 1.5 :
            B[1][j][:] = 0
            B[1][j][3] = 1
    
    while (len(stack) > 0):
        targetIndex = stack.pop()
        tempi, tempj = targetIndex[0], targetIndex[1]
        tmpList = []
        for i in range(4):
            if B[tempi][tempj][i] > 0:
                tmpList.append(i+1)
        #print "len(B[tempi][tempj])",len(B[tempi][tempj]),B[tempi][tempj]
        #print "tmpList", tmpList
        for ptrDir in tmpList:
            #print tempi, tempj ,ptrDir
            if ptrDir == 1:
                stack.append([tempi-1, tempj-1])

            elif ptrDir== 2 :
                stack.append([tempi-1, tempj])
    
            elif ptrDir== 3 :
                stack.append([tempi, tempj-1])
                
            elif ptrDir == 4:
                numPath += 1
        
    
    
    score = score + math.log(numPath)
    #print "numPath", numPath
    return score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath



def SWAlignmentFixRefMPQS2(seq1,seq2,  parameterRobot):
    score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj ,numPath= 0,0,0,0,0,0,0, 0

    p  = parameterRobot.p 
    wts = math.log(p/3*p)
    wti = math.log(p/3) 
    wtd = math.log(p)
    wtm=  2*math.log(1-p)


    #print wts, wti, wtd, wtm
    m = len(seq1) + 1
    n = len(seq2) + 1
    
    H = np.zeros(m*n, dtype = np.double).reshape(m,n)
    #B = np.zeros(m*n, dtype = np.int32).reshape(m,n)
    B = [[ [] for j in range(n)]for i in range(m)]
    
    # Assign weights 
    for i in range(m):
        H[i][0] = wtd*i
        B[i][0].append(4)
        
    for j in range(n):
        H[0][j]  =0
        B[0][j].append(4) 
    
    B[0][0].pop()

    
    for i in range(1, m):
        for j in range(1, n):
            if seq1[i-1] == seq2[j-1]:
                w = wtm 
            else:
                w= wts
            
            H[i][j] = max (H[i-1][j-1] + w, H[i-1][j]+wtd, H[i][j-1] + wti, H[i][0])
            
            if H[i-1][j-1] + w == H[i][j] :
                B[i][j].append(1) 
                
            if H[i-1][j]+wtd == H[i][j] :
                B[i][j].append(2)
                
            if  H[i][j-1] + wti == H[i][j] :
                B[i][j].append(3)
                
            if H[i][0] == H[i][j]:
                B[i][j].append(4) 

    # Backtrack 
    alignedSeq1 = []
    alignedSeq2 = []
    
    bestindex = np.argmax(H[:][m-1])

 
    
    endi = m -1 
    endj = bestindex
    
    #print endi , endj
    score = H[endi][endj]
    print "score" , score
    
    tempi, tempj = endi , endj
    while (B[tempi][tempj][0] != 4):
        if B[tempi][tempj][0] == 1:
            alignedSeq1.append(seq1[tempi -1 ])
            alignedSeq2.append(seq2[tempj-1 ]) 
            tempi = tempi -1 
            tempj = tempj -1           

        elif B[tempi][tempj][0] == 2 :
            alignedSeq1.append(seq1[tempi-1 ])
            alignedSeq2.append(0) 
            tempi = tempi - 1

        elif B[tempi][tempj][0] == 3 :
            alignedSeq1.append(0)
            alignedSeq2.append(seq2[tempj-1 ])
            tempj = tempj -1

    
    
    for dummy in range(tempi-1,-1,-1):
        alignedSeq1.append(seq1[dummy])
        alignedSeq2.append(0)
        
    starti, startj = tempi , tempj 
    
    returnalignedSeq1 = reverseString(alignedSeq1)
    returnalignedSeq2 = reverseString(alignedSeq2)
    
    #printSeq(returnalignedSeq1)
    #printSeq(returnalignedSeq2)
    
    # Count Path
    numPath = 0
    tempi, tempj = endi , endj
    stack = [[tempi, tempj]]
    
    for j in range(n):
        if 4 in B[1][j] and len(B[1][j]) >1 :
            B[1][j] = [4]
    
    while (len(stack) > 0):
        targetIndex = stack.pop()
        tempi, tempj = targetIndex[0], targetIndex[1]
        
        #print "len(B[tempi][tempj])",len(B[tempi][tempj]),B[tempi][tempj]
        for ptrDir in B[tempi][tempj]:
            #print tempi, tempj ,ptrDir
            if ptrDir == 1:
                stack.append([tempi-1, tempj-1])

            elif ptrDir== 2 :
                stack.append([tempi-1, tempj])
    
            elif ptrDir== 3 :
                stack.append([tempi, tempj-1])
                
            elif ptrDir == 4:
                numPath += 1
        
    
    score = score + math.log(numPath)
    print "numPath",numPath
    return score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath





def globalAlignment(seq1, seq2, parameters):
    score, returnalignedSeq1, returnalignedSeq2 = 0, [],[] 
    [wts, wti, wtd, wtm] = parameters    
    #print [wts, wti, wtd, wtm] 
    
    m = len(seq1) + 1
    n = len(seq2) + 1
    
    #sprint " m, n ", m , n
    
    H = np.zeros(m*n, dtype = np.int32).reshape(m,n)
    B = np.zeros(m*n, dtype = np.int32).reshape(m,n)
    # Assign weights 
    
    for i in range(m):
        H[i][0] = wti*i
    for j in range(n):
        H[0][j] = wtd*j
        
    for i in range(1,m):
        for j in range(1,n):
            if seq1[i-1] == seq2[j-1]:
                w = wtm
            else:
                w = wts
                
            H[i][j] = max(H[i-1][j-1]+w,H[i-1][j]+wtd, H[i][j-1]+wti )
    
    score = H[m-1][n-1]
    
    #print H
    
    i, j = m-1,n-1 
    while (i>0 or j >0 ):

        if (i>0 and j >0 ):
            if seq1[i-1] == seq2[j-1]:
                w = wtm
            else:
                w = wts
                
        if i>0 and j >0 and  H[i][j] == H[i-1][j-1]+ w:
            returnalignedSeq1.append(seq1[i-1])
            returnalignedSeq2.append(seq2[j-1])
            
            i = i-1
            j = j-1
                
        elif (i > 0 and H[i][j] == H[i-1][j] + wti):
            returnalignedSeq1.append(seq1[i-1])
            returnalignedSeq2.append(0)
            i =i-1
            
        elif (j > 0 and H[i][j] == H[i][j-1] + wtd):
            returnalignedSeq1.append(0)
            returnalignedSeq2.append(seq2[j-1])
            j = j -1 
        
        
    revalignedSeq1= reverseString(returnalignedSeq1)
    revalignedSeq2= reverseString(returnalignedSeq2)


    #print "After : i,j",i, j 
    return score, revalignedSeq1 , revalignedSeq2


def transformToDesiredForm(alignedSeq1, alignedSeq2):
    newAlignedSeq1, newAlignedSeq2 = [] ,[]
    #print "alignedSeq1", alignedSeq1
    #print  "alignedSeq2", alignedSeq2
    #print len(alignedSeq1), len(alignedSeq2)
    Lseq = len(alignedSeq1)
    #alignedSeq1[1] =0 
    #alignedSeq2[0] = 0
    
    for i in range(Lseq):
        newAlignedSeq1.append(alignedSeq1[i])
        newAlignedSeq2.append(alignedSeq2[i])
    
    
    # Screen the sub and run of same character case 
    i =0 
    while (i <= len(newAlignedSeq1) -2 ):
        if newAlignedSeq1[i] == 0 and newAlignedSeq2[i] == 0 :
            newAlignedSeq1.pop(i)
            newAlignedSeq2.pop(i)
                    
        elif newAlignedSeq1[i] == newAlignedSeq2[i+1] and  newAlignedSeq1[i] == newAlignedSeq2[i+1] and  newAlignedSeq2[i] == 0:
            newAlignedSeq2[i] = newAlignedSeq1[i] 
            newAlignedSeq2[i+1] = 0
            i = i +1 
        elif newAlignedSeq2[i] == newAlignedSeq1[i+1] and  newAlignedSeq2[i] == newAlignedSeq1[i+1] and  newAlignedSeq1[i] == 0:
            newAlignedSeq1[i] = newAlignedSeq2[i] 
            newAlignedSeq1[i+1] = 0
            i = i +1 

        elif  newAlignedSeq1[i] == 0 and newAlignedSeq2[i+1] == 0:
            newAlignedSeq1[i] = newAlignedSeq1[i+1]
            newAlignedSeq1[i+1] =0 
            i = i +1 
        elif  newAlignedSeq2[i] == 0 and newAlignedSeq1[i+1] == 0:
            newAlignedSeq2[i] = newAlignedSeq2[i+1] 
            newAlignedSeq2[i+1] = 0   
            i = i +1          
        else: 
            i = i +1 
   
    #print "newAlignedSeq1", newAlignedSeq1
    #print "newAlignedSeq2", newAlignedSeq2
    
    return newAlignedSeq1, newAlignedSeq2



def updateVote(newAlignedSeq1, newAlignedSeq2,tempVoteTable, istart, jstart, iend, jend ):
    # jstart , jend : for the long reads 
    # (newAlignedSeq1, newAlignedSeq2) : short  , long 
    #print len(newAlignedSeq1), len(newAlignedSeq2) 
    #printSeq(newAlignedSeq1)
    #printSeq(newAlignedSeq2)
    temp = 1
    
    #Lseq = jend - jstart 
    
    Lseq = min(len(newAlignedSeq1), len(newAlignedSeq2))
    cumCount = 0 
    insCumCount = 0
    
    for j in range(Lseq): 
        if newAlignedSeq2[j] == 0 :
            #tempVoteTable.insCount[cumCount + jstart ].append()
            # format : [number of insert, base , count ] 
            subbase = newAlignedSeq1[j] -1 
            currentIns = [insCumCount, subbase, 1 ]
            existInIns = False
            for eachins, dummyindex in zip(tempVoteTable.insCount[cumCount + jstart-1], range(len(tempVoteTable.insCount[cumCount + jstart-1]))):
                #print "eachins" ,eachins
                if eachins[0:2] == currentIns[0:2]:
                    tempVoteTable.insCount[cumCount + jstart-1][dummyindex][2] = eachins[2]+1 
                    existInIns = True
                 #   print "debug",tempVoteTable.insCount[cumCount + jstart-1][dummyindex][2] 
                    
            if not existInIns :
                tempVoteTable.insCount[cumCount -1  +jstart].append(currentIns)
            
            insCumCount = insCumCount + 1 
        elif newAlignedSeq1[j] == 0:
            tempVoteTable.delCount[cumCount+ jstart] += 1 
            
        elif newAlignedSeq2[j] == newAlignedSeq1[j]:
            tempVoteTable.confirmCount[cumCount + jstart] += 1 
            
        elif newAlignedSeq2[j]  != newAlignedSeq1[j]  :
            subbase = newAlignedSeq1[j] -1 
            tempVoteTable.subCount[cumCount+jstart][subbase] += 1 
            
            
        if newAlignedSeq2[j] != 0:
            if (j< len(newAlignedSeq2) -1 and  newAlignedSeq2[j+1] != 0) or j==len(newAlignedSeq2) -1:
                tempVoteTable.confirmNoIns[cumCount + jstart] += 1 
                
            cumCount = cumCount + 1 
            insCumCount = 0 
            
    #printVoteTable(tempVoteTable, jstart)
    #print "updateVote"
def lockCountUpdate(cleanedLongRead, count):
    
    if len(cleanedLongRead) <2 :
        return -1
    else:
        if count >=2 :
            return -1
        
        if count < 0:
            return -1
        
        if cleanedLongRead[-2]!= cleanedLongRead[-1]:
            return count +1
        else:
            return count



    
def consensusForLongReads(tempVoteTable, parameterRobot, roundNum =0):
    L = len(tempVoteTable.longread)
    cleanedLongRead =  [] 
    thresholdForSupport = parameterRobot.thresholdForSupport
    subthreshold = parameterRobot.subthreshold
    
    ratiothresholddel = parameterRobot.thresFordel
    ratiothresholdins = parameterRobot.thresForins
    mindelcount = parameterRobot.delMin
    mininscount = parameterRobot.insMin
    
    countToConfirm  = 0 
    countToDel = 0 
    countToIns = 0
    countToSub = 0


    lock = False
    changeCounter = -1
    
    for indexout in range(L):
        #print indexout, tempVoteTable.SNPlocList
        #if indexout in range(210,220):
        #    print "indexout, changeCounter",indexout, changeCounter
        
        
        insCountSum = 0 
        toAddList = [] 
        # Format toAddList : [placeToAdd, BaseToAdd,count]
        for eachitem,index in zip(tempVoteTable.insCount[indexout], range(len(tempVoteTable.insCount[indexout]))):
            if eachitem[2] >= ratiothresholdins*tempVoteTable.confirmNoIns[indexout] and eachitem[2] >= mininscount:
                toAddList.append(eachitem[0:3])
                
              
        # Filtering insertion
        toAddList = sorted(toAddList)
        runningindex = 0 
        
        potentialSubList = []
        
        while (runningindex < len(toAddList) -1 ):
            if toAddList[runningindex][0] == toAddList[runningindex+1 ][0]:
                toAddList.pop(runningindex) 
                if not toAddList[runningindex][0] in potentialSubList:
                    potentialSubList.append(toAddList[runningindex][0])
            else:
                runningindex = runningindex + 1
        
        #print tempVoteTable.subCount[indexout],np.max(tempVoteTable.subCount[indexout]) 
        # Handle sub, ins, del 
        
        sumOfSub = np.max(tempVoteTable.subCount[indexout]) + tempVoteTable.confirmCount[indexout]
        
        if ( tempVoteTable.delCount[indexout]<=ratiothresholddel*(tempVoteTable.confirmCount[indexout] ) or  tempVoteTable.delCount[indexout] < mindelcount) and (  np.max(tempVoteTable.subCount[indexout]) <= tempVoteTable.confirmCount[indexout])  :
            cleanedLongRead.append(tempVoteTable.longread[indexout])
            
            
            changeCounter = lockCountUpdate(cleanedLongRead, changeCounter)
            
            if np.max(tempVoteTable.subCount[indexout]) >= max(subthreshold, thresholdForSupport*sumOfSub )and tempVoteTable.confirmCount[indexout] >= max(subthreshold, thresholdForSupport*sumOfSub ) :
                tempVoteTable.logSNPloc([len(cleanedLongRead)-1]) 
                #print "Now1"
            
            countToConfirm += 1 
        
            
        elif  np.max(tempVoteTable.subCount[indexout]) > tempVoteTable.confirmCount[indexout] and ratiothresholddel*np.max(tempVoteTable.subCount[indexout])>=tempVoteTable.delCount[indexout]:
            optindex = np.argmax(tempVoteTable.subCount[indexout])
            print tempVoteTable.longread[indexout], tempVoteTable.subCount[indexout][optindex] + 1 
            print indexout
            cleanedLongRead.append( optindex + 1 )
            changeCounter = lockCountUpdate(cleanedLongRead, changeCounter)
            
            if np.max(tempVoteTable.subCount[indexout]) >= max(subthreshold, thresholdForSupport*sumOfSub ) and tempVoteTable.confirmCount[indexout] >= max(subthreshold, thresholdForSupport*sumOfSub ):
                tempVoteTable.logSNPloc([len(cleanedLongRead)-1])        
                #print "Now"
            countToSub += 1 
            
        
        else:
            if roundNum ==0 :       
                countToDel += 1   
            else:
                ### Locking treatment
                if changeCounter < 0 :   
                    countToDel += 1   
                    changeCounter = 0
                    
                else:
                    cleanedLongRead.append(tempVoteTable.longread[indexout])
                    changeCounter = lockCountUpdate(cleanedLongRead, changeCounter)
                    
                    if  np.max(tempVoteTable.subCount[indexout]) > tempVoteTable.confirmCount[indexout] and ratiothresholddel*np.max(tempVoteTable.subCount[indexout])>=tempVoteTable.delCount[indexout]:
                        if np.max(tempVoteTable.subCount[indexout]) >= max(subthreshold, thresholdForSupport*sumOfSub )and tempVoteTable.confirmCount[indexout] >= max(subthreshold, thresholdForSupport*sumOfSub ) :
                            tempVoteTable.logSNPloc([len(cleanedLongRead)-1]) 
                    
                    countToConfirm += 1 
            
                ### End locking    
                
                
        if len(toAddList) > 0 :
            for eachitem in toAddList : 
                
                #for eachpotential in potentialSubList:
                #    tempVoteTable.logSNPloc([len(cleanedLongRead) + eachpotential])  
                
                cleanedLongRead.append(eachitem[1] + 1 ) 
                changeCounter = lockCountUpdate(cleanedLongRead, changeCounter)
                
                countToIns += 1 
    
                print "toAddList", toAddList, indexout, len(cleanedLongRead)
                
    #print "cleanedLongRead:"
    #printSeq(cleanedLongRead)
    #print "tempVoteTable.longread:"
    #printSeq(tempVoteTable.longread)
    print "countToConfirm",  countToConfirm
    print "countToDel ", countToDel 
    print "countToIns", countToIns
    print  "countToSub ", countToSub 
    print "len(cleanedLongRead)",len(cleanedLongRead) 
    
    return cleanedLongRead

def consensusForLongReadsFinalTest(tempVoteTable, parameterRobot, roundNum =1):
    L = len(tempVoteTable.longread)
    cleanedLongRead =  [] 
    thresholdForSupport = parameterRobot.thresholdForSupport
    subthreshold = parameterRobot.subthreshold
    print "subthreshold" ,subthreshold 
    
    ratiothresholddel = parameterRobot.thresFordel
    ratiothresholdins = parameterRobot.thresForins
    mindelcount = parameterRobot.delMin
    mininscount = parameterRobot.insMin
    
    countToConfirm  = 0 
    countToDel = 0 
    countToIns = 0
    countToSub = 0


    lock = False
    changeCounter = -1
    
    for indexout in range(L):
        #print indexout, tempVoteTable.SNPlocList
        #if indexout in range(210,220):
        #    print "indexout, changeCounter",indexout, changeCounter
        
        
        insCountSum = 0 
        toAddList = [] 
        # Format toAddList : [placeToAdd, BaseToAdd,count]
        for eachitem,index in zip(tempVoteTable.insCount[indexout], range(len(tempVoteTable.insCount[indexout]))):
            ### Temp SNP finding Hack
            if roundNum == 2:
                if eachitem[2] > subthreshold:
                    locToAdd = len(cleanedLongRead) - 1 + eachitem[0]
                    if not locToAdd in tempVoteTable.SNPlocList:
                        tempVoteTable.logSNPloc([locToAdd]) 

            
            ### End Hack  
            if eachitem[2] >= ratiothresholdins*tempVoteTable.confirmNoIns[indexout] and eachitem[2] >= mininscount:
                toAddList.append(eachitem[0:3])
                
              
        # Filtering insertion
        toAddList = sorted(toAddList)
        runningindex = 0 
        
        potentialSubList = []
        
        while (runningindex < len(toAddList) -1 ):
            if toAddList[runningindex][0] == toAddList[runningindex+1 ][0]:
                toAddList.pop(runningindex) 
                if not toAddList[runningindex][0] in potentialSubList:
                    potentialSubList.append(toAddList[runningindex][0])
            else:
                runningindex = runningindex + 1
        
        #print tempVoteTable.subCount[indexout],np.max(tempVoteTable.subCount[indexout]) 
        # Handle sub, ins, del 
        
        sumOfSub = np.max(tempVoteTable.subCount[indexout]) + tempVoteTable.confirmCount[indexout]
        
        # Hack 
        if tempVoteTable.confirmCount[indexout] >  subthreshold and np.max(tempVoteTable.subCount[indexout])  > subthreshold:
             if not len(cleanedLongRead) in tempVoteTable.SNPlocList:
                 tempVoteTable.logSNPloc([len(cleanedLongRead) ])

             
        # end hack 
        
        if ( tempVoteTable.delCount[indexout]<= float(ratiothresholddel*(sumOfSub)) or  tempVoteTable.delCount[indexout] < mindelcount) and (  np.max(tempVoteTable.subCount[indexout]) <= tempVoteTable.confirmCount[indexout])  :
            cleanedLongRead.append(tempVoteTable.longread[indexout])
            
            
            changeCounter = lockCountUpdate(cleanedLongRead, changeCounter)
            
            if np.max(tempVoteTable.subCount[indexout]) >= max(subthreshold, thresholdForSupport*sumOfSub )and tempVoteTable.confirmCount[indexout] >= max(subthreshold, thresholdForSupport*sumOfSub ) :
                tempVoteTable.logSNPloc([len(cleanedLongRead)-1]) 
                

            countToConfirm += 1 
           
        elif  np.max(tempVoteTable.subCount[indexout]) > tempVoteTable.confirmCount[indexout] and ratiothresholddel*np.max(tempVoteTable.subCount[indexout])>=tempVoteTable.delCount[indexout]:
            optindex = np.argmax(tempVoteTable.subCount[indexout])
            print tempVoteTable.longread[indexout], tempVoteTable.subCount[indexout][optindex] + 1 
            print indexout
            cleanedLongRead.append( optindex + 1 )
            changeCounter = lockCountUpdate(cleanedLongRead, changeCounter)
            
            if np.max(tempVoteTable.subCount[indexout]) >= max(subthreshold, thresholdForSupport*sumOfSub ) and tempVoteTable.confirmCount[indexout] >= max(subthreshold, thresholdForSupport*sumOfSub ):
                tempVoteTable.logSNPloc([len(cleanedLongRead)-1])        

            countToSub += 1 
            
        
        else:
            if roundNum ==0 :       
                countToDel += 1   
                
            else:
                print "indexout", indexout
                ### Locking treatment
                if changeCounter < 0 :   
                    countToDel += 1   
                    changeCounter = 0
                    tempSNPLoc = len(cleanedLongRead) 
                    
                else:
                    tempVoteTable.logSNPloc([tempSNPLoc])
                    
                    cleanedLongRead.append(tempVoteTable.longread[indexout])
                    changeCounter = lockCountUpdate(cleanedLongRead, changeCounter)
                    
                    if  np.max(tempVoteTable.subCount[indexout]) > tempVoteTable.confirmCount[indexout] and ratiothresholddel*np.max(tempVoteTable.subCount[indexout])>=tempVoteTable.delCount[indexout]:
                        if np.max(tempVoteTable.subCount[indexout]) >= max(subthreshold, thresholdForSupport*sumOfSub )and tempVoteTable.confirmCount[indexout] >= max(subthreshold, thresholdForSupport*sumOfSub ) :
                            tempVoteTable.logSNPloc([len(cleanedLongRead)-1]) 
                            
                    
                    countToConfirm += 1 
            
                ### End locking    
                
                
        if len(toAddList) > 0 :
            for eachitem in toAddList : 
                
                for eachpotential in potentialSubList:
                    tempVoteTable.logSNPloc([len(cleanedLongRead) + eachpotential])  
                
 #               baseIndex = eachitem[1]
 #               if tempVoteTable.subCount[indexout][baseIndex] > subthreshold:
 #                   tempVoteTable.logSNPloc([len(cleanedLongRead)-1])
 #               elif tempVoteTable.subCount[indexout-1][baseIndex] > subthreshold:
 #                   tempVoteTable.logSNPloc([len(cleanedLongRead)-2])
 #               elif tempVoteTable.subCount[indexout+1][baseIndex] > subthreshold:
 #                   tempVoteTable.logSNPloc([len(cleanedLongRead)])
                    
                cleanedLongRead.append(eachitem[1] + 1 )                
                changeCounter = lockCountUpdate(cleanedLongRead, changeCounter)
                
                
                
                
                countToIns += 1 
    
                print "toAddList", toAddList, indexout, len(cleanedLongRead)
                
    #print "cleanedLongRead:"
    #printSeq(cleanedLongRead)
    #print "tempVoteTable.longread:"
    #printSeq(tempVoteTable.longread)
    print "countToConfirm",  countToConfirm
    print "countToDel ", countToDel 
    print "countToIns", countToIns
    print  "countToSub ", countToSub 
    print "tempVoteTable.logSNPloc", len(tempVoteTable.SNPlocList)
    print "len(cleanedLongRead)",len(cleanedLongRead) 
    
    return cleanedLongRead

def printVoteTable(voteTable, jstart):
    
    #[jstart +threshold ]
    for i in range(len(voteTable.longread)):
        print i, voteTable.longread[i], voteTable.delCount[i], voteTable.subCount[i], voteTable.confirmCount[i], voteTable.insCount[i], voteTable.confirmNoIns[i]
        #print "i", i
        #print "Count ", voteTable.longread[i]
        #print "Del Count", voteTable.delCount[i]
       
        #print "Sub Count", voteTable.subCount[i]
        
        #print "Confirm Count", voteTable.confirmCount[i]
            
        #print "Ins Count", voteTable.insCount[i]
        
        #print "confirmNoIns", voteTable.confirmNoIns[i]



def selectivePrint(voteTable, potentialSNPList):
    tmpList = sorted(potentialSNPList)
    extendlength = 3
    
    
    for i in range(len(voteTable.longread)):
        check = False
        for j in range(i-extendlength, i+extendlength + 1):
            if j in potentialSNPList:
                check = True
                
        if check:
            print i, voteTable.longread[i], voteTable.delCount[i], voteTable.subCount[i], voteTable.confirmCount[i], voteTable.insCount[i], voteTable.confirmNoIns[i]
    
def groundTruthViewer(indexlong,readsMapping,longRepeatStat, parameterRobot, motherGenome ):
    copy1, copy2 = [], [] 
    readIndex, readStart, reverse = readsMapping[1][indexlong]
    repeatLoc1, repeatLoc2, lengthOfRep = longRepeatStat[1]
    Llong = parameterRobot.Llong
    
    if readStart  > repeatLoc2:
        delta = -repeatLoc2 + repeatLoc1
        SNP1, SNP2 = repeatLoc2 +lengthOfRep /4 - readStart, repeatLoc2 +3*lengthOfRep /4 - readStart
    else:
        delta = -repeatLoc1 +repeatLoc2
        SNP1, SNP2 = repeatLoc1 +lengthOfRep /4 - readStart, repeatLoc1 +3*lengthOfRep /4 - readStart
    
    copy1 = motherGenome[readStart: readStart+ Llong]
    copy2 = motherGenome[readStart+ delta: readStart+delta+ Llong]
    
    if reverse == 1: 
        copy1 = reverseStrand(copy1)
        copy2 = reverseStrand(copy2)
    
    
    return copy1, copy2 , SNP1, SNP2
#print globalAlignment(np.ones(200), np.ones(200), [-10,-1,-1,1])

### MAP Cleaning



class graphNode(object):
    def __init__(self,id):
        self.id = id
        self.content = ""
        self.childrenList = []
        self.parentList = []
        
    def setChildren(self, childrenList):
        for eachitem in childrenList:
            self.childrenList.append(eachitem)
    
    def setParentList(self,parentList):
        for eachitem in parentList:
            self.parentList.append(eachitem)
            
class choiceGraph(object):
    def __init__(self):
        self.graphNodes = []
        self.tmpList = []
        self.pathList = []
        self.source = []
        self.sink = []
        self.correspondingSNP = -100
        self.pairList = []
    def reset(self):
        self.tmpList = []
            
    def initNodes(self):
        print "init nodes"
    def printAllPaths(self):
        print "\nAll paths List for SNP",  self.correspondingSNP
        for eachitem in self.pathList:
            print eachitem
            
        print "End All paths List : CAUTION: Beyond this point , you may get an A+\n "
    def filterPaths(self):
        
        i =0 
        while (i < len(self.pathList) -1 ):
            check = False
            for j in range(i +1, len(self.pathList)):
                if self.pathList[i] == self.pathList[j]:
                    self.pathList.pop(i)
                    check = True 
                    break
            
            if not check :
                i = i +1 
    def groupPath(self, typeOfMode =0 ):
        residualList = range(len(self.pathList))
        for i in range(len(self.pathList)-1):
            for j in range(i+1, len(self.pathList)):
                #if len(self.pathList[i]) == len(self.pathList[j]) and hammingDistance(self.pathList[i], self.pathList[j]) == 1:
                if len(self.pathList[i]) == len(self.pathList[j]) :
                 
                    self.pairList.append([i,j]) 
                    if i in residualList:
                        residualList.remove(i)
                    if j in residualList:
                        residualList.remove(j)
        if typeOfMode == 0:            
            for eachitem in residualList:
                self.pairList.append([eachitem])
        elif typeOfMode == 1:
            for eachitem in range(len(self.pathList)):
                self.pairList.append([eachitem])
        
        
        self.voteForPairList= [0 for i in range(len(self.pairList))]
               
            
    def findAllPaths(self):
        self.tmpList.append(self.source.content)
        self.DFS(self.source)
        self.printAllPaths()
        
    def DFS(self,x):
        if x.content == "E":
            mylist = []
            for eachitem in self.tmpList:
                mylist.append(eachitem)
            #self.pathList.append(mylist)
            self.pathList.append(mylist[1:-1])
        else:
            for y in x.childrenList:
                self.tmpList.append(y.content)
                #print self.tmpList
                self.DFS(y)
                self.tmpList.pop()
            
def filterSNPList(SNPlist, lookRange):          
    tempSNPList = sorted(SNPlist)
    
    runningindex  = 0 
    olditem = -1000 
    while (runningindex < len(tempSNPList))  :
        if tempSNPList[runningindex] < olditem + lookRange :
            tempSNPList.pop(runningindex)
        else:
            
            olditem = tempSNPList[runningindex]
            runningindex += 1 
    
    return tempSNPList

def identifyChoices(SNPlist, voteTable, parameterRobot):

    trustRange = parameterRobot.lookRange / 3
    lookRange = parameterRobot.lookRange
    
    SNPlist = filterSNPList(SNPlist,trustRange)
    templateGraphList = []
    
    for eachsnp in SNPlist:
        lbdd, upbdd = max(0,eachsnp - lookRange), min(len(voteTable.longread) , eachsnp+ lookRange)
        templateGraph = formAlphabetGraph(lbdd, upbdd, voteTable, parameterRobot)
        templateGraph.correspondingSNP = eachsnp
        templateGraph.findAllPaths()
        ### Hack for computation concern but will increase false negative rate: 
        if len(templateGraph.pathList) <30:
            templateGraphList.append(templateGraph)
        
    return templateGraphList

def setEdge(parent, child):
    parent.childrenList.append(child)
    child.parentList.append(parent)

def formAlphabetGraph(lbdd, upbdd, voteTable, parameterRobot):
    # From Graph
    subTresholdRat = parameterRobot.thresholdForSupport
    subthreshold = parameterRobot.subthreshold
    thresForinsRat , thresFordelRat=parameterRobot.thresForins , parameterRobot.thresFordel
    #print "lbdd, upbdd", lbdd, upbdd
    
    # Backbone
    alGraph = choiceGraph()
    startNode = graphNode(0)
    startNode.content = "S"
    alGraph.source = startNode
    alGraph.graphNodes.append(startNode)

    for i in range(lbdd,upbdd):
        alphabet = voteTable.longread[i]
        nodeIndex = len(alGraph.graphNodes)
        newNode = graphNode(nodeIndex)
        newNode.content = alphabet
        alGraph.graphNodes.append(newNode)
    
    endNode = graphNode(-1)
    endNode.content = "E"
    
    alGraph.graphNodes.append(endNode)
    alGraph.sink = endNode
    
    # Add edges and sides branchees
    for i in range(len(alGraph.graphNodes)-1):
        setEdge(alGraph.graphNodes[i], alGraph.graphNodes[i+1])
        
        
    
    for i in range(1, upbdd - lbdd + 1):
        # Single delete treatment
        if voteTable.delCount[i+lbdd -1] > thresFordelRat*voteTable.confirmCount[i+lbdd -1]:
            setEdge(alGraph.graphNodes[i-1], alGraph.graphNodes[i+1])
            if i+2 <= upbdd - lbdd + 1  and voteTable.delCount[i+1+lbdd -1] > thresFordelRat*voteTable.confirmCount[i+1+lbdd -1]:
                setEdge(alGraph.graphNodes[i-1], alGraph.graphNodes[i+2])
        # Double delete treatment
        
        
        lvlzeroInsert = []
        lvloneInsert = []
        # zero level insert
        for item in voteTable.insCount[i+lbdd -1]:

            if (item[2] >  thresForinsRat*voteTable.confirmCount[i+lbdd -1] or item[2] > subthreshold )and item[0] == 0:
                newNode = graphNode(len(alGraph.graphNodes))
                newNode.content = item[1] + 1 
                alGraph.graphNodes.append(newNode)
                
                setEdge(alGraph.graphNodes[i], newNode)
                setEdge(newNode,alGraph.graphNodes[i+1])
                
                lvlzeroInsert.append(newNode)
        # levvel 1 insert 
        for item in voteTable.insCount[i+lbdd -1]:

            if (item[2] >  thresForinsRat*voteTable.confirmCount[i+lbdd -1] or item[2] > subthreshold )and item[0] == 1:
                newNode = graphNode(len(alGraph.graphNodes))
                newNode.content = item[1] + 1 
                alGraph.graphNodes.append(newNode)
                
                for eachnode in lvlzeroInsert:
                    setEdge(eachnode, newNode)
                
                setEdge(newNode,alGraph.graphNodes[i+1])
                lvloneInsert.append(newNode)
        # levvel 2 insert 
        for item in voteTable.insCount[i+lbdd -1]:

            if (item[2] >  thresForinsRat*voteTable.confirmCount[i+lbdd -1] or item[2] > subthreshold )and item[0] == 2:
                newNode = graphNode(len(alGraph.graphNodes))
                newNode.content = item[1] + 1 
                alGraph.graphNodes.append(newNode)
                
                for eachnode in lvloneInsert:
                    setEdge(eachnode, newNode)
                
                setEdge(newNode,alGraph.graphNodes[i+1])       
                
        # Sub case
        for j in range(4):
            if voteTable.subCount[i+lbdd -1][j] > subTresholdRat*voteTable.confirmCount[i+lbdd -1] or voteTable.subCount[i+lbdd -1][j]  > subthreshold:
                newNode = graphNode(len(alGraph.graphNodes))
                newNode.content= j +1 
                alGraph.graphNodes.append(newNode)
                
                for eachparent in alGraph.graphNodes[i].parentList:
                    setEdge(eachparent, newNode)
                    
                for eachchild in alGraph.graphNodes[i].childrenList:
                    setEdge(newNode,eachchild)
                
                
    return alGraph

def hammingDistance(seq1, seq2):
    count = 0 
    for i in range(min(len(seq2), len(seq1))):
        if seq1[i] != seq2[i]:
            count += 1  
    return count 



def computeSumScore(eachtemplate,eachGraph, shortToLongMap, shortReadsList, parameterRobot, longread):
    sumscore = 0
    lookRange = parameterRobot.lookRange *2 
    editParameters = [parameterRobot.editsub, parameterRobot.editins, parameterRobot.editdel, parameterRobot.editmatch]
    
    targetSNPloc = eachGraph.correspondingSNP
    readSegmentList = []
    
    for eachmap in shortToLongMap:
        indexlong ,   indexshort  ,  jstart  ,  jend  ,  istart   , iend= eachmap
        if jstart <= max(0, targetSNPloc - lookRange) and jend >= min( targetSNPloc +lookRange, parameterRobot.Llong - 1 ) :
            tempRead = []
            sStart , sEnd = max(istart- jstart +targetSNPloc - lookRange, 0)   , min(istart- jstart +targetSNPloc+lookRange, len(shortReadsList[indexshort]) )
            #sStart , sEnd  = sStart  -4 , sEnd -4
            tempRead = shortReadsList[indexshort][sStart:sEnd]
            
            if eachGraph.correspondingSNP == 192: 
                tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj= SWAlignment(tempRead, eachtemplate, parameterRobot) 
                #tempScore, returnalignedSeq1, returnalignedSeq2  =globalAlignment( eachtemplate,tempRead, editParameters)
                
                printSeq(returnalignedSeq1)
                printSeq(returnalignedSeq2)
                print "tempScore", tempScore
                
                gdscore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj =  SWAlignment(longread , shortReadsList[indexshort], parameterRobot)
                
                printSeq(returnalignedSeq1)
                printSeq(returnalignedSeq2)
                print "gdscore", gdscore
                
                print "tempScore, indexshort, sStart , sEnd ",tempScore,indexshort,  sStart , sEnd  
            readSegmentList.append(tempRead)
    
    
    for eachsegment in readSegmentList:
        
        #print len(eachsegment), len(eachtemplate)
        tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj= SWAlignment(eachsegment, eachtemplate, parameterRobot)
        #tempScore, returnalignedSeq1, returnalignedSeq2  =globalAlignment(eachsegment, eachtemplate, editParameters)
        #printSeq(revalignedSeq1)
        #printSeq(revalignedSeq2)
                
        sumscore += tempScore
        
    #print " eachGraph.correspondingSNP :", eachGraph.correspondingSNP
    if len(readSegmentList) == 0:
        return 0 
    
    averageScore = sumscore/ len(readSegmentList)

    return averageScore

def findBestSingletonScore(voteAns, pairList):
    scoreSingletonMax = -1000
    tmpIndex = -1
    
    for i in range(len(pairList)):
        if len(pairList[i]) == 1:
            item = pairList[i][0]
            if voteAns[item] > scoreSingletonMax:
                tmpIndex = i
                scoreSingletonMax = voteAns[item]
    return tmpIndex
    
def establishVoteOnTemplate(eachGraph,shortToLongMap,shortReadsList,parameterRobot):
    voteAns  = []
    
    voteAns = [0 for i in range(len(eachGraph.pathList))]
    
    lookRange = parameterRobot.lookRange 
    editParameters = [parameterRobot.editsub, parameterRobot.editins, parameterRobot.editdel, parameterRobot.editmatch]
    
    targetSNPloc = eachGraph.correspondingSNP
    readSegmentList = []
    
    for eachmap in shortToLongMap:
        indexlong ,   indexshort  ,  jstart  ,  jend  ,  istart   , iend= eachmap
        #print "jstart, jend", jstart, jend
        if jstart <= max(0, targetSNPloc - lookRange) and jend >= min( targetSNPloc +lookRange, parameterRobot.Llong - 1 ) :
            tempRead = []
            sStart , sEnd = max(istart- jstart +targetSNPloc - lookRange*2, 0)   , min(istart- jstart +targetSNPloc+lookRange*2, len(shortReadsList[indexshort]) )
            #sStart , sEnd  = sStart  -4 , sEnd -4
            tempRead = shortReadsList[indexshort][sStart:sEnd]
            
            readSegmentList.append(tempRead)
    
    print "len(readSegmentList)", len(readSegmentList)
    for eachsegment in readSegmentList:
        scoreList = []
        for eachtemplate,dummyi in zip(eachGraph.pathList,range(len(eachGraph.pathList))):
            tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj= SWAlignment(eachsegment, eachtemplate, parameterRobot)
            scoreList.append(tempScore)
        maxScore = max(scoreList)
        for j in range(len(scoreList)):
            if maxScore == scoreList[j]:
                voteAns[j] += 1 
        
        
            
    return voteAns

def carefulConsensus3(templateGraphList,shortToLongMap,shortReadsList,parameterRobot, cleanedLongRead,longReadIndex ):
    ultimateSNPList, segList = [], [] 
    
    for i in range(len(templateGraphList)):
        scoreMatrix = formScoreMatrix(i, templateGraphList,shortToLongMap,shortReadsList,parameterRobot, cleanedLongRead,longReadIndex )
        hypothesisChoices = findHypothesisChoice(templateGraphList[i])        
        finalAns = decideVote(scoreMatrix, hypothesisChoices, parameterRobot, templateGraphList[i])
        print "finalAns"
        
        if len(finalAns) == 2:
      #      print "templateGraphList[i].pathList[finalAns[0:2]]"
      #      print templateGraphList[i].pathList[finalAns[0]]
      #      print templateGraphList[i].pathList[finalAns[1]]
            ultimateSNPList.append(templateGraphList[i].correspondingSNP)
            segList.append([templateGraphList[i].pathList[finalAns[0]],templateGraphList[i].pathList[finalAns[1]]])
        elif len(finalAns) == 1:
            print "templateGraphList[i].pathList[finalAns[0]]"
            print templateGraphList[i].pathList[finalAns[0]]
    
    return ultimateSNPList, segList 
            
def formScoreMatrix(i, templateGraphList,shortToLongMap,shortReadsList,parameterRobot, cleanedLongRead,longReadIndex ):
    scoreMatrix = []
    T = len(templateGraphList[i].pathList)
    eachGraph = templateGraphList[i]
    perr = parameterRobot.p
    
    
    lookRange = parameterRobot.lookRange 
    editParameters = [parameterRobot.editsub, parameterRobot.editins, parameterRobot.editdel, parameterRobot.editmatch]
    
    targetSNPloc = eachGraph.correspondingSNP
    readSegmentList = []
    
    for eachmap in shortToLongMap:
        indexlong ,   indexshort  ,  jstart  ,  jend  ,  istart   , iend= eachmap
        if indexlong == longReadIndex:
            if jstart <= max(0, targetSNPloc - lookRange) and jend >= min( targetSNPloc +lookRange, parameterRobot.Llong - 1 ) :
                tempRead = []
                sStart , sEnd = max(istart- jstart +targetSNPloc - lookRange*2, 0)   , min(istart- jstart +targetSNPloc+lookRange*2, len(shortReadsList[indexshort]) )
                tempRead = shortReadsList[indexshort][sStart:sEnd]
                
                readSegmentList.append(tempRead)
    
    D = len(readSegmentList)
    
    scoreMatrix = [[ [] for j in range(D)] for i in range(T)]
    
    for eachtemplate, t in zip(eachGraph.pathList, range(len(eachGraph.pathList))):
        for eachsegment,d  in zip(readSegmentList, range(len(readSegmentList))):
            tempScore, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj, numPath= SWAlignmentFixRefMPQS(eachtemplate,eachsegment,  parameterRobot)
            #countConfirm, countDel, countIns, countSub = countEdits(returnalignedSeq2, returnalignedSeq1) # care for the 1,2 and its fcn
            #scoreMatrix[t][d] = [countConfirm, countDel, countIns, countSub, len(eachtemplate),numPath]          
            scoreMatrix[t][d] = tempScore
        
    #scoreMatrix[t][d] = [countConfirm, countDelIns, countSub]          
    return scoreMatrix

def countEdits(returnalignedSeq1, returnalignedSeq2):
    countConfirm, countDel, countIns= 0, 0, 0 
    countSub = 0 
    for i in range(len(returnalignedSeq1)):
        if returnalignedSeq1[i] == 0  :
            countDel += 1 
        elif returnalignedSeq2[i] == 0:
            countIns += 1 
        elif returnalignedSeq1[i] == returnalignedSeq2[i]:
            countConfirm += 1 
        else :
            countSub += 1 
    return countConfirm, countDel, countIns, countSub
            
def findHypothesisChoice(eachGraph):
    hypothesisChoices = []
    eachGraph.groupPath(1)
    hypothesisChoices = eachGraph.pairList
    
    return hypothesisChoices        
 
def decideVote(scoreMatrix, hypothesisChoices, parameterRobot, eachGraph ):
    finalAns = []  # Being the Associated SNP/bases 
    perr = parameterRobot.p
    thresholdForLength = 0.08
    SNPrate = 0.01
    
    scoreForTemplate = [ 0 for i in range(len(hypothesisChoices))]

    for j  in range(len(hypothesisChoices)):
        # Approximate the solution 
        if len(hypothesisChoices[j]) == 1:
            indext1 = hypothesisChoices[j][0]
            tmpScore = math.log(1-SNPrate)
            
            for d in range(len(scoreMatrix[indext1])):
                #countConfirm, countDel, countIns, countSub, lengthSeed, numPath =  scoreMatrix[indext1][d]
                #confirmNoDelete,confirmNoInsert = countConfirm ,  lengthSeed - countIns
                 
                #tmpScore += math.log(1-perr)*confirmNoDelete + math.log(1-perr)*confirmNoInsert + math.log(perr)*countDel + math.log(perr*perr/3)*countSub + math.log(perr/3)*countIns + math.log(numPath)
                tmpScore += scoreMatrix[indext1][d]
            scoreForTemplate[j]  = tmpScore
        elif len(hypothesisChoices[j]) == 2:
            indext1, indext2 = hypothesisChoices[j][0], hypothesisChoices[j][1]
            indexForUse = 0
            tmpScore = math.log(SNPrate)
                
            for d in range(len(scoreMatrix[0])):
                #countConfirm, countDel, countIns, countSub, lengthSeed, numPath =  scoreMatrix[indext1][d]
                #confirmNoDelete,confirmNoInsert = countConfirm ,  lengthSeed - countIns
                #s1 = math.log(1-perr)*confirmNoDelete + math.log(1-perr)*confirmNoInsert + math.log(perr)*countDel + math.log(perr*perr/3)*countSub + math.log(perr/3)*countIns+ math.log(numPath)
                s1  = scoreMatrix[indext1][d]
                
                #countConfirm, countDel, countIns, countSub, lengthSeed, numPath =  scoreMatrix[indext2][d]
                #confirmNoDelete,confirmNoInsert = countConfirm ,  lengthSeed - countIns
                #s2 = math.log(1-perr)*confirmNoDelete + math.log(1-perr)*confirmNoInsert + math.log(perr)*countDel + math.log(perr*perr/3)*countSub + math.log(perr/3)*countIns+ math.log(numPath)
                s2 = scoreMatrix[indext2][d]
                
                if s1 >s2:
                    tmpScore = tmpScore +  s1 + math.log(0.5)
                elif s1 < s2:
                    tmpScore =tmpScore + s2 + math.log(0.5)
                else:
                    tmpScore = tmpScore + s1 
                    
            scoreForTemplate[j]  = tmpScore
        
    for i in range(len(hypothesisChoices)):
        print "i", i
        print "scoreForTemplate", scoreForTemplate[i]
        print "hypothesisChoices",hypothesisChoices[i]
    
        
    maxScoreIndex = np.argmax(scoreForTemplate)
    maxScore = max(scoreForTemplate)
    finalAns = hypothesisChoices[maxScoreIndex]

    # Filtering steps 1 ( length filter for ties) :
    
    if len(finalAns) > 1:
        singletonTop = -10000
        singletonIndex = -1
        for i in range(len(hypothesisChoices)):
            if len(hypothesisChoices[i]) == 1:
                if scoreForTemplate[i] > singletonTop:
                    singletonTop = scoreForTemplate[i]
                    singletonIndex = hypothesisChoices[i][0]
                    
        #print singletonTop, singletonIndex
       # if abs(maxScore- singletonTop) < thresholdForLength*abs(maxScore):
       #     if len(eachGraph.pathList[singletonIndex]) > len(eachGraph.pathList[finalAns[0]]):
       #         finalAns = [singletonIndex]
                
                    
        # Filtering steps 2 (compute filter ) 
        finalAns = filterLongRuns(finalAns, eachGraph, singletonIndex)

    return finalAns 

def filterLongRuns(finalAns, eachGraph, singletonIndex):
    lookRun = 5
    if len(finalAns) > 1 : 
        path1Index, path2Index = finalAns[0], finalAns[1]
        path1, path2 = eachGraph.pathList[path1Index], eachGraph.pathList[path2Index]
        j = -1 
        for i in range(len(path1)):
            if path1[i] != path2[i]:
                j = i
                
        check1 = True

        if j-lookRun <0 :
            check1 = False
        else:
            tempx = path1[j-1]
            for i in range(j-lookRun, j):
                if tempx != path1[i]:
                    check1 = False
        
        check2 = True
        if j + lookRun >= len(path1):
            check2 = False
        else:
            tempy =  path1[j+1]
            for i in range(j+1, j+lookRun+1):
                if tempy != path1[i]:
                    check2 = False
        
            
        if check1 or check2:

            searchItem = path1[0:j] + path1[j+1:len(path1)]
  
            found = False
            
            for item in eachGraph.pathList:
                if item == searchItem:
                    found = True
            if not found :
                finalAns = [singletonIndex]
                
    return finalAns


###


def dataTypeConversion(noisyReads):
    shortReads = noisyReads[0]
    longReads = noisyReads[1]

    returnShortReads = []
    returnLongReads = []
    
    for eachitem in shortReads:
        endNoisy1= 0 
        while (eachitem[endNoisy1] != 0):
            endNoisy1 += 1
        returnShortReads.append(eachitem[0:endNoisy1])
        
    
    for eachitem in longReads:
        endNoisy1= 0 
        while (eachitem[endNoisy1] != 0):
            endNoisy1 += 1
        returnLongReads.append(eachitem[0:endNoisy1])
        
    
    
    return returnShortReads, returnLongReads 
def cleaning(noisyReads,shortToLongMap, toProcessListInput,parameterRobot, debug = "init"):
### Main CLeaning Function
    print "Cleaning"
    print "debug", debug
    
    # Parameters and setting 
    Nlong, Nshort, Llong, Lshort = parameterRobot.Nlong, parameterRobot.Nshort, parameterRobot.Llong, parameterRobot.Lshort
    liid = parameterRobot.liid
    thresForRandom = parameterRobot.thresForRandom
    overhang = parameterRobot.liid 
    editParameters = [parameterRobot.editsub, parameterRobot.editins, parameterRobot.editdel, parameterRobot.editmatch]
    
    in1IndexList, in2IndexList, out1IndexList, out2IndexList, commonIndexList = toProcessListInput[0],toProcessListInput[1],toProcessListInput[2],toProcessListInput[3],toProcessListInput[4]
    toProcessList =  in1IndexList+ in2IndexList+ out1IndexList+ out2IndexList+ commonIndexList

    numberOfRound = 2
    longReadToUse = commonIndexList[0]
    print "longReadToUse", longReadToUse
    #print "toProcessList", toProcessList
    
    ### Proprocess Dtype : 
    
    ''' 
    Need some hacks to change into desired form 
    '''
    
    #noisyReads[0], noisyReads[1]  = dataTypeConversion(noisyReads)
    
    
    longReads = noisyReads[1]
    shortReadsList = noisyReads[0]
    
    #shortReadsList = noisyReads[0]  + transformToReverseStrand(noisyReads[0] )
    '''
    end hacking the noisy reads
    '''
    
    # Format :  "indexlong", "indexshort",  "jstart", "jend", "istart", "iend"

    tempMapping = [[] for i in range(Nlong)]
    for i in range(len(shortToLongMap)):
        indexToAdd = shortToLongMap[i][0]
        indexshort = shortToLongMap[i][1]
        tempMapping[indexToAdd].append(indexshort)
    
    ### Debugging
    for item in tempMapping:
        item = sorted(item)
    #    if len(item)!= 0:
    #        print len(item)
    #        print item
            
    # End Debug
    
    #print "tempMapping", tempMapping # need filter same item
    
    # Obtain Ideal seed reads from short reads and "known" mapping 
    
    longSeedReads = []
    voteTableList = []
    cleanedLongRead = []
    
    maxTrialNum = min(3, len(commonIndexList))
    
    defaultList = commonIndexList[0:maxTrialNum]
    
    if debug == "init":

        for indexlong in toProcessList:
            
            #eachlongread = longReads[indexlong][0:-liid]
            eachlongread = longReads[indexlong]
            
            tempVoteTable = common.voteTable(indexlong, eachlongread) 
            tempVoteTable.initVoteTable()            
            
            if indexlong in defaultList:
                for roundNum in range(numberOfRound):
                    if roundNum == 0:
                        #eachlongread = longReads[indexlong][0:-liid]
                        eachlongread = longReads[indexlong][0:-int(liid/2)]
                    else:
                        eachlongread = cleanedLongRead
               
                    tempVoteTable = common.voteTable(indexlong, eachlongread) 
                    tempVoteTable.initVoteTable()    
                     
                    for indexshort,ckindex in zip(tempMapping[indexlong], range(len(tempMapping[indexlong]))):
                        
                        eachshortreaad = shortReadsList[indexshort]
                        score, alignedSeq1, alignedSeq2, istart, jstart, iend, jend = SWAlignment(eachshortreaad, eachlongread, parameterRobot )
                        check = meetRequirement(score, alignedSeq1, alignedSeq2,eachshortreaad, eachlongread,istart, jstart, iend, jend, thresForRandom, liid ,overhang )
                        
                        if check: 
                            newAlignedSeq1, newAlignedSeq2 = transformToDesiredForm(alignedSeq1, alignedSeq2)
                            updateVote(newAlignedSeq1, newAlignedSeq2,tempVoteTable,istart, jstart, iend, jend )
                            
                    cleanedLongRead = consensusForLongReads(tempVoteTable, parameterRobot,roundNum)
                
            else:
                
                cleanedLongRead= eachlongread

            longSeedReads.append(cleanedLongRead)
            voteTableList.append(tempVoteTable)

        loggingIndel.logVoteTable(parameterRobot.defaultFolder,voteTableList)    
        
    elif debug == "map" or debug == "vote":
        parameterRobot.thresFordel = 0.3
        parameterRobot.thresForins = 0.2
        
        voteTableList = loggingIndel.loadVoteTable(parameterRobot.defaultFolder)
        print "defaultList", defaultList
        maxScoreForChoice, longReadToUse = 10000 , commonIndexList[0]
        
        trial = 0 
        
        for  dummy in range(len(voteTableList)):
            
            tempVoteTable = voteTableList[dummy]
            longReadIndex = tempVoteTable.index
            
            if np.sum(tempVoteTable.confirmCount[:]) == 0:
                #cleanedLongRead = consensusForLongReadsFinalTest(tempVoteTable,parameterRobot,2)
                cleanedLongRead = tempVoteTable.longread
                longSeedReads.append(cleanedLongRead)
                print "Hi"
            elif  np.sum(tempVoteTable.confirmCount[:])> 0 :
                
                print "Doing the final vote"
                indexlong = tempVoteTable.index
                
                tempVoteTable.SNPlocList = []
                
                cleanedLongRead = consensusForLongReadsFinalTest(tempVoteTable,parameterRobot,2)
                longSeedReads.append(cleanedLongRead)
                tempVoteTable.filterSNP()
                
                eachlongread = cleanedLongRead
                
                printVoteTable(tempVoteTable, 0)
                print "indexlong: Potential SNP", indexlong, tempVoteTable.SNPlocList

                
                templateGraphList= identifyChoices(tempVoteTable.SNPlocList, tempVoteTable, parameterRobot)

                selectivePrint(tempVoteTable, tempVoteTable.SNPlocList)          
                    
                  
                ultimateSNPList, segList = carefulConsensus3(templateGraphList,shortToLongMap,shortReadsList,parameterRobot, cleanedLongRead,longReadIndex )

                tempVoteTable.segmentList = segList
                tempVoteTable.SNPlocList = ultimateSNPList
                tempVoteTable.filterSegmentList(parameterRobot)
                ultimateSNPList = tempVoteTable.SNPlocList
                
                if tempVoteTable.dist() < maxScoreForChoice  and (len(ultimateSNPList) > 1 or 10000 == maxScoreForChoice ):
                    maxScoreForChoice, longReadToUse = tempVoteTable.dist() , tempVoteTable.index
                    loggingIndel.logSegList(parameterRobot.defaultFolder ,tempVoteTable.leftSegList, tempVoteTable.rightSegList)
                    
                    print "BADDDD", maxScoreForChoice, longReadToUse,tempVoteTable.leftSegList, tempVoteTable.rightSegList, len(ultimateSNPList)

                print  "ultimateSNPList" ,ultimateSNPList       
                ftmp = open(parameterRobot.defaultFolder + "finalRound_40_49.txt", 'a')
                ftmp.write(str(ultimateSNPList))
                ftmp.write("\n")
                ftmp.close()         
                trial += 1 
            else:
                print tempVoteTable.confirmCount[11]
        
        #loggingIndel.logVoteTable(parameterRobot.defaultFolder,voteTableList)  

    elif debug == "bridge":
        leftSegList, rightSegList = loggingIndel.loadSegList(parameterRobot.defaultFolder)
        voteTableList = loggingIndel.loadVoteTable(parameterRobot.defaultFolder)
        
        print "len(leftSegList), len(rightSegList)" , len(leftSegList), len(rightSegList)
        for  dummy in range(len(voteTableList)):
            tempVoteTable = voteTableList[dummy]   
            cleanedLongRead = tempVoteTable.longread 
            longSeedReads.append(cleanedLongRead)
            if tempVoteTable.index in defaultList:
                #print "Loading ddd"
                tempVoteTable.leftSegList = leftSegList 
                tempVoteTable.rightSegList = rightSegList
            
    
    # Format OutputReads with SNP location (Process longSeedReads, and voteTalbeList)
    in1List, in2List, out1List, out2List, commonList = [],[],[],[],[]
    

     
    for i in range(len(longSeedReads)):
        indexlong = voteTableList[i].index
        eachlongread = longSeedReads[i]
        
        tempVoteTable = common.voteTable(indexlong, eachlongread) 
        tempVoteTable.initVoteTable()
        tempVoteTable.SNPlocList = voteTableList[i].SNPlocList
        tempVoteTable.leftSegList, tempVoteTable.rightSegList = voteTableList[i].leftSegList, voteTableList[i].rightSegList
        
        if i < len(in1IndexList): 
            in1List.append(tempVoteTable)
            
        elif i < len(in1IndexList)+ len(in2IndexList):
            in2List.append(tempVoteTable)
            
        elif i < len(in1IndexList)+ len(in2IndexList) + len(out1IndexList):
            out1List.append(tempVoteTable)
        
        elif i < len(in1IndexList)+ len(in2IndexList) + len(out1IndexList) +  len(out2IndexList):
            out2List.append(tempVoteTable)
            
        elif i <  len(in1IndexList)+ len(in2IndexList) + len(out1IndexList) + len(out2IndexList)+ len(commonIndexList): 
            commonList.append(tempVoteTable)
            
        else:
            print "error"
                
    
    return in1List, in2List, out1List, out2List, commonList , longReadToUse


    
import sys 
from max import *
from matrix import *

def loadSeq(fileName):
    header=None
    dna=''
    file=open(fileName,'r')
    
    lines=file.readlines()
    
    for line in lines:
        if line.startswith('>'):
            header=line 
        else:
            dna=dna+line.strip()
        
    file.close()

    return dna

def globalalignment(seq1, seq2, outputfile, scoringmatrix):

    #Dynamic Programming Solution

    lenv=len(seq1)
    lenw=len(seq2)
    s = [[0 for x in range(lenw+1)] for x in range(lenv+1)]  # +1 to account for first row and column which has gaps
    
    #Scoring Matrix
    gp=int(scoringmatrix[-1][0])
    print("gp: ", gp)
    

    for j in range(lenw+1): #seq w
        s[0][j]=0-j #first row

    for i in range(lenv+1): #seq v
        s[i][0]= 0-i #first column
    
    for i in range(1,lenv+1):
        for j in range(1, lenw+1):
            max=-999
            cur=-999

            #Case 1 (Diagonal)
            idx1=scoringmatrix[0].index(seq1[i-1])
            idx2=scoringmatrix[0].index(seq2[j-1])
            if seq1[i-1]==seq2[j-1]:
                cur=s[i-1][j-1] + int(scoringmatrix[idx1+1][idx2+1]) #match score
            else:
                cur=s[i-1][j-1] + int(scoringmatrix[idx1+1][idx2+1]) #mismatch score
            
            if cur>max:
                max=cur
            
            #Case 2 (Vertical)
            cur=s[i-1][j] + gp #gap penalty

            if cur>max:
                max=cur

            #Case 3 (Horizontal)
            cur= s[i][j-1] + gp #gap penalty

            if cur>max:
                max=cur

            #Putting the value in the table
            s[i][j]=max 
    print("--------------------------------------------------------------")
    #printing table
    for i in s:
        print('\t'.join(map(str, i)))
    print("--------------------------------------------------------------")

    #Trackback
    identity=0
    result= s[lenv][lenw]
    print("Result: ", result)
    print("--------------------------------------------------------------")
    sequence1=""
    sequence2=""

    i=lenv-1
    j=lenw-1
    
    while i>=0 and j>=0:
        print("i: ", i, "j: ", j)
        score=s[i+1][j+1] 
        print("Score: ", score) 

        #Diagonal Score
        idx1=scoringmatrix[0].index(seq1[i])
        idx2=scoringmatrix[0].index(seq2[j])
        if seq1[i]==seq2[j]:
            score1=s[i][j]+ int(scoringmatrix[idx1+1][idx2+1]) #match score
            identity=identity+1 
        else:
            score1=s[i][j]+ int(scoringmatrix[idx1+1][idx2+1]) #mismatch score
        print("Score1: ", score1)

        #Vertical Score:
        score2= s[i][j+1] + gp #gap penalty 
        print("Score2: ", score2)

        #Horizontal Score:
        score3= s[i+1][j] + gp #gap penalty
        print("Score3: ",score3)
        
        if score==score1: #Diagonal
            sequence1=sequence1 + seq1[i]
            sequence2=sequence2 + seq2[j]
            i=i-1
            j=j-1
            
        
        elif score==score2: #Vertical
            sequence1=sequence1+ seq1[i]
            sequence2=sequence2+ "_"
            i=i-1

        elif score==score3: #Horizontal
            sequence1=sequence1+"_"
            sequence2=sequence2+seq2[j]
            j=j-1
    
        print("=============================")

    text1= "seq1: ", "1 ", sequence1[::-1], " ", str(maximum(lenv, lenw)), "\n"
    text2= "seq2: ", "2 ", sequence2[::-1], " ", str(maximum(lenv, lenw)), "\n"
    text3= "Score: ", str(result), "\n"
    text4="Identities: ", str(identity),"/",str(maximum(lenv, lenw)), "  ", "(", str("{:.0%}".format(identity/maximum(lenv, lenw))), ")"

    #writing on a file
    f = open(outputfile, "w")
    f.writelines(text1)
    f.writelines(text2)
    f.writelines(text3)
    f.writelines(text4)
    f.close()

def semiglobalalignment(seq1, seq2, outputfile, scoringmatrix):

    #Dynamic Programming Solution

    lenv=len(seq1)
    lenw=len(seq2)
    s = [[0 for x in range(lenw+1)] for x in range(lenv+1)]  # +1 to account for first row and column which has gaps

    #Scoring Matrix
    gp=int(scoringmatrix[-1][0])
    print("gp: ", gp)

    
    for j in range(lenw+1): #seq w
        s[0][j]=0 #first row, no penalty for terminal gaps

    for i in range(lenv+1): #seq v
        s[i][0]= 0 #first column, no penalty for terminal gaps 
    
    for i in range(1,lenv+1):
        for j in range(1, lenw+1):
            max=-999
            cur=-999

            #Case 1 (Diagonal)
            idx1=scoringmatrix[0].index(seq1[i-1])
            idx2=scoringmatrix[0].index(seq2[j-1])
            if seq1[i-1]==seq2[j-1]:
                cur=s[i-1][j-1] + int(scoringmatrix[idx1+1][idx2+1]) #match score
    
            else:
                cur=s[i-1][j-1] + int(scoringmatrix[idx1+1][idx2+1]) #mismatch score
            
            if cur>max:
                max=cur
            
            #Case 2 (Vertical), vertical moves allowed in last column
            if j==lenw:
                cur=s[i-1][j] #no gap penalty
            else: 
                cur=s[i-1][j] + gp #gap penalty

            if cur>max:
                max=cur

            #Case 3 (Horizontal), horizontal moves allowed in bottom row
            if i==lenv:
                cur=s[i][j-1]  #no gap penalty 
            else: 
                cur= s[i][j-1] + gp #gap penalty

            if cur>max:
                max=cur

            #Putting the value in the table
            s[i][j]=max 
    print("--------------------------------------------------------------")
    #printing table
    for i in s:
        print('\t'.join(map(str, i)))
    print("--------------------------------------------------------------")

    #Trackback
    identity=0
    result= s[lenv][lenw]
    print("Result: ", result)
    print("--------------------------------------------------------------")
    sequence1=""
    sequence2=""

    i=lenv-1
    j=lenw-1
    
    while i>=0 and j>=0:
        print("i: ", i, "j: ", j)
        score=s[i+1][j+1] 
        print("Score: ", score) 

        #Diagonal Score
        idx1=scoringmatrix[0].index(seq1[i])
        idx2=scoringmatrix[0].index(seq2[j])
        if seq1[i]==seq2[j]:
            score1=s[i][j]+ int(scoringmatrix[idx1+1][idx2+1]) #match score
            identity=identity+1 
        else:
            score1=s[i][j] + int(scoringmatrix[idx1+1][idx2+1]) #mismatch score
        print("Score1: ", score1)

        #Vertical Score:
        if j==lenw-1:
            score2=s[i][j+1] #no gap penalty
        else: 
            score2= s[i][j+1] + gp #gap penalty 
        print("Score2: ", score2)

        #Horizontal Score:
        if i==lenv-1:
            score3=s[i+1][j]  #no gap penalty 
        else: 
            score3= s[i+1][j] + gp #gap penalty
        print("Score3: ",score3)
        
        if score==score1: #Diagonal
            sequence1=sequence1 + seq1[i]
            sequence2=sequence2 + seq2[j]
            i=i-1
            j=j-1
                   
        elif score==score2: #Vertical
            sequence1=sequence1+ seq1[i]
            sequence2=sequence2+ "_"
            i=i-1

        elif score==score3: #Horizontal
            sequence1=sequence1+"_"
            sequence2=sequence2+seq2[j]
            j=j-1
    
        print("=============================")

    text1= "seq1: ", str(i+1), " ", sequence1[::-1], " ", str(maximum(lenv, lenw)), "\n"
    text2= "seq2: ", str(j+1), " ", sequence2[::-1], " ", str(maximum(lenv, lenw)), "\n"
    text3= "Score: ", str(result), "\n"
    text4="Identities: ", str(identity),"/",str(maximum(lenv, lenw)), "  ", "(", str("{:.0%}".format(identity/maximum(lenv, lenw))), ")"

    #writing on a file
    f = open(outputfile, "w")
    f.writelines(text1)
    f.writelines(text2)
    f.writelines(text3)
    f.writelines(text4)
    f.close()

#Bonus Problem
def localalignment(seq1, seq2, outputfile, scoringmatrix):

    lenv=len(seq1)
    lenw=len(seq2)
    s = [[0 for x in range(lenw+1)] for x in range(lenv+1)]  # +1 to account for first row and column which has gaps
    highscore=-999999
    highesti=0
    highestj=0
    
    #Scoring Matrix
    gp=int(scoringmatrix[-1][0])
    print("gp: ", gp)
    
    for j in range(lenw+1): #seq w
        s[0][j]=0-j #first row

    for i in range(lenv+1): #seq v
        s[i][0]= 0-i #first column
    
    for i in range(1,lenv+1):
        for j in range(1, lenw+1):
            max=-99999
            cur=-99999

            #Case 1 (Diagonal)
            idx1=scoringmatrix[0].index(seq1[i-1])
            idx2=scoringmatrix[0].index(seq2[j-1])
            if seq1[i-1]==seq2[j-1]:
                cur=s[i-1][j-1] + int(scoringmatrix[idx1+1][idx2+1]) #match score
            else:
                cur=s[i-1][j-1] + int(scoringmatrix[idx1+1][idx2+1]) #mismatch score
            
            if cur>max:
                max=cur
            
            #Case 2 (Vertical)
            cur=s[i-1][j] + gp #gap penalty

            if cur>max:
                max=cur

            #Case 3 (Horizontal)
            cur= s[i][j-1] + gp #gap penalty

            if cur>max:
                max=cur

            #Putting the value in the table
            if max>=0:
                s[i][j]=max
            elif max<0:
                s[i][j]=0  

            if max>highscore:
                highscore=max
                highesti=i
                highestj=j
    
    #Trackback from highest score to 0
    result=highscore
    identity=0
    sequence1=""
    sequence2=""
    i=highesti
    j=highestj
    print(highesti, highestj)
    
    while s[i][j]!=0:
        score=s[i][j] 

        #Diagonal Score
        idx1=scoringmatrix[0].index(seq1[i-1])
        idx2=scoringmatrix[0].index(seq2[j-1])
        if seq1[i-1]==seq2[j-1]:
            score1=s[i-1][j-1]+ int(scoringmatrix[idx1+1][idx2+1]) #match score
            identity=identity+1 
        else:
            score1=s[i-1][j-1]+ int(scoringmatrix[idx1+1][idx2+1]) #mismatch score

        #Vertical Score:
        score2= s[i-1][j] + gp #gap penalty 
    
        #Horizontal Score:
        score3= s[i][j-1] + gp #gap penalty
        
        if score==score1: #Diagonal
            sequence1=sequence1 + seq1[i-1]
            sequence2=sequence2 + seq2[j-1]
            i=i-1
            j=j-1
            
        elif score==score2: #Vertical
            sequence1=sequence1+ seq1[i-1]
            sequence2=sequence2+ "_"
            i=i-1

        elif score==score3: #Horizontal
            sequence1=sequence1+"_"
            sequence2=sequence2+seq2[j-1]
            j=j-1

    text1= "seq1: ", str(i), " ", sequence1[::-1], " ", str(highesti), "\n"
    text2= "seq2: ", str(j), " ", sequence2[::-1], " ", str(highestj), "\n"
    text3= "Score: ", str(result), "\n"
    text4="Identities: ", str(identity),"/",str(maximum(lenv, lenw)), "  ", "(", str("{:.0%}".format(identity/maximum(lenv, lenw))), ")"

    #writing on a file
    f = open(outputfile, "w")
    f.writelines(text1)
    f.writelines(text2)
    f.writelines(text3)
    f.writelines(text4)
    f.close()
            
               

def main():
    seq1= loadSeq(sys.argv[1]) #this is v, iterated by i
    seq2= loadSeq(sys.argv[2]) #this is w, iterated by j 
    outputfile= sys.argv[3]
    inp=sys.argv[4] #T for Protein and F for DNA
    if inp=="F":
        scoringmatrix=matrix("dnaMatrix")
    elif inp=="T":
        scoringmatrix=matrix("BLOSUM45")
    else:
        print("Error in command line argument! should be in format: python3 align.py filename1.txt filename2.txt outputfilename.txt F/T G/S")
        return

    print("seq1 in FASTA format: ", seq1)
    print("seq2 in FASTA format: ", seq2)
    print("length of seq1: ", len(seq1))
    print("length of seq2: ", len(seq2))

    atype=sys.argv[5] #G for global and S for semi-global

    if atype=="G":
        globalalignment(seq1, seq2, outputfile, scoringmatrix)
    elif atype=="S":
        semiglobalalignment(seq1, seq2, outputfile, scoringmatrix)
    elif atype=="L":
        localalignment(seq1, seq2, outputfile, scoringmatrix)
    else:
        print("Error in command line argument! should be in format: python3 align.py filename1.txt filename2.txt outputfilename.txt F/T G/S")
        return

if __name__=="__main__":
    main()



    
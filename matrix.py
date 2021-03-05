def matrix(fileName):
    header=None
    scoringmatrix=list()
    file=open(fileName,'r')
    
    lines=file.readlines()
    
    for line in lines:
        if line.startswith('#'):
            header=line 
        else:
            line=line.split()
            scoringmatrix.append(line)

    #print(scoringmatrix)

    gp=scoringmatrix[-1]
    #print("Gap penalty: ", gp)
        
    file.close()

    return scoringmatrix

def main():
    print(matrix("dnaMatrix"))

if __name__=="__main__":
    main()

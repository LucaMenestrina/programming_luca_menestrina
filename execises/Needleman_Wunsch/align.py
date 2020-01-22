# two input sequences
seq1="ACTGG"
seq2="ACCA"

#AA or bases list
AA=["A","C","T","G"]

#gets the scoring matrix from a square matrix file
def get_scoring_matrix(file):
    matrix={}
    with open(file+".txt","r") as inputfile:
        AA=inputfile.readline().rstrip().split()
        n=0
        for line in inputfile:
            for i in range(len(AA)):
                matrix[AA[n]+AA[i]]=int(line.rstrip().split()[i+1])
                matrix[AA[i]+AA[n]]=int(line.rstrip().split()[i+1])
            n+=1
    return(matrix)

#aligns and scores two sequences on the basis of a specific scoring matrix with a defined gap penalty (linear)
def scorem(sequence1,sequence2,scoring_matrix=get_scoring_matrix("scoring_matrix"),d=-2):
    matrix=[[0]*(len(sequence1)+1) for i in range(len(sequence2)+1)]
    matrix_traceback=[[0]*(len(sequence1)+1) for i in range(len(sequence2)+1)]
    for y in range(len(sequence2)+1):
        for x in range(len(sequence1)+1):
            if y==0:
                matrix[0][x]=d*x
                matrix_traceback[0][x]="l"
            else:
                if x==0:
                    matrix[y][0]=d*y
                    matrix_traceback[y][0]="t"
                else:
                    matrix[y][x]=max(matrix[y-1][x]+d,matrix[y][x-1]+d,matrix[y-1][x-1]+get_scoring_matrix("scoring_matrix")[sequence2[y-1]+sequence1[x-1]])
                    if matrix[y][x] == matrix[y-1][x-1]+get_scoring_matrix("scoring_matrix")[sequence2[y-1]+sequence1[x-1]]:
                        matrix_traceback[y][x]="d"
                    elif matrix[y][x] == matrix[y-1][x]+d:
                        matrix_traceback[y][x]="t"
                    else:
                        matrix_traceback[y][x]="l"
    matrix_traceback[0][0]="O"
    return(matrix,matrix_traceback)

#partial test
#scorem(seq1,seq2)

# reconstruct the original sequences from the alignment matrix keeping in account gaps
def align(sequence1,sequence2):
    matrix=(scorem(sequence1,sequence2)[1])
    alignment=["",""]
    X=len(matrix[0])
    Y=len(matrix)
    while matrix[Y-1][X-1]!="O":
        direction=matrix[Y-1][X-1]
        if direction == "l":
            X-=1
            alignment[0]+=sequence1[X-1]
            alignment[1]+="-"
        elif direction == "t":
            Y-=1
            alignment[0]+="-"
            alignment[1]+=sequence2[Y-1]
        else:
            Y-=1
            X-=1
            alignment[0]+=sequence1[X-1]
            alignment[1]+=sequence2[Y-1]
    alignment[0]=alignment[0][::-1]
    alignment[1]=alignment[1][::-1]
    return(alignment)

#function test
align(seq1,seq2)

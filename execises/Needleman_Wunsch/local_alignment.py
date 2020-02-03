# two input sequences
seq1="GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKL"
seq2="GSGYLVGDSLTFVDLLVAQHTADLLAANAALLDEFPQFKAHQE"

#AA or bases list
AA=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","B","Z","X"]

#gets the scoring matrix from a square matrix file
def get_scoring_matrix(file="scoring_matrix_blosum"):
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
def scorem(sequence1,sequence2,scoring_matrix=get_scoring_matrix("scoring_matrix_blosum"),d=-2):
    # creates a scoring matrix
    matrix=[[0]*(len(sequence1)+1) for i in range(len(sequence2)+1)]
    # creates a parallel matrix for the traceback
    matrix_traceback=[[0]*(len(sequence1)+1) for i in range(len(sequence2)+1)]
    # matrix initialization y=0 or x=0, score = 0
    for y in range(len(sequence2)+1):
        for x in range(len(sequence1)+1):
            if y==0:
                matrix[0][x]=0
                matrix_traceback[0][x]="done"
            else:
                if x==0:
                    matrix[y][0]=0
                    matrix_traceback[y][0]="done"
                # fills both the matrices
                else:
                    matrix[y][x]=max(matrix[y-1][x]+d,matrix[y][x-1]+d,matrix[y-1][x-1]+scoring_matrix[sequence2[y-1]+sequence1[x-1]],0)
                    if matrix[y][x] == matrix[y-1][x-1]+scoring_matrix[sequence2[y-1]+sequence1[x-1]]:
                        matrix_traceback[y][x]="d"
                    elif matrix[y][x] == matrix[y-1][x]+d:
                        matrix_traceback[y][x]="t"
                    elif matrix[y][x] == 0:
                        matrix_traceback[y][x]="done"
                    else:
                        matrix_traceback[y][x]="l"
    matrix_traceback[0][0]="done"
    return(matrix,matrix_traceback)

##to visualize matrices
#for matr in [0,1]:
#    for line in range(len(seq2)+1):
#        print(scorem(seq1,seq2)[matr][line])

#partial test
#scorem(seq1,seq2)

# reconstruct the original sequences from the alignment matrix keeping in account gaps
def align(sequence1,sequence2,threshold):
    alignments=[]
    matrix_traceback=scorem(sequence1,sequence2)[1]
    matrix=scorem(sequence1,sequence2)[0]
    scores_dict={}
    for row in range(len(matrix)):
        for col in range(len(matrix[0])):
            scores_dict[row,col]=matrix[row][col]
    for starting_pos in scores_dict.keys():
        alignment=[0,"",""]
        if scores_dict[starting_pos]>=threshold:
            Y=starting_pos[0]+1
            X=starting_pos[1]+1
            score=matrix[Y-1][X-1]
            if matrix_traceback[Y-1][X-1] == "d":
                while matrix_traceback[Y-1][X-1]!="done":
                    direction=matrix_traceback[Y-1][X-1]
                    if direction == "l":
                        X-=1
                        alignment[1]+=sequence1[X-1]
                        alignment[2]+="-"
                    elif direction == "t":
                        Y-=1
                        alignment[1]+="-"
                        alignment[2]+=sequence2[Y-1]
                    elif direction == "d":
                        Y-=1
                        X-=1
                        alignment[1]+=sequence1[X-1]
                        alignment[2]+=sequence2[Y-1]
                alignment[1]=alignment[1][::-1]
                alignment[2]=alignment[2][::-1]
                alignment[0]=score
                alignments.append(alignment)
    # in order to sort all the alignments on the basis of their scores
    alignments.sort(reverse=True)
    return(alignments)


##function test
#align(seq1,seq2,4)

# simply to print in a nice way the results
def format_print(sequence1,sequence2,threshold):
    for allineamento in range(len(align(sequence1,sequence2,threshold))):
        for sequenza in [1,2]:
            print(align(sequence1,sequence2,threshold)[allineamento][sequenza])
        print("Score: %d\n" %align(sequence1,sequence2,threshold)[allineamento][0])


#----------------------------------------------------------------------------------------------------------------------

## example test
format_print(seq1,seq2,57)

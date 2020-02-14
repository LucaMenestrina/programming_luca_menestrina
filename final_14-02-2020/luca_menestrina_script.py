import input_data

def calc_matrices(sequence1,sequence2,scoring_matrix,gap):
    F=[[0]*(len(sequence1)+1) for i in range(len(sequence2)+1)]
    P=[[0]*(len(sequence1)+1) for i in range(len(sequence2)+1)]
    for y in range(len(sequence2)+1):
        for x in range(len(sequence1)+1):
            if y==0:
                F[0][x]=gap*x
                P[0][x]="l"
            elif x==0:
                F[y][0]=gap*y
                P[y][0]="t"
            else:
                sd=F[y-1][x-1]+scoring_matrix[sequence1[x-1]+sequence2[y-1]]
                st=F[y-1][x]+gap
                sl=F[y][x-1]+gap
                max_value=max(sd,st,sl)
                F[y][x]=max_value
                if max_value==sd:
                    P[y][x]="d"
                elif max_value==st:
                    P[y][x]="t"
                elif max_value==sl:
                    P[y][x]="l"
    P[0][0]="done"
    return(F,P)

def traceback(F,P,sequence1,sequence2):
    aligned1=""
    aligned2=""
    X=len(sequence1)
    Y=len(sequence2)
    score=F[Y][X]
    while P[Y][X] != "done":
        if P[Y][X]=="t":
            Y-=1
            aligned1+="-"
            aligned2+=sequence2[Y]
        elif P[Y][X]=="l":
            X-=1
            aligned1+=sequence1[X]
            aligned2+="-"
        elif P[Y][X]=="d":
            X-=1
            Y-=1
            aligned1+=sequence1[X]
            aligned2+=sequence2[Y]
    aligned1=aligned1[::-1]
    aligned2=aligned2[::-1]
    return(aligned1,aligned2,score)

print(traceback(calc_matrices(input_data.seq1,input_data.seq2,input_data.BLOSUM52,-2)[0],calc_matrices(input_data.seq1,input_data.seq2,input_data.BLOSUM52,-2)[1],input_data.seq1,input_data.seq2)[0])
print(traceback(calc_matrices(input_data.seq1,input_data.seq2,input_data.BLOSUM52,-2)[0],calc_matrices(input_data.seq1,input_data.seq2,input_data.BLOSUM52,-2)[1],input_data.seq1,input_data.seq2)[1])
print("Score: "+str(traceback(calc_matrices(input_data.seq1,input_data.seq2,input_data.BLOSUM52,-2)[0],calc_matrices(input_data.seq1,input_data.seq2,input_data.BLOSUM52,-2)[1],input_data.seq1,input_data.seq2)[2]))

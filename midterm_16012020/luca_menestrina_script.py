# MENESTRINA LUCA
# Mid term exam programming for bioinformatics module 2
# Wrote for python 3

# Gets the alignments from the fasta file
def get_fasta(fasta_file):
    general_list=[]
    with open (fasta_file,"r") as alignments:
        for line in alignments:
            if ">" not in line:
                general_list.append(line)
    separated_list=[]
    for item in range(int(len(general_list)/2)):
        separated_list.append([])
    n=0
    i=0
    for item in range(len(general_list)):
        separated_list[i].append(general_list[item].rstrip())
        if n==1:
            i+=1
            n=0
        else:
            n+=1
    return(separated_list)

# Gets the scoring_matrix
# It requires the name as input as "PAM250" or "BLOSUM62" WITH quotes (it needs a string)
def get_scoring_matrix(scoring_matrix):
    dict={}
    with open (scoring_matrix+".txt","r") as matrix:
        AA=matrix.readline().split()[3]
        n=0
        for line in matrix:
            for i in range(n+1):
                dict[AA[i]+AA[n]]=int(line.split()[i][:-1])
                dict[AA[n]+AA[i]]=int(line.split()[i][:-1])
            n+=1
    return(dict)

# Scores two input sequences with the assigned scoring matrix
def scoring(sequence1,sequence2,scoring_matrix):
    score=0
    for i in range(len(sequence1)):
        if sequence1[i]=="-" or sequence2[i]=="-":
            score+=-2
        else:
            score+=get_scoring_matrix(scoring_matrix)[sequence1[i]+sequence2[i]]
    return(sequence1,sequence2,scoring_matrix,score)

# Scores two input sequences with the assigned scoring matrix with affine gap penalty
def scoring_affine(sequence1,sequence2,scoring_matrix):
    gap_length=1
    score=0
    for i in range(len(sequence1)):
        if sequence1[i]=="-" or sequence2[i]=="-":
            if sequence1[i+1]=="-" or sequence2[i+1]=="-":
                gap_length+=1
            else:
                score+=(-2-((gap_length-1)*0.5))
                gap_length=1
        else:
            score+=get_scoring_matrix(scoring_matrix)[sequence1[i]+sequence2[i]]
    return(sequence1,sequence2,scoring_matrix,score)

scoring_affine(get_fasta("alignments.fasta")[1][0],get_fasta("alignments.fasta")[1][1],"BLOSUM62")

# Iterates through the alignments and the scoring matrices and prints a results
def score_everything(fasta_file):
    matrices_list=["BLOSUM62","PAM250"]
    for alignments in range(len(get_fasta(fasta_file))):
        for n in range(2):
            if n==0:
                print("\nAlignment"+str(alignments+1)+":\n"+"First sequence:  "+scoring(get_fasta(fasta_file)[alignments][0],get_fasta(fasta_file)[alignments][1],matrices_list[n])[0]+"\nSecond sequence: "+scoring(get_fasta(fasta_file)[alignments][0],get_fasta(fasta_file)[alignments][1],matrices_list[n])[1])
            print("With "+matrices_list[n]+" the score is: "+str(scoring(get_fasta(fasta_file)[alignments][0],get_fasta(fasta_file)[alignments][1],matrices_list[n])[3]))

def score_everything_affine(fasta_file):
    matrices_list=["BLOSUM62","PAM250"]
    for alignments in range(len(get_fasta(fasta_file))):
        for n in range(2):
            if n==0:
                print("\nAlignment"+str(alignments+1)+":\n"+"First sequence:  "+scoring_affine(get_fasta(fasta_file)[alignments][0],get_fasta(fasta_file)[alignments][1],matrices_list[n])[0]+"\nSecond sequence: "+scoring_affine(get_fasta(fasta_file)[alignments][0],get_fasta(fasta_file)[alignments][1],matrices_list[n])[1])
            print("With "+matrices_list[n]+" the score is: "+str(scoring_affine(get_fasta(fasta_file)[alignments][0],get_fasta(fasta_file)[alignments][1],matrices_list[n])[3]))

score_everything("alignments.fasta")

score_everything_affine("alignments.fasta")

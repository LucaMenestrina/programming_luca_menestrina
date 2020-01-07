import math

#sequenza1="TCA"
#sequenza2="GA"

# creates a dictionary to store blosum data and fills it (it needs a file calles blosum.txt with the blosum scores in the same folder)
blosum={}
with open("blosum.txt","r") as data_structure:
    # creates a list of all AA
    AA=data_structure.readline().split()
    # loops through every line in the file and every AA in the list, creating every possible couple
    for line in data_structure:
        for i in range(len(AA)):
            blosum[line.split()[0]+AA[i]]=int(line.split()[i+1])

bases_score={}
with open("base_score.txt","r") as data_structure:
    # creates a list of all AA
    AA=data_structure.readline().split()
    # loops through every line in the file and every AA in the list, creating every possible couple
    for line in data_structure:
        for i in range(len(AA)):
            bases_score[line.split()[0]+AA[i]]=int(line.split()[i+1])

# test sequences
sequenza1="ALAS"
sequenza2="AS"

# defines a scoring function that returns the sum of the blosum values for all the couples in the two strings
def score(seq1,seq2,data_structure):
    score=0
    for base in range(len(seq1)):
        score+=data_structure[seq1[base]+seq2[base]]
    return score



def align_and_score(seq1,seq2,data_structure):
    length=len(seq1+seq2)
    possibilities={}
    scores={}
    for i in range(length+1):
        possibilities[i]=["-"*(len(seq2)-i)+seq1+"-"*(i-(len(seq1)-len(seq2))-len(seq2)),"-"*(i+(len(seq1)-len(seq2))-len(seq1))+seq2+"-"*(len(seq1)-i)]
    for i in range(len(possibilities)):
        scores[i]=score(possibilities[i][0],possibilities[i][1],data_structure)
    print("Best alignment: "+str(possibilities[max(scores,key=scores.get)])+"\nWith a score of: "+str(scores[max(scores,key=scores.get)]))
    #return(possibilities[max(scores,key=scores.get)],scores[max(scores,key=scores.get)])

align_and_score(sequenza1,sequenza2,blosum)

sequenza3="GA"
sequenza4="TCA"
align_and_score(sequenza3,sequenza4,bases_score)

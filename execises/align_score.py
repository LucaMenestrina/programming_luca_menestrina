import math

# creates a dictionary to store blosum data and fills it (it needs a file calles blosum.txt with the blosum scores in the same folder)
blosum={}
with open("blosum.txt","r") as data_structure:
    # creates a list of all AA
    AA=data_structure.readline().split()
    # loops through every line in the file and every AA in the list, creating every possible couple
    for line in data_structure:
        for i in range(len(AA)):
            blosum[line.split()[0]+AA[i]]=int(line.split()[i+1])

# test sequences
sequenza1="ALASVLIRLITRLYP"
sequenza2="ASAVHLNRLITRLYP"

# defines a scoring function that returns the sum of the blosum values for all the couples in the two strings
def score(seq1,seq2):
    score=0
    for base in range(len(seq1)):
        score+=blosum[seq1[base]+seq2[base]]
    return score

#test function
score(sequenza1,sequenza2)

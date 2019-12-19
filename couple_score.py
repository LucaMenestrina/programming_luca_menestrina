import math

sequenza1="ACAGGTGGACCT"
sequenza2="ACTGGTCGACTT"

# function that calculates the score for a specific base couple
def score(seq1,seq2,base1,base2):
    s=seq1+seq2
    base_dict={}
    couple_dict={}
    # counts all bases in the sequences
    for base in s:
        base_dict[base]=base_dict.get(base,0)+1
    # updates the dictionary with the frequencies of the bases
    for key in base_dict:
        base_dict[key]=base_dict[key]/len(s)
    # counts all the bases' couples (it creates a dictionary only for the existing couples, it is not general)
    # every time it counts a couple it updates also the value of the opposite couple
    for i in range(len(seq1)):
        couple_dict[seq1[i]+seq2[i]]=couple_dict.get(seq1[i]+seq2[i],0)+1
        couple_dict[seq2[i]+seq1[i]]=couple_dict.get(seq2[i]+seq1[i],0)+(seq1[i]!=seq2[i])
    # updates the dictionary with the frequences of the couples
    for key in couple_dict:
        couple_dict[key]/=len(seq1)
    # computes the score for the input couple
    score=math.log(couple_dict[base1+base2]/(base_dict[base1]*base_dict[base2]))
    # returns the dictionary with the bases' frequencies, the one with the couples' frequencies and the score for the input couple
    return base_dict,couple_dict,score

# example with the two input sequences and AT as couple (the same as TA)
score(sequenza1,sequenza2, "A","T")

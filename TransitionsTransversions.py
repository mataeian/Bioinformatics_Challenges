# Problem
# For DNA strings s1 and s2 having the same length, their transition/transversion ratio R(s1,s2) is the ratio of the total number of transitions to the total number of transversions, where symbol substitutions are inferred from mismatched corresponding symbols as when calculating Hamming distance (see “Counting Point Mutations”).

# Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).

# Return: The transition/transversion ratio R(s1,s2)

#to solve:
1. read the fasta file
2. compare both strings
3. count transistions
4. count transversions
5. calculate ratio

from Bio import SeqIO

def Rratio(s1, s2):
    transitions = 0
    transversions = 0
    for i in range(len(s1)):
        if 



def main():
    fasta_file = "rosalind_sseq.txt"
    sequences = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))
    s = sequences[0]
    t = sequences[1]
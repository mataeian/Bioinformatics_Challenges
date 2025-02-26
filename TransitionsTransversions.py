# Problem
# For DNA strings s1 and s2 having the same length, their transition/transversion ratio R(s1,s2) is the ratio of the total number of transitions to the total number of transversions, where symbol substitutions are inferred from mismatched corresponding symbols as when calculating Hamming distance (see “Counting Point Mutations”).

# Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).

# Return: The transition/transversion ratio R(s1,s2)

#to solve:
# 1. read the fasta file
# 2. compare both strings
# 3. count transistions
# 4. count transversions
# 5. calculate ratio
# 6. return results

from Bio import SeqIO

def Rratio(s1, s2):
    transitions = 0
    transversions = 0
    
    transition_pairs = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
    transversion_pairs = [('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'),
                         ('G', 'C'), ('C', 'G'), ('G', 'T'), ('T', 'G')]
    
    for i in range(len(s1)):
        if s1[i]== s2[i]:
            continue
        elif (s1[i], s2[i]) in transition_pairs:
            transitions += 1
        elif (s1[i], s2[i]) in transversion_pairs:
            transversions += 1
    
    if transversions == 0:
        return 0
    return transitions/transversions



def main():
    fasta_file = "rosalind_tran.txt"
    sequences = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))
    s1 = sequences[0]
    s2 = sequences[1]

    ratio = Rratio(s1,s2)
    print(ratio)

if __name__ == "__main__":
    main()
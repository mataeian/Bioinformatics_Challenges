# Problem
# A subsequence of a string is a collection of symbols contained in order (though not necessarily contiguously) in the string (e.g., ACG is a subsequence of TATGCTAAGATC). The indices of a subsequence are the positions in the string at which the symbols of the subsequence appear; thus, the indices of ACG in TATGCTAAGATC can be represented by (2, 5, 9).

# As a substring can have multiple locations, a subsequence can have multiple collections of indices, and the same index can be reused in more than one appearance of the subsequence; for example, ACG is a subsequence of AACCGGTT in 8 different ways.

# Given: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.

# Return: One collection of indices of s in which the symbols of t appear as a subsequence of s. If multiple solutions exist, you may return any one.

# To solve:
# Read the FASTA file
# Find the indices
# Return the indices

from Bio import SeqIO

            
def find_subsequence_indices (s, t):
    indeces = []
    t_index = 0
    for i in range(len(s)):
        if t_index >= len(t):
            break
        if s[i] == t[t_index]:
            indeces.append(i + 1)  
            t_index += 1
    return indeces


def main():
    fasta_file = "rosalind_sseq.txt"
    sequences = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))
    s = sequences[0]
    t = sequences[1]

indeces = find_subsequence_indices (s,t)
print(" ".join(map(str, indices)))

if __name__ == "__main__":
    main()

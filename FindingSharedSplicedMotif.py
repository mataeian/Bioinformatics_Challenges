# Problem
# A string u is a common subsequence of strings s and t if the symbols of u appear in order as a subsequence of both s and t. For example, "ACTG" is a common subsequence of "AACCTTGG" and "ACACTGTGA".

# Analogously to the definition of longest common substring, u is a longest common subsequence of s and t if there does not exist a longer common subsequence of the two strings. Continuing our above example, "ACCTTG" is a longest common subsequence of "AACCTTGG" and "ACACTGTGA", as is "AACTGG".

# Given: Two DNA strings s and t (each having length at most 1 kbp) in FASTA format.

# Return: A longest common subsequence of s and t. (If more than one solution exists, you may return any one.)

#To solve:
# 1. read the 2 sequences
# 2. compare the 2 strings
# 3. collect shared motifs
# 4. report the longest



from Bio import SeqIO
from Bio.Seq import Seq








def main():
    fasta_file = "input.fasta"  
    sequences = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))
    
    s = sequences[0]  
    t = sequences[1]  

    lcs = longest_common_subsequence(s, t)

    print(lcs)

if __name__ == "__main__":
    main()
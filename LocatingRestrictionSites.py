# Problem

# A DNA string is a reverse palindrome if it is equal to its reverse complement. For instance, GCATGC is a reverse palindrome because its reverse complement is GCATGC. See Figure 2.

# Given: A DNA string of length at most 1 kbp in FASTA format.

# Return: The position and length of every reverse palindrome in the string having length between 4 and 12. You may return these pairs in any order.

# To solve:
# Read the DNA sequence
# Generate reverse complements
# Check for palindromes
# Record positions and lengths




from Bio import SeqIO
from Bio.Seq import Seq


def find_reverse_palindromes(dna, min_length=4, max_length=12):
    results = []
    length = len(dna)
    
    for l in range(min_length, max_length + 1):
        for i in range(length - l + 1):
            substring = dna[i:i+l]
            if substring == str(Seq(substring).reverse_complement()):
                results.append((i + 1, l))  # +1 to convert to 1-based indexing
    
    return results


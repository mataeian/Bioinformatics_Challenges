# Problem For a fixed positive integer k, order all possible k-mers taken from an underlying alphabet lexicographically.

# Then the k-mer composition of a string s can be represented by an array A for which A[m] denotes the number of times that the mth k-mer (with respect to the lexicographic order) appears in s.

# Given: A DNA string s in FASTA format (having length at most 100 kbp).

# Return: The 4-mer composition of s.

# To solve:
# for DNA the alphabet is {"A", "C","G","T"}, so there are 4k possible k=mers. for k=4, there are 44 =256 possible 4mers.
# 1. generate all possible 4-mers 
# 2. order lexicographically
# 3. count frequencies (m)
# 4. list m

from itertools import product # It generates all possible combinations of elements from the input iterables
from collections import defaultdict #it provides a default value for a nonexistent key, so that the program doesn't raise a KeyError

def generate_kmers(k):
    """Generate all possible k-mers in lexicographical order."""
    alphabet = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in product(alphabet, repeat=k)]


def count_kmers(s,k):
    """Count the frequency of each k-mer in the string s."""
    kmer_counts = defaultdict(int)
    for i in range(len(s) - k + 1):
        kmer = s[i:i+k]
        kmer_counts[kmer] += 1
    return kmer_counts

def main():
    fasta_file = "rosalind_kmer.txt"  
    with open(fasta_file, "r") as handle:
        # Skip the header line
        next(handle)
        # Read the DNA sequence
        s = ''.join(line.strip() for line in handle)

    # Generate all possible 4-mers
    kmers = generate_kmers(4)

    # Count the frequency of each 4-mer in the DNA sequence
    kmer_counts = count_kmers(s, 4)

    # Output the 4-mer composition in lexicographical order
    composition = [kmer_counts[kmer] for kmer in kmers]
    print(' '.join(map(str, composition)))
   
   
if __name__ == "__main__":
    main()
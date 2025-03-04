# Problem For a fixed positive integer k, order all possible k-mers taken from an underlying alphabet lexicographically.

# Then the k-mer composition of a string s can be represented by an array A for which A[m] denotes the number of times that the mth k-mer (with respect to the lexicographic order) appears in s.

# Given: A DNA string s in FASTA format (having length at most 100 kbp).

# Return: The 4-mer composition of s.
To solve:
for DNA the alphabet is {"A", "C","G","T"}, so there are 4k possible k=mers. for k=4, there are 44 =256 possible 4mers.
find 4-mers
order lexicographically
count A and number as m
list m

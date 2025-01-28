# Problem
# Assume that an alphabet A has a predetermined order; that is, we write the alphabet as a permutation A=(a1,a2,…,ak)
# , where a1<a2<⋯<ak . For instance, the English alphabet is organized as (A,B,…,Z)

# Given two strings s and t having the same length n, we say that s precedes t in the lexicographic order (and write s<Lext) if the first symbol s[j] that doesn't match t[j] satisfies sj<tj in A.

# Given: A collection of at most 10 symbols defining an ordered alphabet, and a positive integer n (n≤10).

# Return: All strings of length n that can be formed from the alphabet, ordered lexicographically (use the standard order of symbols in the English alphabet).

import itertools

def generate_lexicographic_strings(alphabet, n):
    all_strings = itertools.product(alphabet, repeat=n)
    
    result = ["".join(string) for string in all_strings]
    return result

alphabet = ["A", "B","C", "D", "E", "F", "G", "H", "I", "J"]
n  = 2

lexicographic_strings = generate_lexicographic_strings(alphabet, n)

for string in lexicographic_strings:
    print(string)
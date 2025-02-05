# Problem
# A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGTATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".

# Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

# Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.

# Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)


def parse_fasta(file):
    sequences = []
    with open ("rosalind_lcsm.txt", "r") as f:
        seq = ""
        for line in f:
            if line.startswith(">"):  
                if seq:  
                    sequences.append(seq)
                seq = ""  
            else:
                seq += line.strip()  
        sequences.append(seq) 
    return sequences


def longest_common_substring(dna_strings):
    # Sort the strings by length to start with the shortest one
    dna_strings.sort(key=len)
    shortest_str = dna_strings[0]
    other_strings = dna_strings[1:]
    longest_common = ""

    # Check all substrings of the shortest string
    for i in range(len(shortest_str)):
        for j in range(i + 1, len(shortest_str) + 1):
            candidate = shortest_str[i:j]
            if all(candidate in dna for dna in other_strings): #check if the substring exists in other DNA strings
                if len(candidate) > len(longest_common):
                    longest_common = candidate

    return longest_common

dna_strings = parse_fasta("rosalind_lcsm.txt")
result = longest_common_substring(dna_strings)
print(result)
    
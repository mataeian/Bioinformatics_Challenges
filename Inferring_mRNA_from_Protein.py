# Problem
# For positive integers a and n, a modulo n (written amodn in shorthand) is the remainder when a is divided by n. For example, 29mod11=7
#  because 29=11×2+7.
# Modular arithmetic is the study of addition, subtraction, multiplication, and division with respect to the modulo operation. We say that a and b are congruent modulo n if amodn=bmodn; in this case, we use the notation a≡bmodn.

# Given: A protein string of length at most 1000 aa.

# Return: The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)

with open("rosalind_mrna.txt", "r") as file:
    protein_string= file.read().strip()

codon_counts = {
    'F': 2, 'L': 6, 'S': 6, 'Y': 2, 'C': 2, 'W': 1, 'P': 4, 'H': 2,
    'Q': 2, 'R': 6, 'I': 3, 'M': 1, 'T': 4, 'N': 2, 'K': 2, 'V': 4,
    'A': 4, 'D': 2, 'E': 2, 'G': 4, 'Stop': 3
}

def rna_possibilities(protein_string):
    total_possibilities = 1
    for amino_acid in protein_string:
        total_possibilities *= codon_counts[amino_acid]
        total_possibilities %= 1000000

    total_possibilities*= codon_counts["Stop"]
    total_possibilities %= 1000000

    return total_possibilities

result = rna_possibilities(protein_string)
print(result)
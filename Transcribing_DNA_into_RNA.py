# Problem
# An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.
# Given a DNA string to corresponding to a coding strand, its transcribed RNA string u is formed by replacing all occurrences of 'T' in t with 'U' in u
# Return: The transcribed RNA string of t
with open("rosalind_rna.txt", "r") as file:
    t = file.read()

u= t.replace("T", "U")

with open("rosalind_rna_output.txt", "w") as output:
    output.write(u)

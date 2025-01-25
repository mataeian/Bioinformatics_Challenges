#problem
#The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s , then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

#Given: A DNA string s of length at most 1000 bp.

#Return: The reverse complement sc of s

with open("rosalind_revc.txt", "r") as file:
    s=file.read().strip() #strip removes any extra spaces in the file

complement = {"A":"T", "T":"A","C":"G", "G":"C"}

sc = "".join(complement[base] for base in reversed(s))

with open("rosalind_revc_output.txt", "w") as output:
    output.write(sc)

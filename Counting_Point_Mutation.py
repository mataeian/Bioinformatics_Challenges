#Problem
# Given two strings s and t of equal length, the Hamming distance between s and t , denoted dH(s,t), is the number of corresponding symbols that differ in s and t
# Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).
# Return: The Hamming distance dH(s,t)

with open("rosalind_hamm.txt", "r") as file:
    s = file.readline().strip()
    t = file.readline().strip()

def calculate_Hamming_distance (s,t):
    dh = 0
    for chr in range(len(s)):
        if s[chr] != t[chr]:
            dh += 1
    return dh

result = calculate_Hamming_distance(s, t)
print(result)



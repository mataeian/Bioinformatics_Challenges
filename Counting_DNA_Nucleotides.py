#problem:
#Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s

s = "TGTGTCCTACATTGAGGTTCTAGTTCTCCCTTTGCATCTGGCGTGTTGGTTATCGGCAACCATTTTGGTTCCGTGAGGCTGCATGACCTGAGTGATAAGCACGTAAGTTGACGGTCCCAACCGAAGCGTCTCTAGCTCCCTCCCTCTAGTAATGTAGGCACCTCGTGTCACGGTGCGGTAAGCAGTCTTGGTTCGATGTTTTCCGAATTAGACTGCAGTATTTACCACCTTCGCGGGATTATCACCGACAAATCTCGGAGTCGAAACTGTACTGGGCAGTTACTTCAAACCGTGGAAGACGGGAGTACGCTCCTGAAGGGGCTACGATCCCGTTACGATTCGGACTATTGCCTCAAGAATATTTGGCCAATGGCTCTGAGCCACACGGAAAGATGAGATATTGCGATAGTGGGAAGACGCGTTTCCAGTAGGCAAAACTACGAGAAGTACGGCAAAGCTTGGTCAACTCAAAACAGACGAGAACGCACGCGGCAGGGAAAGTGCCCCCCTAGGACGACACAGTAAGGCGGTTACGGGTGTCCGAATGGGACTCACAATTTATATATTGATCCTCGCCTGGATTGATTAGCTTCGCCCCAACGAGCTCGCGTCGATGTCCAACTGTATGCACGTGAGGTCGCGAAATAAAGGTGACACACCCGGCCCTGCCACGTTTAGCACCACTACGGCACTTATCGGCAGCTGTCACCGAAAGTGACGTAATCAGGTGATGTGCCCTCTATCCGCAGCAGTAATGGCACTAGTATAAGCCCGGAAGGAATTTCCGCTAACTTGCTCCCTTAGAC"
count_A= s.count("A")
count_C= s.count("C")
count_G= s.count("G")
count_T= s.count("T")

counts = (count_A, count_C, count_G, count_T)
print(counts)


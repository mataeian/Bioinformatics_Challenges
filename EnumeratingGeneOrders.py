#Problem
# A permutation of length n is an ordering of the positive integers {1,2,…,n}. For example, π=(5,3,2,1,4) is a permutation of length 5
# Given: A positive integer n≤7.
# Return: The total number of permutations of length n, followed by a list of all such permutations (in any order).

import itertools

def generate_permutations(n):
    numbers = list(range(1, n + 1))
    
    all_permutations = list(itertools.permutations(numbers))
    
    total_permutations = len(all_permutations)
    print(total_permutations)
    
    with open("permutation_output", "w") as file:
        file.write(f"{total_permutations}\n")
        for perm in all_permutations:
            file.write(" ".join(map(str, perm)) + "\n") #The map() function applies the str function to each element of the perm tuple, converting each number to a string.


generate_permutations(5)



# Problem
# A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences can be finite or infinite. Two examples are the finite sequence (π,−2–√,0,π) and the infinite sequence of odd numbers (1,3,5,7,9,…). We use the notation an to represent the n-th term of a sequence.
# A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal to the number of rabbits that were alive two months prior. As a result, if Fn represents the number of rabbit pairs alive after the n
# -th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence relation Fn=Fn−1+Fn−2 (with F1=F2=1 to initiate the sequence). Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.
# When finding the n-th term of a sequence defined by a recurrence relation, we can simply use the recurrence relation to generate terms for progressively larger values of n. This problem introduces us to the computational technique of dynamic programming, which successively builds up solutions by using the answers to smaller cases.

# Given: Positive integers n≤40 and k≤5.

# Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).

# To solve:
# 1. Initialization: Start with 1 pair of rabbits in the first month.

# 2. Recurrence Relation: For each subsequent month, the number of rabbit pairs is the sum of the rabbit pairs from the previous month and the rabbit pairs from two months ago multiplied by k (since each pair of reproduction-age rabbits produces k pairs of offspring).

# 3. Base Cases: The first two months have 1 pair of rabbits each.

# 4. Iteration: Use a loop to compute the number of rabbit pairs for each month up to n.

def rabbit_pairs(n, k):
    # Initialize the base cases
    if n == 1 or n == 2:
        return 1
    
    # Initialize an array to store the number of rabbit pairs for each month
    dp = [0] * (n + 1)
    dp[1] = 1
    dp[2] = 1
    
    # Fill the array using the recurrence relation
    for i in range(3, n + 1):
        dp[i] = dp[i - 1] + k * dp[i - 2]
    
    # Return the number of rabbit pairs after n months
    return dp[n]

n = 31
k = 2

result = rabbit_pairs(n, k)

print(result)  

#helper method that returns if the bases are pairs
def is_pair(base1, base2):    
    return (base1 == 'A' and base2 == 'U') or (base1 == 'U' and base2 == 'A') or (base1 == 'G' and base2 == 'C') or (base1 == 'C' and base2 == 'G') or (base1 == 'G' and base2 == 'U') or (base1 == 'U' and base2 == 'G')


# rna folding algorithm to find the maximum number of pairs with a minimum pairing distance which is set to 4
def build_table(seq, min_distance=4):
    n = len(seq) # length of the input sequence
    
    #initialize the dp table - fill with 0's
    dp = [[0] * n for _ in range(n)]
    
    #now fill in the dp table
    for length in range(2, n + 1):  # length of subsequence (starting from length 2)
        for i in range(n - length + 1):
            j = i + length - 1

            # Option 1: If xi and xj form a valid pair and j - i >= min_distance
            if is_pair(seq[i], seq[j]) and (j - i >= min_distance):
                option1 = dp[i + 1][j - 1] + 1 
            else: option1 = 0

            # Option 2: Skip xj (take the value from dp(i + 1, j))
            option2 = dp[i + 1][j]

            # Option 3: Skip xi (take the value from dp(i, j - 1))
            option3 = dp[i][j - 1]

            # Option 4: Split at every possible k and calculate the sum of bp(i, k) and bp(k + 1, j)
            option4 = 0
            for k in range(i, j):
                option4 = max(option4, dp[i][k] + dp[k + 1][j])

            # Store the maximum of the 4 options
            dp[i][j] = max(option1, option2, option3, option4)

    return dp



#traceback function to get the optimal structure
def traceback(dp, seq, i, j, structure, min_distance=4):
    if i < j:
        if dp[i][j] == dp[i][j - 1]:  # No pairing at position j
            traceback(dp, seq, i, j - 1, structure, min_distance)
        elif is_pair(seq[i], seq[j]) and dp[i][j] == dp[i + 1][j - 1] + 1 and (j - i >= min_distance):  # i and j are paired
            structure[i] = "("
            structure[j] = ")"
            traceback(dp, seq, i + 1, j - 1, structure, min_distance)
        else:
            # Try splitting at some k
            for k in range(i, j):
                if dp[i][j] == dp[i][k] + dp[k + 1][j]:
                    traceback(dp, seq, i, k, structure, min_distance)
                    traceback(dp, seq, k + 1, j, structure, min_distance)
                    break





# main function to run the algorithm and print the structure
def predict_rna_structure(seq, min_distance=4):
    n = len(seq)
    dp = build_table(seq, min_distance) 
    for rows in dp:
        print(rows)
    structure = ['.'] * n  #initialize the structure - unpaired bases - "."
    traceback(dp, seq, 0, n - 1, structure, min_distance)
    final_structure = "".join(structure)
    num_pairs = dp[0][n - 1]
    return final_structure, num_pairs


# input RNA sequence for which we want to compute the optimal base pairing structure

sequence = 'GUGAGAC'
#call to function to run program on input sequence
structure, max_pairs = predict_rna_structure(sequence)
print(f"RNA Sequence:      {sequence}")
print(f"Optimal Structure: {structure}")
print(f"Max Number of Pairs: {max_pairs}")

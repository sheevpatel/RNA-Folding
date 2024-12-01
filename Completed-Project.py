# This is a Helper method to check ONLY if a pair is valid
def is_pair(base1, base2):
    return (base1 == 'A' and base2 == 'U') or (base1 == 'U' and base2 == 'A') or \
           (base1 == 'C' and base2 == 'G') or (base1 == 'G' and base2 == 'C') or \
           (base1 == 'G' and base2 == 'U') or (base1 == 'U' and base2 == 'G')

# This is a helper method that returns the score of the base pair
def get_pair_score(base1, base2):
    if (base1 == 'A' and base2 == 'U') or (base1 == 'U' and base2 == 'A'):
        # I let Watson-Crick Pair C-G = 3
        return 2
    elif (base1 == 'C' and base2 == 'G') or (base1 == 'G' and base2 == 'C'):
        # I let Watson-Crick Pair A-U = 2
        return 3
    elif (base1 == 'G' and base2 == 'U') or (base1 == 'U' and base2 == 'G'):
        # I let Watson-Crick Pair G-U = 1
        return 1
    return 0

# DP Table Building Algorithm to maximize pairs and points, with priority to pairs
# This table is 3D - it really is 2D but each cell in the 2D table is an array of length 2
# The 1st value in the array is the max number of pairs for the subsequence
# The 2nd value in the array is the point total from the respective pairs for that subsequence
# min_distance is the minimum distance parameter
def build_table(seq, min_distance=4):
    n = len(seq)  # length of the input sequence
    
    # initialize the 3D dp table - fill with [0,0] for each cell
    dp = [[[0, 0] for _ in range(n)] for _ in range(n)]
    
    # now fill in the dp table
    for length in range(2, n + 1):  # length of subsequence starts at 2, n+1 bc range is exclusive for parameter 2
        for i in range(n - length + 1): 
            j = i + length - 1

            # Objective Function is implemented to get max base pairs
            # Points are accumulated as well for the respective base pairs

            # option 1: if xi and xj form a valid pair and j - i >= min_distance
            if is_pair(seq[i], seq[j]) and (j - i >= min_distance):
                option1 = dp[i + 1][j - 1][0] + 1  # add a pair
                option1_points = dp[i + 1][j - 1][1] + get_pair_score(seq[i], seq[j])  # add points
            else:
                option1 = 0
                option1_points = 0

            # option 2: skip xj - take the value from dp(i + 1, j)
            option2 = dp[i + 1][j][0]
            option2_points = dp[i + 1][j][1]

            # option 3: skip xi - take the value from dp(i, j - 1) 
            option3 = dp[i][j - 1][0]
            option3_points = dp[i][j - 1][1]

            # option 4: split at every possible k and calculate the sum of dp(i, k) and dp(k + 1, j)
            # from research - Bifurcation
            option4 = 0
            option4_points = 0
            for k in range(i, j):
                if dp[i][k][0] + dp[k + 1][j][0] > option4:
                    option4 = dp[i][k][0] + dp[k + 1][j][0]
                    option4_points = dp[i][k][1] + dp[k + 1][j][1]
                elif dp[i][k][0] + dp[k + 1][j][0] == option4:
                    option4_points = max(option4_points, dp[i][k][1] + dp[k + 1][j][1])

            # store the best option in the dp table - prioritize pairs first, then points
            options = [
                [option1, option1_points],
                [option2, option2_points],
                [option3, option3_points],
                [option4, option4_points]
            ]
            
            # select the best option by comparing first pairs, then points
            best_option = max(options, key=lambda x: (x[0], x[1]))  # First by pairs, then by points
            dp[i][j][0] = best_option[0]  # max pairs
            dp[i][j][1] = best_option[1]  # associated points

    return dp

# traceback function to get the optimal structure and score
def traceback(dp, seq, i, j, structure, min_distance=4):
    #base case is if i ≥ j then break
    if i < j:
        # option 1: No pairing at position j
        if dp[i][j] == dp[i][j - 1]:
            traceback(dp, seq, i, j - 1, structure, min_distance)
        # option 2: if i and j are paired, choose this option
        elif is_pair(seq[i], seq[j]) and dp[i][j] == [dp[i + 1][j - 1][0] + 1, dp[i + 1][j - 1][1] + get_pair_score(seq[i], seq[j])] and (j - i >= min_distance):
            structure[i] = "("
            structure[j] = ")"
            traceback(dp, seq, i + 1, j - 1, structure, min_distance)
        else:
            # option 3: try splitting at some k
            for k in range(i, j):
                if dp[i][j] == [dp[i][k][0] + dp[k + 1][j][0], dp[i][k][1] + dp[k + 1][j][1]]:
                    traceback(dp, seq, i, k, structure, min_distance) #traceback 1st part of split
                    traceback(dp, seq, k + 1, j, structure, min_distance) #traceback 2nd part of split
                    break

# main function to run the algorithm and print the structure with scores
def predict_rna_structure(seq, min_distance=4):
    dp = build_table(seq, min_distance) # call build table function

    print("DP Table")
    for row in dp:
        print(row)
    
    # initialize the structure - unpaired bases - "."
    structure = ['.'] * len(seq)  

    # call traceback
    traceback(dp, seq, 0, len(seq) - 1, structure, min_distance)
    final_structure = "".join(structure) # structure into string format

    # max number of pairs is the value at dp[0][len(seq) - 1][0]
    max_pairs = dp[0][len(seq) - 1][0]
    
    # max score are calculated from the dp table at end as well, index 1
    max_score = dp[0][len(seq) - 1][1]
    
    return final_structure, max_pairs, max_score

# input RNA sequence of C, G, U, A
sequence = 'GUGAGAC'

# call to function to run program on input sequence
structure, max_pairs, max_score = predict_rna_structure(sequence)
print(f"RNA Sequence:        {sequence}")
print(f"Optimal Structure:   {structure}")
print(f"Max Number of Pairs: {max_pairs}")
print(f"Max Score:           {max_score}")

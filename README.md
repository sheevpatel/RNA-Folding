# RNA Secondary Structure Predictor

## Overview
This program predicts the optimal secondary structure of an RNA sequence by maximizing base pairings and associated scores using a dynamic programming approach. The algorithm also considers a minimum distance constraint between paired bases.

The prediction includes:
- The RNA secondary structure in dot-parenthesis notation.
- The maximum number of base pairs.
- The maximum score based on the pairings.

## Features
1. **Base Pair Validation:** Checks whether two bases can pair (A-U, C-G, G-U).
2. **Pair Scoring:** Assigns scores to base pairs based on their type:
   - A-U or U-A: 2 points
   - C-G or G-C: 3 points
   - G-U or U-G: 1 point
3. **Dynamic Programming (DP) Algorithm:**
   - Builds a DP table to compute the maximum number of base pairs and their score.
   - 3D table however size is n x n x 2 (not n^3)
   - Considers multiple options for pairing or skipping bases to ensure optimality.
4. **Traceback:** Recovers the optimal RNA structure from the DP table.
5. **User Input:** Allows the user to input any RNA sequence consisting of A, U, C, and G bases.

## Requirements
- Python 3.x

No additional libraries are required.

### Output
The program outputs:
- The dynamic programming table (for debugging purposes).
- The RNA sequence.
- The predicted secondary structure in dot-parenthesis notation.
- The maximum number of base pairs.
- The maximum score.


## Key Functions

### `is_pair(base1, base2)`
Checks if two RNA bases can form a valid pair.

### `get_pair_score(base1, base2)`
Returns the score for a given base pair.

### `build_table(seq, min_distance=4)`
Builds the DP table to find the optimal RNA structure based on the number of pairs and their scores.

### `traceback(dp, seq, i, j, structure, min_distance=4)`
Recovers the RNA structure by backtracking through the DP table.

### `predict_rna_structure(seq, min_distance=4)`
The main function that coordinates the prediction process and outputs the results.

## Customization
- **Minimum Distance Constraint:** Modify the `min_distance` parameter in the function calls to change the minimum allowable distance between paired bases.


## Notes
- This implementation prioritizes maximizing the number of base pairs before maximizing the score.
- The dynamic programming table is printed to facilitate debugging and understanding of the algorithm.


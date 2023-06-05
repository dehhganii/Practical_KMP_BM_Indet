**Overview**
 
The goal of the KMP-INDET algorithm is identifying matches of an indeterminate pattern in an indeterminate test. The KMP-INDET algorithm extends the Knuth-Morris-Pratt algorithm for improved space and time efficiency.

**Pseudocode** 

1. Define necessary headers and libraries.
2. Implement utility functions:
    *isPrime*: Checks if a number is prime.
    *gen_prime*: Generates a vector of prime numbers.
    *fill_check_indet*: Fills the *check_indet* vector with values indicating whether a number is indeterminate or not.
    *gcd*: Computes the greatest common divisor of two numbers.
    *borderarray*: Computes the border array for a given pattern.
    *z_array*: Computes the Z array for a given sequence.
    *match*: Finds the position in the text where the gcd condition fails while comparing elements.
    *prefixarray*: Computes the prefix array for a given text.
3. Implement the *indet* function to check if a number is indeterminate.
4. Implement the *compute_shift* function to calculate the shift value based on the given conditions.
5. Implement the *KMP_Indet* function:
    Initialize variables and containers.
    Compute the longest prefix m_ell in the pattern that is not indeterminate.
    Compute the border array for the pattern.
    Iterate over the text:
        Compare elements in the text and pattern based on gcd conditions.
        If a match is found:
            Check for complete pattern match.
            Compute the shift value based on indeterminate conditions.
        If a mismatch occurs:
            Compute the shift value based on indeterminate conditions.
    Return the number of occurrences found.
6. Implement the *read_file* function to read comma-separated integers from a file and return them as a vector.
7. Implement the main function:
    Parse command-line arguments.
    Read the text and pattern from files.
    Fill the *check_indet* vector.
    Measure the execution time of the KMP algorithm.
    Print the execution time.

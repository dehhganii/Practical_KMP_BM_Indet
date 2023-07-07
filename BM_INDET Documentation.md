**Overview**
 
The goal of the BM-INDET algorithm is identifying matches of an indeterminate pattern in an indeterminate test. The BM-INDET algorithm extends the Boyer-Moore algorithm for improved space and time efficiency.

**Pseudocode**

1. Define necessary headers and libraries.
2. Implement utility functions:
    *gcd*: Computes the greatest common divisor of two numbers.
    *copy*:Copies elements from one array to another based on certain conditions involving pos,myequal,and x.
    *mymatch*: Compares elements at two positions and updates array accordingly.
    *poslen*: Computes the length of a pattern in a given sequence and updating relevent arrays.
    *isPrime*: Checks if a number is prime.
    *gen_prime*: Generates a vector of prime numbers.
    *fill_check_indet*: Fills the *check_indet* vector with values indicating whether a number is indeterminate or not.
    *costum_map*: Creates custom mapping between prime numbers and their corresponding indices.
    *roots*:Finds the prime factors of a number.
    *z_array*: Computes the Z array for a given sequence.
    *n_array*:Computes the N array for a given sequence.
    *big_l_prime_array*:Calculates big l prime array for a given sequence.
    *big_l_array*:Calculates big l array for a given sequence.
    *small_l_prime_array*:Calculates small l prime array for a given pattern sequence.
    *create_dense_bad_char_table*:Create a dense bad character table object based on the given sequence.
    *indet_good_suffix*:Handles mismatches ocurring during string matching. 
3. Implement the *dense_bad_char_shift* and *good_suffix_rule_fun* to calculate shift distances based on certain conditions. 
4. Implement the *boyer_more_indet* function:
    Initialize variables and containers.
    Compute *n, big_l_prime, big_l, small_l_prime, and amap* arrays
    Compute *m_ell_suf*, the length of the longest suffix in p that contains indeterminate characters
    Calculate value of period using *len_p* and *small_l_prime*
    Iterate over the text,backwards:
        Check if element is indeterminate. 
        Compare elements in the text and pattern based on gcd conditions.
            if gcd=1:
                Computes bad character shift.
                Checks if good suffix shift determined based on indeterminate charcters or good suffix rule
                Max value between two methods assigned to shift
                Mismatch is set to true
        if mismatch !=True:
            Increment number of occurrences 
            Check if  indeterminate characters need to be updated based on  current position.
            Determine skip shift for the good suffix based on indeterminate characters or period.
            Update shift with the maximum value of shift and the skip shift.
    Free allocated memory
    Return the number of occurrences found.

5. Implement the *read_file* function to read comma-separated integers from a file and return them as a vector.
6. Implement the main function:
    Parse command-line arguments.
    Read the text and pattern from files.
    Fill the *check_indet* vector.
    Measure the execution time of the BM algorithm.
    Print the execution time.
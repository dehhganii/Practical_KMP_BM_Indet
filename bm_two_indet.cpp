#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <array>
#include <vector>
#include <cassert>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

/**
 * @brief Return minimum of two numbers
 * 
 */
#define MIN(a,b) ((a) < (b) ? (a) : (b))
/**
 * @brief Return maximum of two numbers
 * 
 */ 
#define MAX(a,b) ((a) > (b) ? (a) : (b))
/**
 * @brief Determine if number is indeterminate or not
 * 
 */
#define INDET(n) (!check_indet[(n)])


int bit_AND(int a, int b) {
    return a & b;
}
/**
 * @brief Copies elements from one array(len) to another based on certain condition involving pos,myequal,and x
 * 
 * @param J Index value used in loop 
 * @param lim Limit value to determine if loop should keep iterating
 * @param pos Pointer to an array that holds integer values. Used in comparisons and computations within the loop
 * @param len Pointer to an arrary used to modify elements of the array based on certain conditions
 * @param x Pointer to an array used to access elements for computation within the loop
 * @param m Integer value
 * @param myequal Pointer to integer array used to determine a condition within the loop
 * @param n Integere value
 */
void copy(int J, int lim, int *pos, int *len, int *x, int m, int *myequal, int n) {
  int Jprime, I;
  int pos1, pos2;

  Jprime = 0;
  while (J > n && pos[J] > lim) {
    I = pos[J];
    if (myequal[Jprime]) {
      len[J] = MIN(len[Jprime], I-lim+1);
    }
    else {
      pos1 = m-2;
      pos2 = I-1;
      while (pos2 >= 0 && bit_AND(x[pos1], x[pos2]) != 0) {
        --pos1;
        --pos2;
      }
      len[J] = m-1-pos1;
    }
    ++J;
    ++Jprime;
  }
}

/**
 * @brief Compare elements at pos1 and pos2 and update myequal array accordingly
 * 
 * @param pos1 Index value in x array
 * @param pos2 Another index value in x array
 * @param j Index value in myequal array
 * @param x Pointer to x array
 * @param myequal Pointer to myequal array
 * @return pos1 Final value of pos1
 */
int mymatch(int pos1, int pos2, int j, int *x, int *myequal) {
  while (pos2 >= 0 && bit_AND(x[pos1], x[pos2]) != 0) {
    if (x[pos1] != x[pos2]) {
      myequal[j] = 0;
    }
    --pos1;
    --pos2;
  }
  return(pos1);
}

/**
 * @brief Computes the length of a pattern in a given sequence and updating relevent arrays
 * 
 * @param x Pointer to x array
 * @param len_p Length of pattern
 * @param matched_t Length of matched substring
 * @param pos Pointer to array that stores positions within x array, based on certain conditions
 * @param len Pointer to array that stores lengths
 * @param myequal Pointer to myequal array
 * @return Length of the pattern or the difference between len_p and position in the pos array
 */
int poslen(int *x, int len_p, int matched_t, int *pos, int *len, int *myequal) {
  int lambda, n, j, lim, suff;
  int i, end;
  int m = len_p+matched_t+1;

  lambda = x[m-1];
  n = 0;
  for (i = len_p-2; i >= 0; --i) {
    if (bit_AND(x[i], lambda) != 0) {
      pos[n] = i;
      len[n] = 1;
      myequal[n] = x[i] == lambda;
      ++n;
    }
  }
  j = 0;
  lim = len_p-1;
  while (j < n) {
    i = pos[j];
    suff = len[j];
    if (i-suff < lim) {
      suff = m-1-mymatch(m-1-suff, i-suff, j, x, myequal);
      len[j] = suff;
      //if (pos[j] < len_p-1 && len[j] == matched_t) {
      if (len[j] == matched_t) {
        return(len_p-1-pos[j]);
      }
      if (pos[j] < matched_t && pos[j] == len[j]-1) {
        return(len_p-1-pos[j]);
      }
      end = i-suff+1;
      if (end < lim) {
        lim = end;
        copy(j+1, lim, pos, len, x, m, myequal, n);
      }
    }
    ++j;
  }
  return(len_p);
}

/**
 * @brief Initilizes variables used for measuring execution time
 * 
 */
clock_t start_BM, finish_BM;

/**
 * @brief Global vector declared with a size of 30 elements and each element is initialized to a value of -1
 * 
 */
vector<int> check_indet;

void init_indet(int sigma){
    int max = pow(2,sigma);
    check_indet.resize(max, 0);
}

/**
 * @brief Checks if number is power
 * 
 * @param n Number that is being checked
 * @return true If n is prime 
 * @return false If n is not prime
 */
bool isPowerOfTwo(int n) {
    if (n == 0)
        return false;
    while (n % 2 == 0)
        n /= 2;
    return n == 1;
}

/**
 * @brief Generates an n number of prime numbers and stores them in the vector prime
 * 
 * @param n Number of prime numbers generated
 * @return vector<int> primes Vector with all generated prime numbers
 */
vector<int> gen_twos(int n) {
    vector<int> twos;
    int i = 1;
    while (twos.size() < n) {
        if (isPowerOfTwo(i)) {
            twos.push_back(i);
        }
        i++;
    }
    return twos;
}

/**
 * @brief Fills check_indet vector with a value of 1 corresponding to the indices of the values found in primes
 * 
 * @param primes Vector with prime number values 
 */
void fill_check_indet(const vector<int>& twos) {
    for (int i = 0; i < twos.size(); i++) {
        check_indet[twos[i]] = 1;
    }
}


/**
 * @brief Creates custom mapping between prime numbers and their corresponding indices
 * 
 * @param sigma Size of alphabet
 * @param prime Vector madeup of prime integers
 * @return amap Cutsom mapping of prime numbers to respective indices
 */
vector<int> costum_map(int sigma, const vector<int> &two)
{
    vector<int> amap(two[sigma - 1] + 1, -1); 
    for (int i = 0; i < sigma; i++)
    {
        amap[two[i]] = i;
    }
    return amap;
}

/**
 * @brief Finds the prime factors(roots) of a number
 * 
 * @param n Integer
 * @param sigma Size of alphabet
 * @param prime Vector madeup of prime integers
 * @return root_vector All prime factors of n 
 */
vector<int> roots(int n, int sigma, const vector<int> &two)
{
    vector<int> root_vector;
    for (int i = 0; i < sigma; i++)
    {
        if (n % two[i] == 0)
        {
            root_vector.push_back(two[i]);
        }
    }
    return root_vector;
}

/**
 * @brief  Z-array of a sequence
 * 
 * @param s Sequence
 * @param len Length of the sequence
 * @return z Z-array  
 */
vector<int> z_array(const vector<int> &s, int len)
{
    assert(len > 1);
    vector<int> z(len);
    z[0] = len;
    for (int i = 1; i < len; i++)
    {
        if (bit_AND(s[i], s[i - 1]) > 0)
        {
            z[1] += 1;
        }
        else
        {
            break;
        }
    }
    int r = 0;
    int l = 0;
    if (z[1] > 0)
    {
        r = z[1];
        l = 1;
    }
    for (int k = 2; k < len; k++)
    {
        assert(z[k] == 0);
        if (k > r)
        {
            for (int j = k; j < len; j++)
            {
                if (bit_AND(s[j], s[j - k]) > 0)
                {
                    z[k] += 1;
                }
                else
                {
                    break;
                }
            }
            r = k + z[k] - 1;
            l = k;
        }
        else
        {
            int nbeta = r - k + 1;
            int zkp = z[k - 1];
            if (nbeta > zkp)
            {
                z[k] = zkp;
            }
            else
            {
                int nmatch = 0;
                for (int i = r + 1; i < len; i++)
                {
                    if (bit_AND(s[i], s[i - k]) > 0)
                    {
                        nmatch += 1;
                    }
                    else
                    {
                        break;
                    }
                }
                l = k;
                r = r + nmatch;
                z[k] = r - k + 1;
            }
        }
    }
    return z;
}

/**
 * @brief Calculates n-array of a given sequence
 * 
 * @param s Input sequence
 * @param len Length of s
 * @return z_tmp_rev N-array 
 */
vector<int> n_array(const vector<int> &s, int len)
{
    vector<int> tmp_s(len);
    for (int i = 0; i < len; i++)
    {
        tmp_s[i] = s[len - 1 - i];
    }
    vector<int> z_tmp = z_array(tmp_s, len);
    vector<int> z_tmp_rev(len);
    for (int i = 0; i < len; i++)
    {
        z_tmp_rev[i] = z_tmp[len - 1 - i];
    }
    return z_tmp_rev;
}

/**
 * @brief Calculates big l prime array for a given sequence
 * 
 * @param p Input sequence
 * @param len Length of p 
 * @return big_l_prime_array_s Big L prime array of p 
 */
vector<int> big_l_two_array(const vector<int> &p, int len)
{
    vector<int> big_l_two_array_s(len);
    vector<int> n = n_array(p, len);
    int i = 0;
    for (int j = 0; j < len - 1; j++)
    {
        i = len - n[j];
        if (i < len)
        {
            big_l_two_array_s[i] = j + 1;
        }
    }
    return big_l_two_array_s;
}

/**
 * @brief Calculated big l array for a given sequence
 * 
 * @param pattern Pattern sequence
 * @param big_l_prime Big l prime array of the pattern sequence
 * @param len Length of pattern sequence
 * @return big_l_array_s Big l array
 */
vector<int> big_l_array(const vector<int> &pattern, const vector<int> &big_l_two, int len)
{
    vector<int> big_l_array_s(len);
    big_l_array_s[1] = big_l_two[1];
    for (int j = 2; j < len; j++)
    {
        big_l_array_s[j] = MAX(big_l_array_s[j - 1], big_l_two[j]);
    }
    return big_l_array_s;
}
/**
 * @brief Calculates small l prime array for a given pattern sequence
 * 
 * @param n N-array
 * @param len Length of pattern sequence
 * @return small_l_prime_array_s Small l prime array 
 */
vector<int> small_l_two_array(vector<int> n, int len)
{
    vector<int> small_l_two_array_s(len);
    for (int j = 0; j < len; j++)
    {
        if (n[j] == j + 1)
        {
            small_l_two_array_s[len - j - 1] = j + 1;
        }
    }
    for (int j = len - 2; j > -1; j--)
    {
        if (small_l_two_array_s[j] == 0)
        {
            small_l_two_array_s[j] = small_l_two_array_s[j + 1];
        }
    }
    return small_l_two_array_s;
}


/**
 * @brief Create a dense bad character table object based on the given sequence
 * 
 * @param pattern Pattern sequence
 * @param pat_len Length of pattern sequence 
 * @param sigma Size of alphabet
 * @param amap Mapping of characters from corresponding indices
 * @param primes Vector of prime integers
 * @return tab Dense bad character shift table 
 */
vector<vector<int> > create_dense_bad_char_table(const vector<int> &pattern, int pat_len, int sigma, const vector<int> &amap, const vector<int> &twos) {
    vector<vector<int> > tab;
    vector<int> nxt(sigma, 0);
    for(int i=0; i < pat_len; i++) {
        int c = pattern[i];
        if(INDET(c)) {
            vector<int> root_vector = roots(c, sigma,twos);
            for(int j=0; j < root_vector.size(); j++) {
                nxt[amap[root_vector[j]]] = i + 1;
            }
        }
        else {
            nxt[amap[c]] = i + 1;
        }
        tab.push_back(nxt);
    }
    return tab;
}

/**
 * @brief Calculates shift distance for given index and character
 * 
 * @param i Current index
 * @param c Character 
 * @param bad_char_table Bad character shift table
 * @param amap Custom mapping of characters to respective indices
 * @param sigma Size of alphabet
 * @param primes Vector of prime numbers
 * @return Shift distance
 */
int dense_bad_char_shift(int i, int c, const vector<vector<int> > &bad_char_table, const vector<int> &amap, int sigma, vector<int> twos) {
    if(!INDET(c)) {
        return i - (bad_char_table[i][amap[c]]-1);
    }
    else {
        
        vector<int> c_roots = roots(c, sigma, twos);
        int shift = 0;
        int bound = c_roots.size();
        for(int k = 0; k < bound; k++) {
            shift = MAX(shift, bad_char_table[i][amap[c_roots[k]]]-1);
        }
        
        return i - shift;
    }
}

/**
 * @brief Calculates shift distance based on current index and values in big l array and small prime array
 * 
 * @param i Current index
 * @param big_l_array_s Big l array 
 * @param small_l_prime_array_s Small l prime array
 * @param len_big_l_array Length of big l array
 * @return Shift distance
 */
int good_suffix_rule_fun(int i,const vector<int> &big_l_array_s, const vector<int> &small_l_two_array_s, int len_big_l_array){
    if(i >= len_big_l_array){
        return -1;
    }
    if(i == len_big_l_array - 1){
        return 0;
    }
    int j = i + 1;
    if(big_l_array_s[j] > 0){
        return len_big_l_array - big_l_array_s[j];
    }
    return len_big_l_array - small_l_two_array_s[j];

}

/**
 * @brief Handles mismatches ocurring during string matching
 * 
 * @param pattern Pattern sequence
 * @param len_p Length of pattern
 * @param text Text sequence
 * @param i Current index
 * @param j Current index within pattern
 * @param stringForSuff Pointer to modified pattern sequence array
 * @param pos Pointer to array storing position information
 * @param len Pointer to array stroing length of substring
 * @param myequal Pointer to integer storing a flag 
 * @return poslen Length of the pattern or the difference between len_p and position in the pos array
 */
int indet_good_suffix(const vector<int> &pattern, int len_p, const vector<int> &text, int i, int j,
   int *stringForSuff, int *pos, int *len, int *myequal) {   
    int matched_t = 0;
    if( j==0) {
        matched_t = len_p;
    }
    if(j == len_p) 
    {
        return 1;
    }
    else {
        //j holds the index position of the mismatched letter in the pattern
        matched_t = len_p - j;

    }
    
    for (int tmp = 0; tmp < matched_t ; ++tmp ) {
      // use memcpy?
      stringForSuff[len_p+1+tmp] = text[i+j+1+tmp];
    }

    return poslen(stringForSuff, len_p, matched_t, pos, len, myequal);

}


/**
 * @brief Implements BM algorithm on indeterminate strings
 * 
 * @param t Text sequence
 * @param p Pattern sequence
 * @param len_t Length of text
 * @param len_p Length of pattern
 * @param sigma Size of alphabet
 * @param primes Vector of prime numbers
 * @return occurrences_len Number of pattern occurrences found in the text 
 */
int boyer_more_indet(const vector<int> &t, const vector<int> &p, int len_t, int len_p, int sigma, const vector<int> &twos)
{
    vector<int> n = n_array(p, len_p);
    vector<int> big_l_two = big_l_two_array(p, len_p);
    vector<int> big_l = big_l_array(p, big_l_two, len_p);
    vector<int> small_l_two = small_l_two_array(n, len_p);
    vector<int> amap = costum_map(sigma, twos);
    vector<vector<int> > bad_char_table = create_dense_bad_char_table(p, len_p, sigma, amap, twos);
    int *stringForSuff;
    int *pos;
    int *len;
    int *myequal;
    
    int i;
    int shift;
    int m_ell_suf = 0;
    bool mismatched;
    int right_most_indet = -1;
    int l;
    int skip_gs;
    int gs_shift;

    l = (2*len_p+1)*sizeof(int);
    stringForSuff = (int *)malloc(l);
    if (stringForSuff == NULL) exit(1);
    pos = (int *)malloc(l);
    if (pos == NULL) exit(1);
    len = (int *)malloc(l);
    if (len == NULL) exit(1);
    myequal = (int *)malloc(l);
    if (myequal == NULL) exit(1);

    for (i = 0; i < len_p; ++i) {
      stringForSuff[i] = p[i];
    }
    stringForSuff[len_p] = 1;
    
    l = len_p-1;
    //TODO: the m_ell value should be longest 
    while (l >= 0)
    {
        if (INDET(p[l]))
        {
            m_ell_suf = len_p-1-l;
            break;
        }
        --l;
    }
    
    if (l == -1) { m_ell_suf = len_p; }
    int occurrences_len = 0;
    i = 0;
    int bound = len_t - len_p + 1;
    int len_pMinus1 = len_p-1;
    int period = len_p - small_l_two[1];
    while (i < bound)
    {
        shift = 1;
        mismatched = false;

        for (int j = len_pMinus1; j >= 0; --j)
        {
            if (INDET(t[i + j]))
            {
                right_most_indet = i+j;
            }
            if (bit_AND(p[j], t[i + j]) == 0)
            {
                int bc_shift = dense_bad_char_shift(j, t[i + j], bad_char_table, amap, sigma, twos);

                if ((i + j >= right_most_indet && right_most_indet <= i +len_pMinus1) || len_pMinus1 - j >= m_ell_suf)
                {
                    gs_shift = indet_good_suffix(p, len_pMinus1, t, i, j, stringForSuff, pos, len, myequal);
                }
                else
                {
                    gs_shift = good_suffix_rule_fun(j, big_l, small_l_two, len_p);

                }
                shift = MAX(shift, bc_shift);
                shift = MAX(shift, gs_shift);
                mismatched = true;
                break;
            }
        }
        if (!mismatched)
        {
            occurrences_len++;
            //cout << " the pattern found in this place of text: " << i << ", " << endl;
            //TODO: add code to update indet
            if ((i >= right_most_indet && right_most_indet <= i +len_pMinus1) || m_ell_suf < len_p)
            {
                // printf("wrong indet");
                skip_gs = indet_good_suffix(p, len_pMinus1, t, i, 0, stringForSuff, pos, len, myequal);
            }
            else
            {
                skip_gs = period;
            }
            shift = MAX(shift, skip_gs);
        }

        i += shift;
    }
    free(stringForSuff);
    free(pos);
    free(len);
    free(myequal);
    return occurrences_len;
}


/**
 * @brief Reads file and extracts necessary text 
 * 
 * @param filename The filename for the text being searched 
 * @return result Vector of the integers read from the file 
 */
vector<int> read_file(string filename) {
    vector<int> result;
    ifstream file(filename);
    if (file) {
        string line;
        while (getline(file, line)) {
            stringstream ss(line);
            string number;
            while (getline(ss, number, ',')) {
                result.push_back(stoi(number));
            }
        }
    } else {
        cerr << "Error: Could not open file '" << filename << "'" << endl;
    }
    return result;
}

    int main(int argc, char *argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: ./boyer_moore_out <alphabet_size> <text file> <pattern file>\n";
        return -1;
    }
    
    int sigma = atoi(argv[1]);
    init_indet(sigma);
    vector<int> twos = gen_twos(sigma);
    fill_check_indet(twos);
    string text_file_name = argv[2];
    vector<int> test_text = read_file(text_file_name);
    string pattern_file_name = argv[3];
    vector<int> test_pattern = read_file(pattern_file_name);
    
    long int text_length = test_text.size();
    long int pattern_length = test_pattern.size();
    start_BM = clock();
    vector<int> amap = costum_map(sigma, twos);
    int answers = boyer_more_indet(test_text, test_pattern, text_length, pattern_length, sigma, twos);
    //vector<vector<int>> tab = create_dense_bad_char_table(test_pattern, pattern_length, sigma, amap, primes);
    //dense_bad_char_shift(3, 21, tab, amap, sigma, primes);
    finish_BM = clock();
    double BM_time = (double)((finish_BM - start_BM) / (double)CLOCKS_PER_SEC);
    //cout << "number_of_matches" << answers << endl;
    printf("%.7f", BM_time);
    
    cout << "number of occurences" << answers << endl;
    return 0;
}

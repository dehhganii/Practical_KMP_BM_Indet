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

using namespace std;


#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define INDET(n) ((n) >= 30 || !check_indet[(n)])

int gcd(int a, int b) {
  if (a == 0) return b;
  return gcd(b % a, a);
}

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
      while (pos2 >= 0 && gcd(x[pos1], x[pos2]) != 1) {
        --pos1;
        --pos2;
      }
      len[J] = m-1-pos1;
    }
    ++J;
    ++Jprime;
  }
}


int mymatch(int pos1, int pos2, int j, int *x, int *myequal) {
  while (pos2 >= 0 && gcd(x[pos1], x[pos2]) != 1) {
    if (x[pos1] != x[pos2]) {
      myequal[j] = 0;
    }
    --pos1;
    --pos2;
  }
  return(pos1);
}


int poslen(int *x, int len_p, int matched_t, int *pos, int *len, int *myequal) {
  int lambda, n, j, lim, suff;
  int i, end;
  int m = len_p+matched_t+1;

  lambda = x[m-1];
  n = 0;
  for (i = len_p-2; i >= 0; --i) {
    if (gcd(x[i], lambda) != 1) {
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


clock_t start_BM, finish_BM;
vector<int> check_indet(30, 0);
bool isPrime(int n)
{
    if (n <= 1)
    {
        return false;
    }
        for (int i = 2; i <= n / 2; i++)
    {
        if (n % i == 0)
        {
            return false;
        }
    }

    return true;
}

vector<int> gen_prime(int n)
{
    vector<int> primes;
    int i = 2;

    while (primes.size() < n)
    {
        if (isPrime(i))
        {
            primes.push_back(i);
        }
        i++;
    }

    return primes;
}

void fill_check_indet(const vector<int> &primes)
{
    for(int i=0; i < primes.size(); i++)
    {
        check_indet[primes[i]] = 1;
    }

}



vector<int> costum_map(int sigma, const vector<int> &prime)
{
    vector<int> amap(prime[sigma - 1] + 1, -1); 
    for (int i = 0; i < sigma; i++)
    {
        amap[prime[i]] = i;
    }
    return amap;
}

vector<int> roots(int n, int sigma, const vector<int> &prime)
{
    vector<int> root_vector;
    for (int i = 0; i < sigma; i++)
    {
        if (n % prime[i] == 0)
        {
            root_vector.push_back(prime[i]);
        }
    }
    return root_vector;
}


vector<int> z_array(const vector<int> &s, int len)
{
    assert(len > 1);
    vector<int> z(len);
    z[0] = len;
    for (int i = 1; i < len; i++)
    {
        if (gcd(s[i], s[i - 1]) > 1)
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
                if (gcd(s[j], s[j - k]) > 1)
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
                    if (gcd(s[i], s[i - k]) > 1)
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

vector<int> big_l_prime_array(const vector<int> &p, int len)
{
    vector<int> big_l_prime_array_s(len);
    vector<int> n = n_array(p, len);
    int i = 0;
    for (int j = 0; j < len - 1; j++)
    {
        i = len - n[j];
        if (i < len)
        {
            big_l_prime_array_s[i] = j + 1;
        }
    }
    return big_l_prime_array_s;
}

vector<int> big_l_array(const vector<int> &pattern, const vector<int> &big_l_prime, int len)
{
    vector<int> big_l_array_s(len);
    big_l_array_s[1] = big_l_prime[1];
    for (int j = 2; j < len; j++)
    {
        big_l_array_s[j] = MAX(big_l_array_s[j - 1], big_l_prime[j]);
    }
    return big_l_array_s;
}

vector<int> small_l_prime_array(vector<int> n, int len)
{
    vector<int> small_l_prime_array_s(len);
    for (int j = 0; j < len; j++)
    {
        if (n[j] == j + 1)
        {
            small_l_prime_array_s[len - j - 1] = j + 1;
        }
    }
    for (int j = len - 2; j > -1; j--)
    {
        if (small_l_prime_array_s[j] == 0)
        {
            small_l_prime_array_s[j] = small_l_prime_array_s[j + 1];
        }
    }
    return small_l_prime_array_s;
}



vector<vector<int> > create_dense_bad_char_table(const vector<int> &pattern, int pat_len, int sigma, const vector<int> &amap, const vector<int> &primes) {
    vector<vector<int> > tab;
    vector<int> nxt(sigma, 0);
    for(int i=0; i < pat_len; i++) {
        int c = pattern[i];
        if(INDET(c)) {
            vector<int> root_vector = roots(c, sigma,primes);
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


int dense_bad_char_shift(int i, int c, const vector<vector<int> > &bad_char_table, const vector<int> &amap, int sigma, vector<int> primes) {
    if(!INDET(c)) {
        return i - (bad_char_table[i][amap[c]]-1);
    }
    else {
        
        vector<int> c_roots = roots(c, sigma, primes);
        int shift = 0;
        int bound = c_roots.size();
        for(int k = 0; k < bound; k++) {
            shift = MAX(shift, bad_char_table[i][amap[c_roots[k]]]-1);
        }
        
        return i - shift;
    }
}


int good_suffix_rule_fun(int i,const vector<int> &big_l_array_s, const vector<int> &small_l_prime_array_s, int len_big_l_array){
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
    return len_big_l_array - small_l_prime_array_s[j];

}


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



int boyer_more_indet(const vector<int> &t, const vector<int> &p, int len_t, int len_p, int sigma, const vector<int> &primes)
{
    vector<int> n = n_array(p, len_p);
    vector<int> big_l_prime = big_l_prime_array(p, len_p);
    vector<int> big_l = big_l_array(p, big_l_prime, len_p);
    vector<int> small_l_prime = small_l_prime_array(n, len_p);
    vector<int> amap = costum_map(sigma, primes);
    vector<vector<int> > bad_char_table = create_dense_bad_char_table(p, len_p, sigma, amap, primes);
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
    int period = len_p - small_l_prime[1];
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
            if (gcd(p[j], t[i + j]) == 1)
            {
                int bc_shift = dense_bad_char_shift(j, t[i + j], bad_char_table, amap, sigma, primes);

                if ((i + j >= right_most_indet && right_most_indet <= i +len_pMinus1) || len_pMinus1 - j >= m_ell_suf)
                {
                    gs_shift = indet_good_suffix(p, len_pMinus1, t, i, j, stringForSuff, pos, len, myequal);
                }
                else
                {
                    gs_shift = good_suffix_rule_fun(j, big_l, small_l_prime, len_p);

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
    vector<int> primes = gen_prime(sigma);
    fill_check_indet(primes);
    string text_file_name = argv[2];
    vector<int> test_text = read_file(text_file_name);
    string pattern_file_name = argv[3];
    vector<int> test_pattern = read_file(pattern_file_name);
    
    long int text_length = test_text.size();
    long int pattern_length = test_pattern.size();
    start_BM = clock();
    vector<int> amap = costum_map(sigma, primes);
    int answers = boyer_more_indet(test_text, test_pattern, text_length, pattern_length, sigma, primes);
    //vector<vector<int>> tab = create_dense_bad_char_table(test_pattern, pattern_length, sigma, amap, primes);
    //dense_bad_char_shift(3, 21, tab, amap, sigma, primes);
    finish_BM = clock();
    double BM_time = (double)((finish_BM - start_BM) / (double)CLOCKS_PER_SEC);
    //cout << "number_of_matches" << answers << endl;
    printf("%.7f", BM_time);
    
    //cout << "number of occurences" << answers << endl;
    return 0;
}

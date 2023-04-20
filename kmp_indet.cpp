#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

using namespace std;

clock_t start_KMP, finish_KMP;


#define MAX(a,b) (((a)>(b))?(a):(b))

vector<int> check_indet(30, -1);
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

vector<int> borderarray(const vector<int>& pat, int n) {
    vector<int> border_array(n);
    border_array[0] = 0;
    int j = 0;
    for (int i = 1; i < n; i++) {
        while (j > 0 && pat[i] != pat[j]) {
            j = border_array[j - 1];
        }
        if (pat[i] == pat[j]) {
            j++;
        }
        border_array[i] = j;
    }
    return border_array;
}

int gcd(int a, int b)
{
    if (a == 0)
        return b;
    return gcd(b%a, a);
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

int match(const vector<int> &text, int pos1, int pos2, int j, int n){
    while(gcd(text[pos1], text[pos2])> 1 && pos1 < n && pos2 < n){
        pos1 = pos1 + 1;
        pos2 = pos2 + 1;
    }
    return pos1;
}

vector<int> prefixarray(vector<int> text, int n){
    int lambda = text[0];
    text.push_back(-1);
    int i, j, lim, pref; //end;
    vector<int> prefix(n, 0);
    for(i=0;i<n;i++){
        if(gcd(text[i],lambda) > 1){
            prefix[i] = 1;
        }
    }
    prefix[0]=n;
    j = 1; lim = 1;
    while(j < n){
        if(prefix[j] != 0){
            pref = prefix[j]; 
            pref = match(text, pref, j+pref, j, n);
            prefix[j]=pref;
        }
        j = j +1;
    }
    
    return prefix;
}

bool indet(int n) {
    if(n > 30)
    {
        return true;
    }
    if(check_indet[n] == 1)
    {
        return false;
    }
    return true;
}


int compute_shift(bool indettext, const vector<int> &text, const vector<int> &pattern, int i, int j, vector<int> border, int m_ell){
    int mx = -1;
    vector<int> q_prime;
    int counter_1 = 0;
    int counter_2 = 0;
    if(indettext || j + 1 > m_ell)
    {
        for(int tmp = 0; tmp < j; tmp++) {
            counter_1 ++;
            q_prime.push_back(pattern[tmp]);
        }
        for(int tmp = i - j + 1; tmp <= i; tmp ++) {
            counter_2++;
            q_prime.push_back(text[tmp]);
        }
    
        if(q_prime.size() <= 0){
            return -1;
        }
    
        vector<int> prefix = prefixarray(q_prime, q_prime.size());
        int half = int(q_prime.size()/2);
        for(int k = half; k < q_prime.size() ;k++){
            if( prefix[k] == q_prime.size() - k && mx < prefix[k]){
                mx=prefix[k] - 1;
            }
        }
        j = mx;
    }
    else{
        j = border[j] - 1;
    }
    return j;
}

int KMP_Indet(int n, int m, int sigma, const vector<int> &text, const vector<int> &pattern){

    vector<int> prime = gen_prime(sigma);
    int i, j;
    int m_ell; 
    int occurences_len = 0;
    bool indettext = false;

    m_ell=m; 
    i = 0;
    while(i<m){
        if(indet(pattern[i])){


            m_ell=i;
            break;
        }
        i=i+1;
    }
    vector<int> border;
    if (m_ell <= 1) {
        border = {0};
    }
    else {
        border = borderarray(pattern, m_ell); //Border array of the longest regular prefix of pattern.
    }
    i=-1;
    j=-1;
    //i and j are the index positions in the text and pattern that match.
    //As a result the substrings text[i-j..i] and pattern[0..j] match.
    int right_most_indet = -1;
    while(i < n - 1){
        if(gcd(text[i+ 1], pattern[j + 1]) > 1)
        {
            if(indet(text[i+1])){
                indettext=true;
                right_most_indet = i + 1;
            }
            j = j+1; //j is the index positions the pattern such that the prefix of length j+1 has match with a substsring of text.
            i = i+1;
            if(j == m-1 ){
                //cout << "the pattern found in this position of text: " <<  i-j << endl;
                occurences_len++; 
                j = compute_shift(indettext, text, pattern, i, j, border, m_ell);
                if( i - j <= right_most_indet) {
                    indettext = true;
                }
                else {
                    indettext = false;
                }
            }
        }
        else
        {
            if(j <= 0)
            {
                if( j < 0)
                    i = i+1;
                j = -1;
                
            }
            else{
                if (j == 1)
                {
                    if(gcd(pattern[0],text[i]) > 1) {
                        j = 0;
                    }
                    else 
                    {
                        j = -1;
                    }
                }
            
                else
                {
                
                    j = compute_shift(indettext, text, pattern, i, j,  border, m_ell);
                }
                
            }
            if( i - j <= right_most_indet) 
            {
                indettext = true;
            }
            else 
                {
                indettext = false;
            }
        }
    }
    return occurences_len;
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
        std::cerr << "Usage: ./first_test_kmp <alphabet_size> <text file> <pattern file>\n";
        return -1;
    }
    
    int sigma = atoi(argv[1]);
    vector<int> primes = gen_prime(sigma);
    string text_file_name = argv[2];
    vector<int> test_text = read_file(text_file_name);
    string pattern_file_name = argv[3];
    vector<int> test_pattern = read_file(pattern_file_name);
    fill_check_indet(primes);
    long int text_length = test_text.size();
    long int pattern_length = test_pattern.size();
    start_KMP = clock();
    
    int answer = KMP_Indet(text_length, pattern_length, sigma, test_text, test_pattern);
    finish_KMP = clock();
    double KMP_TIME = (double)((finish_KMP - start_KMP) / (double)CLOCKS_PER_SEC);
    printf("%.7f", KMP_TIME);
    //cout << "number of occurences" << answer << endl;

    return 0;
}




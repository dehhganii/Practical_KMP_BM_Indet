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

using namespace std;
#include <time.h>

clock_t  start_BF, finish_BF;

vector<int> gen_prime(int sigma){
    vector<int> prime(sigma, 0);
    int n = sigma;
    int temp=1;
    int i=3, count, c;
    if (n>=1){
        prime[0] = 2;
    }
    for(count=2; count<=n; i++){
        for(c=2; c<i; c++){
            if(i%c== 0)
                break;
        }

        if(c==i){
            prime[temp] = i;
            count++;
            temp++;
        }
    }
    return prime;
}

int gcd(int a, int b)
{
    if (a == 0)
        return b;
    return gcd(b%a, a);
}


int bruteforce(vector<int> text, vector<int> pattern, int n, int m){
    int number_of_occurences = 0;
    //int x;
    //int y;//variables for gcdExtended function.
    //When the pattern is algined at 'index' position text[i-j..i] match with pattern[0..j]
    int i, j;
    //is the position in the text for which we are checking wheather the pattern matches or not when aligned at this position.
    int index = 0;
    j=-1;
    i= index-1;
    while(index <= n-m){
        if(gcd(text[i+1], pattern[j+1]) > 1){
            j = j+1;
            i = i+1;
            if(j==m-1){
                number_of_occurences ++;
                //cout << "found the pattern in this index: " << index;
                index = index+1;
                i = index-1;
                j = -1;
            }
        }
        else{
            if(j==-1){
                index = index+1;
                i = index-1;
            }
            else{
                index = index+1;
                i = index-1;
                //i = i-j;
                j = -1;
            }
        }
    }
    return number_of_occurences;
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
        std::cerr << "Usage: ./brute_force_out <alphabet_size> <text file> <pattern file>\n";
        return -1;
    }
    
    int sigma = atoi(argv[1]);
    vector<int> primes = gen_prime(sigma);
    string text_file_name = argv[2];
    vector<int> test_text = read_file(text_file_name);
    string pattern_file_name = argv[3];
    vector<int> test_pattern = read_file(pattern_file_name);
    long int text_length = test_text.size();
    long int pattern_length = test_pattern.size();
    start_BF = clock();
    int answers = bruteforce(test_text, test_pattern, text_length, pattern_length);
    finish_BF = clock();

    double BF_time = (double)((finish_BF - start_BF) / (double)CLOCKS_PER_SEC);
    //cout << "number_of_matches" << answers << endl;
    printf("%.7f", BF_time);
}

import random
from itertools import chain, combinations
import sys


def gen_twos(alphabet_size):
    twos = []
    power = 0
    while power < alphabet_size:
        twos.append(2**power)
        power += 1
    return twos

def gen_prime(alphabet_size):
    primes = []
    count = 0
    for num in range(2, 100):
        if num > 1:
            for i in range(2, num):
                if (num % i) == 0:
                    break
            else:
                primes.append(num)
                count += 1
                if count == alphabet_size:
                    return list(primes)

def create_indet_letter(alphabet_two, alphabet_prime):
    indet_letter_two = 0
    indet_letter_prime = 1
    subset_size = random.randint(2, len(alphabet_two) - 1)
    subset = random.sample(range(len(alphabet_two)), subset_size)
    for i in subset:
        indet_letter_two += alphabet_two[i]
        indet_letter_prime *= alphabet_prime[i]
    return indet_letter_two, indet_letter_prime


def gen_solid_string(length, alphabet_size):
    text_two = []
    text_prime = []
    alphabet_two = gen_twos(alphabet_size)
    alphabet_prime = gen_prime(alphabet_size)
    for i in range(length):
        rand = random.randint(0, alphabet_size - 1)
        text_two.append(alphabet_two[rand])
        text_prime.append(alphabet_prime[rand])
    return text_two, text_prime

def gen_indet_string(text_two, text_prime, alphabet_size, num):
    alphabet_two = gen_twos(alphabet_size)
    alphabet_prime = gen_prime(alphabet_size)
    indet_list = random.sample(range(len(text_two)), num)
    for i in indet_list:
        indet_two, indet_prime = create_indet_letter(alphabet_two, alphabet_prime)
        text_two[i] = indet_two
        text_prime[i] = indet_prime
    return text_two, text_prime

def generate_text_pattern(
    text_length, text_indet_letters, pattern_length, pattern_indet_letters, alphabet_size
):
    text_two, text_prime = gen_solid_string(text_length, alphabet_size)
    gen_indet_string(text_two, text_prime, alphabet_size, text_indet_letters)
    pattern_two, pattern_prime = gen_solid_string(pattern_length, alphabet_size)
    gen_indet_string(pattern_two, pattern_prime, alphabet_size, pattern_indet_letters)

    with open("demofile_pattern_prime.txt", "w") as f:
        f.write(",".join(map(str, pattern_prime)))
    with open("demofile_pattern_twos.txt", "w") as f:
        f.write(",".join(map(str, pattern_two)))
    with open("demofile_text_prime.txt", "w") as f:
        f.write(",".join(map(str, text_prime)))
    with open("demofile_text_twos.txt", "w") as f:
        f.write(",".join(map(str, text_two)))


args = sys.argv[1:]
if len(args) != 5:
    print("Usage: python3 script_name text_length k1 pattern_length k2 sigma")
    sys.exit(1)

val1 = int(args[0])
val2 = int(args[1])
val3 = int(args[2])
val4 = int(args[3])
val5 = int(args[4])

generate_text_pattern(val1, val2, val3, val4, val5)

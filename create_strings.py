import random
from itertools import chain, combinations
import sys


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


def create_indet_letter(alphabet):
    indet_letter = 1
    subset_size = random.randint(2, len(alphabet) - 1)
    subset = random.sample(range(len(alphabet)), subset_size)
    for i in subset:
        indet_letter *= alphabet[i]
    return indet_letter


def gen_solid_string(length, alphabet_size):
    text = []
    alphabet = gen_prime(alphabet_size)
    for i in range(length):
        text.append(alphabet[random.randint(0, alphabet_size - 1)])
    return text


def gen_indet_string(text, alphabet_size, num):
    alphabet = gen_prime(alphabet_size)
    indet_list = random.sample(range(len(text)), num)
    for i in indet_list:
        text[i] = create_indet_letter(alphabet)
    return text


def generate_text_pattern(
    text_length, text_indet_letters, pattern_length, pattern_indet_letters, alphabet
):
    text = gen_solid_string(text_length, alphabet)
    text = gen_indet_string(text, alphabet, text_indet_letters)
    pattern = gen_solid_string(pattern_length, alphabet)
    pattern = gen_indet_string(pattern, alphabet, pattern_indet_letters)
    f = open(
        "demofile_pattern", "w"
    )
    for i in pattern:
        f.write(str(i) + ",")
    f.close()

    f2 = open("demofile_text", "w")
    for i in text:
        f2.write(str(i) + ",")
    f2.close()


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



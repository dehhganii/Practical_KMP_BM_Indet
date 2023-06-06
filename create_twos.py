import random
from itertools import chain, combinations
import sys


def gen_twos(alphabet_size):
    twos = []
    count = 0
    power = 0
    while count < alphabet_size:
        twos.append(2 ** power)
        power += 1
        count += 1
    return twos

def gen_power_set(twos):
    s=list[twos]
    sums=[]
    power_set = list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))
    for i in range(2 ** len(twos)):
        if len(power_set[i]) > 1:
            total_sum = sum(power_set[i])
            sums.append(total_sum)
    return sums

def gen_solid_string(length, alphabet_size):
    text = []
    alphabet = gen_twos(alphabet_size)
    for i in range(length):
        text.append(alphabet[random.randint(0, alphabet_size - 1)])
    return text


def gen_indet_string(text, alphabet_size, num):
    alphabet = gen_twos(alphabet_size)
    indet_list = gen_power_set(alphabet)
    for i in range(0, num):
        pos = random.randint(0, len(text) - 1)
        if text[pos] in alphabet:
            text[pos] = indet_list[random.randint(0, len(indet_list) - 1)]
        else:
            pos = random.randint(0, len(text) - 1)
            text[pos] = indet_list[random.randint(0, len(indet_list) - 1)]
    return text


def gen_indet_pattern(text, alphabet_size, num):
    alphabet = gen_twos(alphabet_size)
    indet_list = gen_power_set(alphabet)
    #text[len(text) - 1] = indet_list[random.randint(0, len(indet_list) - 1)]
    for i in range(0, num):
        pos = random.randint(0, len(text) - 2)
        if text[pos] in alphabet:
            text[pos] = indet_list[random.randint(0, len(indet_list) - 1)]
        else:
            pos = random.randint(0, len(text) - 1)
            text[pos] = indet_list[random.randint(0, len(indet_list) - 1)]
    return text


def generate_text_pattern(
    text_length, text_indet_letters, pattern_length, pattern_indet_letters, alphabet
):
    text = gen_solid_string(text_length, alphabet)
    text = gen_indet_string(text, alphabet, text_indet_letters)
    pattern = gen_solid_string(pattern_length, alphabet)
    pattern = gen_indet_pattern(pattern, alphabet, pattern_indet_letters)
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
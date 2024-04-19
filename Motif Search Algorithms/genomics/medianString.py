import random
import time

def generate_dna_string(length):
    return ''.join(random.choice('ATGC') for _ in range(length))

def mutate_10mer(dna_string):
    index = random.randint(0, len(dna_string) - 10)
    ten_mer = ''.join(random.choices('ATGC', k=10))

    mutated_positions = random.sample(range(10), 4)
    mutated_ten_mer = list(ten_mer)
    for pos in mutated_positions:
        mutated_ten_mer[pos] = random.choice('ATGC'.replace(ten_mer[pos], ''))
    mutated_ten_mer = ''.join(mutated_ten_mer)
    mutated_dna_string = dna_string[:index] + mutated_ten_mer + dna_string[index+10:]
    return mutated_dna_string

dna_strings = [generate_dna_string(500) for _ in range(10)]
mutated_dna_strings = [mutate_10mer(dna_string) for dna_string in dna_strings]


with open('input_file.txt', 'w') as f:
    for dna_string in mutated_dna_strings:
        f.write(dna_string + '\n')


def read_input_file(filename):
    with open(filename, 'r') as f:
        dna_strings = [line.strip() for line in f]
    return dna_strings


def median_string(dna_strings, k):
    def hamming_distance(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    min_distance = float('inf')
    median = None
    for i in range(4 ** k):
        pattern = ''.join(['ACGT'[(i >> (2 * j)) & 3] for j in range(k - 1, -1, -1)])
        distance = sum(min(hamming_distance(pattern, dna[i:i + k]) for i in range(len(dna) - k + 1)) for dna in dna_strings)
        if distance < min_distance:
            min_distance = distance
            median = pattern
    return median

input_file = "input_file.txt"
dna_strings = read_input_file(input_file)

k_values = [9, 10, 11]
for k in k_values:
    print(f"**************************************************************")
    print(f"Running Median String for k={k}")

    start_time_median = time.time()

    median = median_string(dna_strings, k)

    end_time_median = time.time()

    elapsed_time_median = end_time_median - start_time_median

    print(f"Median String: {median}")
    print(f"Time (Median String): {elapsed_time_median:.2f} seconds")

import random
import time
import numpy as np

#rastgele dna oluşturur
def generate_dna_string(length):
    return ''.join(random.choice('ATGC') for _ in range(length))

#10'lu dizi 4'lü mutasyon
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


#konsensüs dizisine olan uzaklığı hesaplar
def score_motifs(motifs):
    consensus = ''
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for base in column:
            count[base] += 1
        consensus += max(count, key=count.get)
    score = sum(1 for motif in motifs for i in range(len(motif)) if motif[i] != consensus[i])
    return consensus, score

#rms algoritması
def randomized_motif_search(dna_strings, k):
    best_motifs = [random.choice([string[i:i + k] for i in range(len(string) - k + 1)]) for string in dna_strings]
    while True:
        profile = form_profile(best_motifs)
        motifs = [profile_most_probable_kmer(string, k, profile) for string in dna_strings]
        if score_motifs(motifs)[1] < score_motifs(best_motifs)[1]:
            best_motifs = motifs
        else:
            return best_motifs


def form_profile(motifs):
    profile = {'A': [], 'C': [], 'G': [], 'T': []}
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        for base in 'ACGT':
            profile[base].append((column.count(base) + 1) / (len(column) + 4))
    return profile

#en iyi k-mer'i bulur
def profile_most_probable_kmer(text, k, profile):
    max_prob = -1
    most_probable_kmer = ''
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        prob = 1
        for j in range(k):
            prob *= profile[pattern[j]][j]
        if prob > max_prob:
            max_prob = prob
            most_probable_kmer = pattern
    return most_probable_kmer

#gibbs algoritması
def gibbs_sampler(dna_strings, k, iterations=1000):
    best_motifs = [random.choice([string[i:i + k] for i in range(len(string) - k + 1)]) for string in dna_strings]
    best_score = score_motifs(best_motifs)[1]
    for _ in range(iterations):
        i = random.randint(0, len(dna_strings) - 1)
        profile = form_profile(best_motifs[:i] + best_motifs[i + 1:])
        motifs_except_i = best_motifs[:i] + best_motifs[i + 1:]
        new_motif = profile_randomly_generated_kmer(dna_strings[i], k, profile)
        motifs_except_i.insert(i, new_motif)
        new_score = score_motifs(motifs_except_i)[1]
        if new_score < best_score:
            best_motifs = motifs_except_i
            best_score = new_score
    return best_motifs


#rastgele k-mer oluşturur
def profile_randomly_generated_kmer(text, k, profile):
    probs = []
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        prob = 1
        for j in range(k):
            prob *= profile[pattern[j]][j]
        probs.append(prob)
    probs = [prob / sum(probs) for prob in probs]
    idx = np.random.choice(len(text) - k + 1, p=probs)
    return text[idx:idx + k]


input_file = "input_file.txt"
dna_strings = read_input_file(input_file)

k_values = [9, 10, 11]
for k in k_values:
    print(f"**************************************************************")
    print(f"Running for k={k}")
    best_motifs_rms = None
    best_motifs_gs = None
    best_score_rms = float('inf')
    best_score_gs = float('inf')
    total_score_rms = 0
    total_score_gs = 0
    num_runs = 5
    for _ in range(num_runs):
        motifs_rms = randomized_motif_search(dna_strings, k)
        motifs_gs = gibbs_sampler(dna_strings, k)
        score_rms = score_motifs(motifs_rms)[1]
        score_gs = score_motifs(motifs_gs)[1]
        total_score_rms += score_rms
        total_score_gs += score_gs
        if score_rms < best_score_rms:
            best_score_rms = score_rms
            best_motifs_rms = motifs_rms
        if score_gs < best_score_gs:
            best_score_gs = score_gs
            best_motifs_gs = motifs_gs
    avg_score_rms = total_score_rms / num_runs
    avg_score_gs = total_score_gs / num_runs


    print(f"--------------------------------------------------------------")
    print(f"*Randomized Motif Search*")
    print(f"Best score (Randomized Motif Search): {best_score_rms}")
    print(f"Average score (Randomized Motif Search): {avg_score_rms}")
    print(f"--------------------------------------------------------------")
    print(f"*Gibbs Sampler*")
    print(f"Best score (Gibbs Sampler): {best_score_gs}")
    print(f"Average score (Gibbs Sampler): {avg_score_gs}")
    print(f"--------------------------------------------------------------")

    #Print the 10 best motifs for Randomized Motif Search
    print("Randomized Motif Search - Best Motifs:")
    for i, motif in enumerate(best_motifs_rms, 1):
        print(f"Motif {i}: {motif}")

    consensus_rms, _ = score_motifs(best_motifs_rms)
    print(f"Consensus string (Randomized Motif Search): {consensus_rms}")

    #Print the 10 best motifs for Gibbs Sampler
    print("\nGibbs Sampler - Best Motifs:")
    for i, motif in enumerate(best_motifs_gs, 1):
        print(f"Motif {i}: {motif}")


    consensus_gs, _ = score_motifs(best_motifs_gs)

    print(f"Consensus string (Gibbs Sampler): {consensus_gs}")
    print()


#Run Randomized Motif Search and Gibbs Sampler algorithms
k_values = [9, 10, 11]
for k in k_values:
    print(f"Running for k={k}")

    start_time_rms = time.time()

    best_motifs_rms = None
    best_score_rms = float('inf')
    total_score_rms = 0
    num_runs = 5
    for _ in range(num_runs):
        motifs_rms = randomized_motif_search(dna_strings, k)
        score_rms = score_motifs(motifs_rms)[1]
        total_score_rms += score_rms
        if score_rms < best_score_rms:
            best_score_rms = score_rms
            best_motifs_rms = motifs_rms
    avg_score_rms = total_score_rms / num_runs


    end_time_rms = time.time()

    elapsed_time_rms = end_time_rms - start_time_rms

    print(f"Time (Randomized Motif Search): {elapsed_time_rms:.2f} seconds")


    start_time_gs = time.time()

    best_motifs_gs = None
    best_score_gs = float('inf')
    total_score_gs = 0
    for _ in range(num_runs):
        motifs_gs = gibbs_sampler(dna_strings, k)
        score_gs = score_motifs(motifs_gs)[1]
        total_score_gs += score_gs
        if score_gs < best_score_gs:
            best_score_gs = score_gs
            best_motifs_gs = motifs_gs
    avg_score_gs = total_score_gs / num_runs

    end_time_gs = time.time()

    elapsed_time_gs = end_time_gs - start_time_gs

    print(f"Time (Gibbs Sampler): {elapsed_time_gs:.2f} seconds")

    print()
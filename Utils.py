import random
import copy
import numpy as np

dna_baz = ['A', 'G', 'T', 'C']

dna_strings = []

#for i in range(0, 10):
#    temp = []
#    for k in range(0, 500):
#        temp.append(dna_baz[random.randint(0, 3)])
#    dna_strings.append(temp)


def save_dna_seq_to_tx(dna=dna_strings):
    file = open("dna_seq.txt", "w")
    for seq in dna:
        file.write("".join(seq))
        file.write("\n")
    file.close()

def read_dna_from_file(file_path):
    with open(file_path) as file:
        for seq in file:
            dna_strings.append(list(seq.replace("\n","")))

#mutation_base = ['A',"A",'A',"A",'A',"A",'A',"A",'A',"A"]
mutation_base = ['A',"C",'G',"T",'A',"A",'C',"A",'G',"T"]
#mutation_base = []
#for i in range(0, 10):
#    mutation_base.append(dna_baz[random.randint(0, 3)])

def get_random_mutation(mutation_b, number_of_mutation):
    rand_idx = np.random.randint(0, len(mutation_b) - 1, size=number_of_mutation)
    temp_dna_baz = ['A', 'G', 'T', 'C']
    for idx in rand_idx:
        baz = mutation_b[idx]
        temp_dna_baz.remove(baz)
        mutation_b[idx] = temp_dna_baz[random.randint(0, 2)]
        temp_dna_baz.append(baz)
    return mutation_b


def apply_mutation(mutation_string, dna_seq):
    k = len(mutation_string)
    idx = random.randint(0, len(dna_seq) - k - 1)
    dna_seq[idx:idx + k] = mutation_string

    return dna_seq


def get_random_k_mer(k, dna):
    idx = random.randint(0, len(dna) - k - 1)
    return dna[idx:idx + k]


def get_motif_matrix(k, dna, ):
    motif_m = []
    for seq_num in range(len(dna)):
        motif_m.append(get_random_k_mer(k, dna[seq_num]))
    return motif_m


def get_baz_matrix(baz, motif_m):
    for seq_index in range(len(motif_m)):
        for char_index in range(len(motif_m[seq_index])):
            if motif_m[seq_index][char_index] == baz:
                motif_m[seq_index][char_index] = 1
            else:
                motif_m[seq_index][char_index] = 0
    return np.asarray(motif_m, dtype=np.int64)


def get_profile(motif_m):
    profile_dict = {}
    for baz_type in dna_baz:
        profile_dict[baz_type] = np.sum(get_baz_matrix(baz_type, copy.deepcopy(motif_m)), axis=0)
    return profile_dict


def increase_by_one(profile_m):
    for seq in profile_m.values():
        for i in range(len(seq)):
            seq[i] = seq[i] + 1
    return profile_m


def apply_pseudo_count(profile_m):
    for seq in profile_m.values():
        for count in seq:
            if count == 0:
                return increase_by_one(profile_m)
                break;
    return profile_m


def get_score(k, profile_m):
    max_val = -1
    score = 0
    for i in range(k):
        for seq in profile_m:
            if seq[i] > max_val:
                max_val = seq[i]
        score = score + (k - max_val)
        max_val = 0
    return score


def get_new_motif_with_prob(k, dna_seq,profile):
    k_mer_dict = {}
    for i in range(0, len(dna_seq) - k + 1):
        k_mer = "".join(dna_seq[i:i + k])
        k_mer_dict[k_mer] = get_prob_for_k_mer(k_mer, profile)
    result = random.choices(list(k_mer_dict.keys()), weights=list(k_mer_dict.values()))[0]
    return result


def get_consensus_string(k,profile_dict):
    max_char = "A"
    max_value = -1
    consensus_string = []
    for i in range(k):
        for key,value in profile_dict.items():
            if value[i] > max_value:
                max_value = value[i]
                max_char = key
        consensus_string.append(max_char)
        max_value = -1
    return "".join(consensus_string)

def get_prob_for_k_mer(k_mer,profile):
    prob = 1

    for baz_idx in range(len(k_mer)):
        baz = k_mer[baz_idx]
        prob = prob*profile[baz][baz_idx]
    return prob

def get_max_k_mer(k,dna_seq,profile):
    k_mer_dict={}
    for i in range(0,len(dna_seq)-k+1):
        k_mer = "".join(dna_seq[i:i + k])
        k_mer_dict[k_mer] = get_prob_for_k_mer(k_mer,profile)
    key = max(k_mer_dict, key= lambda x: k_mer_dict[x])
    return list(key)

def update_motif(k,dna_string,profile,motif_matrix):
    for seq_idx in range(len(dna_string)):
        seq = dna_string[seq_idx]
        motif_matrix[seq_idx] = get_max_k_mer(k,seq,profile)
    return motif_matrix
import random
import string
import numpy as np

from mainVS import run_alg
from edit_distance import edit_distance_ops


def generate_char(pi=0.00, pd=0.00, ps=0.00):
    ops = []
    is_ins = np.random.choice([True, False], p=[pi, 1 - pi])
    is_del = np.random.choice([True, False], p=[pd, 1 - pd])
    is_sub = np.random.choice([True, False], p=[ps, 1 - ps])
    if is_ins:
        c = np.random.choice(['A', 'C', 'G', 'T'])
        ops.append('I' + c)
    if is_del:
        if not is_sub:  # Omer's wish
            ops.append('D')
    if is_sub:
        c = np.random.choice(['A', 'C', 'G', 'T'])
        ops.append('R' + c)
    return ops


# FIXME after last index!!
def generate_copy(orig_str, pi, pd, ps):
    s = ''
    for i in range(len(orig_str)):
        ops = generate_char(pi, pd, ps)
        if ops:
            for op in ops:
                if op[0] == 'I':
                    s += op[1]
                    s += orig_str[i]
                elif op[0] == 'D':
                    continue
                elif op[0] == 'R':
                    s += op[1]
                    continue
        else:
            s += orig_str[i]
    return s


def find_match(s1, s2):
    assert len(s1) == len(s2)
    s = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            s += 1
    return (s / len(s1)) * 100


def generate_test(string_size, cluster_size, pi, pd, ps):
    original_string = ''.join(random.choices('ACGT', k=string_size))
    cluster = [generate_copy(original_string, pi, pd, ps) for i in range(cluster_size)]
    return original_string, cluster

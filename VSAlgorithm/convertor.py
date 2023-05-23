from nltk import edit_distance

from main import run_alg


def from_3ary_to_DNA(s):
    r = ''
    d = {'0': 'A', '1': 'C', '2': 'G', '3': 'T'}
    for c in s:
        r += d[c]
    return r


def from_DNA_to_3ary(s):
    r = ''
    d = {'A': '0', 'C': '1', 'G': '2', 'T': '3'}
    for c in s:
        r += d[c]
    return r


def from_dict_to_list(d, n):
    l = []
    for i in range(n):
        l.append(d[i])
    return l




if __name__ == '__main__':
    with open('input.txt') as f:
        file_s = f.read().strip().replace('\n\n', '\n').split('\n')
        original = from_3ary_to_DNA(file_s[1])
        cluster = [from_3ary_to_DNA(file_s[i]) for i in range(3, 9)]
        res = from_3ary_to_DNA(file_s[10])
        edit_dist = int(file_s[11].replace('Edit Dist ', ''))
        print('original: ', original)
        # print(cluster)
        # print(res)
        # print(edit_dist)
        my_res = run_alg(cluster, 100)
        print('my res  : ', my_res)
        print('Omer res: ', res)
        print('edit dist: my: {}, Omer: {}'.format(edit_distance(my_res, original), edit_dist))
        print('edit distance between results: ', edit_distance(res, my_res))

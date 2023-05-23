from nltk import edit_distance
import sys

from mainVS import run_alg
import ast

banned = []
chunk_i = int(sys.argv[1])
chunk_size = 17
for i in range(chunk_size * chunk_i, chunk_size * (chunk_i + 1)):
    if i + 1 in banned:
        continue
    with open('./tests/test' + str(i) + '.txt', 'r') as inf:
        original_string = inf.readline().strip()[len('original string: '):]
        orig_len = len(original_string)
        g_list = ast.literal_eval(inf.readline().strip()[len('cluster: '):])
        #g_list.sort(key=lambda x: abs(orig_len - len(x)))
        cluster = g_list[:]
        res = run_alg(cluster, orig_len)
        if res[0] == "c":
            with open('./results/result' + str(i) + '.txt', 'w') as resf:
                resf.write(res + "\n")
                resf.write("Cluster        : {}\n".format(str(cluster)))
                continue
        ted = edit_distance(original_string, res)
        with open('./results/result' + str(i) + '.txt', 'w') as resf:
            resf.write("Original String: {}\n".format(original_string))
            resf.write("Cluster        : {}\n".format(str(cluster)))
            resf.write("Result String  : {}\n".format(res))
            resf.write("Edit Distance  : {}\n".format(str(ted)))
            resf.write('\n')

import tests_generator
from edit_distance import edit_distance_ops
from mainVS import run_alg
import sys

if __name__ == '__main__':
    string_size = 100
    #cluster_sizes = [int(sys.argv[1])]
    cluster_sizes = [10] #[6, 10, 20]
    #probabilities = [float(sys.argv[2])]
    probabilities = [0.01] #[0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
    tests_num = 10
    for cluster_size in cluster_sizes:
        for prob in probabilities:
            print("Testing cluster size: {} and prob: {}".format(cluster_size, prob), file=sys.stderr)
            tests_sum = 0
            reg_hist = [0] * 100
            acc_hist = [0] * 100
            # acc_hist = {i: 0 for i in range(100)}
            print(20 * " " + "Testing with {} copies and probability {}".format(cluster_size, prob))
            for test_index in range(1, tests_num + 1):
                original_string, cluster = tests_generator.generate_test(string_size, cluster_size, prob, prob, prob)
                res = run_alg(cluster, string_size)
                ted = (edit_distance_ops(original_string, res)[0])
                reg_hist[int(ted)] += 1
                for hindex in range(int(ted), 100):
                    acc_hist[hindex] += 1
                tests_sum += ted
                #print(original_string)
                #print(res)
                print('.', end="", flush=True)
                with open('./tmp/cs {} prob {}.txt'.format(cluster_size, prob), 'w') as tmpf:
                    tmpf.write('average edit distance: ' + str(tests_sum / test_index) + '\n')
                    tmpf.write('reg histogram:')
                    tmpf.write(str(reg_hist))
                    #exit(0)
            with open('./res/cs {} prob {}.txt'.format(cluster_size, prob), 'w') as resf:
                resf.write('\n')
                resf.write("Average edit distance    : {}".format(tests_sum / tests_num))
                resf.write("Regular histogram is     : {}".format(reg_hist))
                resf.write("Accumulated histogram is : {}".format(acc_hist))
                resf.write('\n')

import tests_generator
from edit_distance import edit_distance_ops
from mainVS import run_alg
import sys

def get_clusters(evya_path):
    file_input = open(evya_path, "r")
    clusters = {}
    reads = []
    label = file_input.readline()  # first label
    label = label.rstrip()
    line = file_input.readline()  # first **
    line = file_input.readline()  # first read
    if len(line) > 10:
        reads.append(line.rstrip())
    while line:
        if len(line) > 10:
            #file_output_reads.write(line)
            line = file_input.readline()
            if len(line) > 10:
                reads.append(line.rstrip())
            continue
        else:
            line = file_input.readline()  # space
            #line = file_input.readline()  # space
            if len(label) >10:
                clusters[label]=reads
                reads = []
            else:
                print("Error")
                exit(0)
            line = file_input.readline()  # label
            if len(line) > 10:
                label = line.rstrip()
                #file_output_labels.write(line)
                line = file_input.readline()  # first **
                #file_output_reads.write("===============================\n")
                line = file_input.readline()  # first read
                if len(line) > 10:
                    reads.append(line.rstrip())
                continue
            else:
                return clusters

if __name__ == '__main__':
    string_size = 110
    #clusters= get_clusters("/Users/omersabary/Dropbox/random_dna_mid_error/datasets_table1/perfect_clusters/Pfister/evyaPfisterPerfect.txt")
    #clusters= get_clusters("/Users/omersabary/Desktop/ReinhardData/fullDS/evyaReinhardPerfect_sortedfilter_2")
    #clusters= get_clusters("/Users/omersabary/Dropbox/random_dna_mid_error/datasets_table1/perfect_clusters/Pfister/evyaPfisterPerfect.txt")
    #clusters= get_clusters("/Users/omersabary/PycharmProjects/pythonProject1/res_Omer/Full_Data_Combined_16/evyaDvirPseudo.txt")
    #clusters= get_clusters("/Users/omersabary/Documents/Microsoft_Data_Inbal_Roy/evyaMicroTechnion_5.txt")
    clusters= get_clusters("/Users/omersabary/PycharmProjects/pythonProject/DataSets_TrellisBMA/NItsan10/evya.txt")

    print(len(clusters.keys()))
    string_size = len(list(clusters.keys())[0])
    #exit(0)
    data_name = "Nitsan10"
    prob = 0.01
    #cluster_sizes = [int(sys.argv[1])]
    #cluster_sizes = [10] #[6, 10, 20]
    #probabilities = [float(sys.argv[2])]
    #probabilities = [0.01] #[0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
    tests_num = 0 #len(clusters.keys())
    tests_sum = 0
    reg_hist = [0] * 100
    acc_hist = [0] * 100
    acc_hist = {i: 0 for i in range(100)}
    for original_string, cluster in clusters.items():
        string_size = len(original_string)
        cluster_size = len(cluster)
        tests_num = tests_num+1
        print(tests_num)
        #for prob in probabilities:
        #    print("Testing cluster size: {} and prob: {}".format(cluster_size, prob), file=sys.stderr)

            #
         #   print(20 * " " + "Testing with {} copies and probability {}".format(cluster_size, prob))
            #for test_index in range(1, tests_num + 1):
        #original_string, cluster = tests_generator.generate_test(string_size, cluster_size, prob, prob, prob)
        res = run_alg(cluster, string_size)
        ted = (edit_distance_ops(original_string, res)[0])
        reg_hist[int(ted)] += 1
        for hindex in range(int(ted), 100):
            acc_hist[hindex] += 1
        tests_sum += ted
        #print(original_string)
        #print(res)
        print('.', end="", flush=True)
        with open('./tmp/cs {} prob {}.txt'.format(data_name, prob), 'w') as tmpf:
            tmpf.write('average edit distance: ' + str(tests_sum / tests_num) + '\n')
            tmpf.write('reg histogram:')
            tmpf.write(str(reg_hist))
            #exit(0)
    with open('./res/cs {} prob {}.txt'.format(data_name, prob), 'w') as resf:
        resf.write('\n')
        resf.write("Average edit distance    : {}".format(tests_sum / (tests_num*string_size)))
        resf.write("Regular histogram is     : {}".format(reg_hist))
        resf.write("Accumulated histogram is : {}".format(acc_hist))
        resf.write('\n')

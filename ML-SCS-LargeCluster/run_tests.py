import subprocess

prob_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.25]
# prob_values = [0.01, 0.02]
count = 1
folder = "ForPaper125_0.25/"
for value in prob_values:
    name_run_file = str(count) + "_125"
    name = "./" + name_run_file
    subprocess.call(["g++","-std=c++11" ,"main.cpp", "reconstruct.cpp", "IndexVectorD.cpp", "utils.cpp", "-o", name_run_file])
    output_file = folder + name_run_file + ".out"
    # print(output_file)

    tmp=subprocess.Popen(["nohup",name , str(value), "&"])
    # print ("printing result")
    # print (tmp)
    count += 1

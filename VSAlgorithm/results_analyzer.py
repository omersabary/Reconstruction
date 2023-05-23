import glob, os, sys

p = sys.argv[1] if len(sys.argv) > 1 else "./results"
b = True

banned = [4490, 4756, 4850, 4879, 4896, 4929, 4937, 4940, 4950, 4955, 4969, 4971, 4973, 4977, 4978, 4979, 4981, 4983, 4985, 4986, 4987]

reg_hist = [0] * 100
acc_hist = [0] * 100
s = 0
i = 0
os.chdir("./results")
for f in glob.glob("*.txt"):
    try:
        if b:
            if (int(f[6:-4])+1) in banned:
                print("Skipping " + f)
                continue
        with open(f, 'r') as file:
            next(file)
            next(file)
            next(file)
            num = int(file.readline().strip()[len('Edit Distance  : '):])
            s += num
            i += 1
            if num >= 100:
                print("FUCK")
                exit()
            reg_hist[num] += 1
            for hindex in range(num, 100):
                acc_hist[hindex] += 1
        if i % 1000 == 0:
            print(i)
    except Exception as e:
        print(e)
        print("f:", f)   
        
print('Average Edit Distance is: ', s / i)
print('Regular Histogram is    : ', reg_hist)
print('Accumulated Histogram is: ', acc_hist)

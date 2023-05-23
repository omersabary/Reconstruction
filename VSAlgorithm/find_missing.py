import glob, os

for f in glob.glob("./tests/*.txt"):
    name = f[len('./tests/test'):]
    if not os.path.isfile('./results/result' + name):
        print('file ./results/result' + name + ' not found!')

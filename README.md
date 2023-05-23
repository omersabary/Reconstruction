# Reconstruction Algorithms for DNA Storage Systems
#### Omer Sabary, Alexander Yucovich, Guy Shapira,  Eitan Yaakobi

In the trace reconstruction problem a length-n string x yields a collection of noisy copies, called traces, y1, …, yt where each yi is independently obtained from x by passing through a deletion channel, which deletes every symbol with some fixed probability. The main goal under this paradigm is to determine the required minimum number of i.i.d traces in order to reconstruct x with high probability. The trace reconstruction problem can be extended to the model where each trace is a result of x passing through a deletion-insertion-substitution channel, which introduces also insertions and substitutions. Motivated by the storage channel of DNA, this work is focused on another variation of the trace reconstruction problem, which is referred by the DNA reconstruction problem. A DNA reconstruction algorithm is a mapping which receives t traces y1, …, yt as an input and produces x^, an estimation of x. The goal in the DNA reconstruction problem is to minimize the edit distance between the original string and the algorithm’s estimation. For the deletion channel case, the problem is referred by the deletion DNA reconstruction problem and the goal is to minimize the Levenshtein distance.

In this work, we present several new algorithms for these reconstruction problems. Our algorithms look globally on the entire sequence of the traces and use dynamic programming algorithms, which are used for the shortest common supersequence and the longest common subsequence problems, in order to decode the original sequence. Our algorithms do not require any limitations on the input and the number of traces, and more than that, they perform well even for error probabilities as high as 0.27. The algorithms have been tested on simulated data as well as on data from previous DNA experiments and are shown to outperform all previous algorithms.
## Algorithms

This repository includes our implementation of the following algorithms:
Algorithms for the deletion-channel:
1. BMA Algorithm. [Batu el al.](http://www.cs.umass.edu/~mcgregor/papers/04-soda.pdf)
2. ML-SCS Algorithm. [Our work](https://www.biorxiv.org/content/10.1101/2020.09.16.300186v1.full)
3. ML-SCS Algorithm for Large Clusters. [Our work](https://www.biorxiv.org/content/10.1101/2020.09.16.300186v1.full)

Algorithms for the deletion-insertion-substitution-channel:
1. BMA-Lookahead Algorithm. [Gopalan et al.](https://patents.google.com/patent/US20180211001A1/en)
2. VS Algorithm. [Viswanathan and Swaminathan](https://dl.acm.org/doi/abs/10.5555/1347082.1347126)
3. The Iterative Reconstruction Algortihm. [Our work](https://www.biorxiv.org/content/10.1101/2020.09.16.300186v1.full)
4. Divider BMA Algorithm. [Our work](https://www.biorxiv.org/content/10.1101/2020.09.16.300186v1.full)
5. Hybrid Algorithm. [Our work](https://www.biorxiv.org/content/10.1101/2020.09.16.300186v1.full)


### Compilation


```bash
g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o *.cpp g++ -o DNA *.o
```


### Usage

```bash
./DNA >results.txt
```

The output file presents edit distance histogram. The i-th entry of the histogram shows the number of clusters that their estimated sequence has edit distance of i errors from the original strings.

The output also includes average edit distance rate (the rate are presented multiplied by 10^-3 for convenience).

### Edlib Installation
Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```bash
pip install edlib
```

### DataSets Links
The two data sets are available here: 
https://drive.google.com/drive/folders/1c3kopMcUsW_tYnjfgjPLuMaLv-cDBV8O?usp=sharing



## Contributing
Pull requests are welcome.
For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
TBA


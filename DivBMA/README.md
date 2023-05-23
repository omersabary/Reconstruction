# Reconstruction Algorithms

Source code of the divider BMA reconstruction algorithm.
 

## Compilation

g++ -std=c++11 DividerBMA.cpp -o DivBMA 


## Usage

./DivBMA >results.txt

The output file presents edit distance histogram. The i-th entry of the histogram shows the number of clusters that their estimated sequence has edit distance of i errors from the original strings. 

The output also includes average edit distance rate (the rate are presented multiplied by 10^-3 for convenience).

## 

## License
TBA

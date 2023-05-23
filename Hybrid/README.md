# Reconstruction Algorithms

Source code of the hybrid algorithm - iterative and divider BMA reconstruction algorithm.
 

## Compilation

  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o LCS2.o LCS2.cpp
  
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o EditDistance.o EditDistance.cpp
  
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Clone.o Clone.cpp
  
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Cluster2.o Cluster2.cpp
  
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o LongestPath.o LongestPath.cpp
  
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o CommonSubstring2.o CommonSubstring2.cpp
  
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o DividerBMA.o DividerBMA.cpp
  
  g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o DNA.o DNA.cpp
  
  g++ -o DNA *.o


## Usage

The input file should be named "evyat.txt": 
Each cluster should be presented in 
Original Strings 
**** 
copy 
copy 
copy 
// Blank row - End of cluster 
// Blank row -  End of cluster 
Original Strings 
**** 
copy 
copy 
copy 
// Blank row - End of cluster 
// Blank row -  End of cluster 

(example file is attached). 

./DNA input_file_full_path path_to_output_files >results.txt

The results can be found in the file output-results.txt  


The output file presents edit distance histogram. The i-th entry of the histogram shows the number of clusters that their estimated sequence has edit distance of i errors from the original strings. 

The output also includes average edit distance rate (the rate are presented multiplied by 10^-3 for convenience).

## 

## License
TBA

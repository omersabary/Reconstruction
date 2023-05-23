CXX = g++
CFLAGS = -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0
RM = rm

# ****************************************************


main: DNA.o LCS2.o EditDistance.o Clone.o Cluster2.o LongestPath.o CommonSubstring2.o DividerBMA.o
	$(CXX) -o main DNA.o LCS2.o EditDistance.o Clone.o Cluster2.o LongestPath.o CommonSubstring2.o
# ****************************************************

DNA.o: DNA.cpp
	$(CXX) $(CFLAGS) -c DNA.cpp

LCS2.o: LCS2.cpp
	$(CXX) $(CFLAGS) -c LCS2.cpp

EditDistance.o: EditDistance.cpp
	$(CXX) $(CFLAGS) -c EditDistance.cpp

Clone.o: Clone.cpp
	$(CXX) $(CFLAGS) -c Clone.cpp

Cluster2.o: Cluster2.cpp
	$(CXX) $(CFLAGS) -c Cluster2.cpp

LongestPath.o: LongestPath.cpp
	$(CXX) $(CFLAGS) -c LongestPath.cpp

CommonSubstring2.o: CommonSubstring2.cpp
	$(CXX) $(CFLAGS) -c CommonSubstring2.cpp

DividerBMA.o: DividerBMA.cpp
	$(CXX) $(CFLAGS) -c DividerBMA.cpp


clean:
	$(RM) *.o


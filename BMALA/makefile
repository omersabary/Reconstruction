CXX = g++
CFLAGS = -std=c++0x -O3 -g3 -c -fmessage-length=0
RM = rm

# ****************************************************


main: BMALookahead.o
	$(CXX) -o main BMALookahead.o
# ****************************************************

BMALookahead.o : BMALookahead.cpp
	$(CXX) $(CFLAGS) -c BMALookahead.cpp




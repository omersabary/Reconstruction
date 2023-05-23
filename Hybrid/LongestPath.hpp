/*
 * LongestPath.hpp
 *
 *  Created on: 8 αιεμι 2019
 *      Author: Alex
 */

#ifndef LONGESTPATH_HPP_
#define LONGESTPATH_HPP_

// A C++ program to find single source longest distances
// in a DAG
#include <iostream>
#include <limits.h>
#include <list>
#include <stack>
#include <vector>
#include <unordered_map>
#include <map>
#define NINF INT_MIN
using namespace std;

// Class to represent a graph using adjacency list
// representation
class Graph {
	int V; // No. of vertices'
	std::unordered_map<int, std::vector<pair<int,int> > > graph;
	// A function used by longestPath
	void TopologicalSortUtil(int v, vector<bool>& visited, stack<int>& Stack);

public:
	Graph(int V); // Constructor

	// function to add an edge to graph
	void AddEdge(const int u, const int v, const int weight);
	int Weight(const int u, const int v) const;
	int GetV() const;
	// Finds longest distances from given source vertex
	std::vector<int> LongestPathLen(int s);
	std::vector<int> LongestPath(int u, int v);
};

std::vector<std::string> DictLongestPath(std::vector<std::map<std::string, int> >& listVec, const int patternLen);

#endif /* LONGESTPATH_HPP_ */

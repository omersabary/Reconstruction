#include "LongestPath.hpp"
#include <cassert>

const int NA = -1;

Graph::Graph(int V) :
		V(V), graph() {
}

void Graph::AddEdge(const int u, const int v, const int weight) {
	graph[u].push_back(make_pair(v, weight));
}

int Graph::Weight(const int u, const int v) const {
	unordered_map<int, vector<pair<int, int> > >::const_iterator uIt = graph.find(u);
	if (uIt == graph.end()) {
		return NA;
	}
	else {
		for (vector<pair<int, int> >::const_iterator it = uIt->second.begin(); it != uIt->second.end(); it++) {
			if (it->first == v) {
				return it->second;
			}
		}
		return NA;
	}
}

int Graph::GetV() const {
	return V;
}

// A recursive function used by longestPath. See below
// link for details
// https:// www.geeksforgeeks.org/topological-sorting/
void Graph::TopologicalSortUtil(int v, vector<bool>& visited, stack<int>& Stack) {
	// Mark the current node as visited
	visited[v] = true;

	// Recur for all the vertices adjacent to this vertex
	unordered_map<int, vector<pair<int, int> > >::iterator vIt = graph.find(v);
	if (vIt != graph.end()) {
		for (vector<pair<int, int> >::iterator it = vIt->second.begin(); it != vIt->second.end(); it++) {
			if (not visited[it->first]) {
				TopologicalSortUtil(it->first, visited, Stack);
			}
		}
	}

	// Push current vertex to stack which stores topological sort
	Stack.push(v);
}

// The function to find longest distances from a given vertex.
// It uses recursive topologicalSortUtil() to get topological
// sorting.
vector<int> Graph::LongestPathLen(int s) {
	stack<int> Stack;
	vector<bool> visited(V, false);
	vector<int> dist(V, NINF);
	vector<int> parents(V, NA);

	// Call the recursive helper function to store Topological
	// Sort starting from all vertices one by one
	for (int i = 0; i < V; i++)
		if (visited[i] == false)
			// TODO: s or i here? Answer: i and s are the same. because one call for Topological sort on s is enough.
			TopologicalSortUtil(i, visited, Stack);

	// distance to source as 0
	dist[s] = 0;

	// Process vertices in topological order
	while (Stack.empty() == false) {
		// Get the next vertex from topological order
		int u = Stack.top();
		Stack.pop();

		// Update distances of all adjacent vertices
		if (dist[u] != NINF) {
			for (vector<pair<int, int> >::iterator it = graph[u].begin(); it != graph[u].end(); it++) {
				if (dist[it->first] < dist[u] + it->second) {
					dist[it->first] = dist[u] + it->second;
					parents[it->first] = u;
				}
			}
		}
	}
	return parents;
}

vector<int> Graph::LongestPath(int u, int v) {
	vector<int> path;
	vector<int> parents = LongestPathLen(u);
	if (parents[v] == NA) {
		return path;
	}
	else {
		path.push_back(v);
		int parent = parents[v];
		while (parent != NA) {
			path.insert(path.begin(), parent);
			parent = parents[parent];
		}
		return path;
	}
}

int CountUpper(const string& str) {
	int upperNum = 0;
	for (unsigned index = 0; index < str.size(); index++) {
		if (isupper(str[index])) {
			upperNum++;
		}
	}
	return upperNum;
}

bool HasEdge(const string& fromKey, const string& toKey, const int patternLen) {
	assert(not fromKey.empty());
	assert(not toKey.empty());
	assert(isupper(fromKey[0]));
	assert(isupper(fromKey.back()));
	int fromPatternSize = CountUpper(fromKey);
	int toPatternSize = CountUpper(toKey);
	string fromSuffix, toPrefix;
	if (fromPatternSize < toPatternSize) {
		fromSuffix = fromKey;
	}
	else {
		int fromIndex = 1;
		while (not isupper(fromKey[fromIndex])) {
			fromIndex++;
		}
		fromSuffix = fromKey.substr(fromIndex);
	}
	if (toPatternSize < fromPatternSize) {
		toPrefix = toKey;
	}
	else {
		int toIndex = toKey.size() - 2;
		while (not isupper(toKey[toIndex])) {
			toIndex--;
		}
		toPrefix = toKey.substr(0, toIndex + 1);
	}
	if (fromSuffix == toPrefix) {
		return true;
	}
	else {
		return false;
	}
}

Graph DictToGraph(vector<map<string, int> >& listVec, const int patternLen, vector<pair<int, string> >& vertexToKey) {
	int dictionariesNum = listVec.size();
	vector<vector<string> > keysLists(dictionariesNum);
	vector<vector<int> > keysToVertex(dictionariesNum);
	vertexToKey.push_back(make_pair(NA, "")); // start vertex
	int vIndex = 1;
	for (int dicIndex = 0; dicIndex < dictionariesNum; dicIndex++) {
		for (map<string, int>::const_iterator it = listVec[dicIndex].begin(); it != listVec[dicIndex].end(); it++) {
			keysLists[dicIndex].push_back(it->first);
			vertexToKey.push_back(make_pair(dicIndex, it->first));
			keysToVertex[dicIndex].push_back(vIndex);
			vIndex++;
		}
	}
	vertexToKey.push_back(make_pair(NA, "")); // end vertex
	int vertexNum = vertexToKey.size();
	Graph g(vertexNum);

	//add edges from start vertex(0) to all vertices of first dictionary with weight 0.
	for (unsigned keyIndex = 0; keyIndex < keysToVertex[0].size(); keyIndex++) {
		g.AddEdge(0, keysToVertex[0][keyIndex], 0);
	}

	//from each vertex in level i to each vertex in level i+1 add edge if has_edge is true
	for (int level = 0; level < dictionariesNum - 1; level++) {
		for (unsigned fromKeyIndex = 0; fromKeyIndex < keysToVertex[level].size(); fromKeyIndex++) {
			string fromKey = keysLists[level][fromKeyIndex];
			for (unsigned toKeyIndex = 0; toKeyIndex < keysToVertex[level + 1].size(); toKeyIndex++) {
				string toKey = keysLists[level + 1][toKeyIndex];
				if (HasEdge(fromKey, toKey, patternLen)) {
					g.AddEdge(keysToVertex[level][fromKeyIndex], keysToVertex[level + 1][toKeyIndex],
							listVec[level][fromKey]);
				}
			}
		}
	}

	//edges from last level to end vertex with their respective weights
	for (unsigned keyIndex = 0; keyIndex < keysToVertex[dictionariesNum - 1].size(); keyIndex++) {
		string key = keysLists[dictionariesNum - 1][keyIndex];
		g.AddEdge(keysToVertex[dictionariesNum - 1][keyIndex], vertexNum - 1, listVec[dictionariesNum - 1][key]);
	}
	return g;
}

vector<string> DictLongestPath(vector<map<string, int> >& listVec, const int patternLen) {
	vector<pair<int, string> > vertexToKey;
	Graph g = DictToGraph(listVec, patternLen, vertexToKey);
	vector<int> fullPath = g.LongestPath(0, g.GetV() - 1);
	assert(fullPath.size()>0);
	vector<string> path;
	for (unsigned index = 1; index < fullPath.size() - 1; index++) {
		path.push_back(vertexToKey[fullPath[index]].second);
	}
	assert(path.size() > 0);
	return path;
}

/*// Driver program to test above functions
 int main()
 {
 // Create a graph given in the above diagram.
 // Here vertex numbers are 0, 1, 2, 3, 4, 5 with
 // following mappings:
 // 0=r, 1=s, 2=t, 3=x, 4=y, 5=z
 Graph g(6);
 g.addEdge(0, 1, 5);
 g.addEdge(0, 2, 3);
 g.addEdge(1, 3, 6);
 g.addEdge(1, 2, 2);
 g.addEdge(2, 4, 4);
 g.addEdge(2, 5, 2);
 g.addEdge(2, 3, 7);
 g.addEdge(3, 5, 1);
 g.addEdge(3, 4, -1);
 g.addEdge(4, 5, -2);

 int s = 1;
 cout << "Following are longest distances from "
 "source vertex "
 << s << " \n";
 g.longestPath(s);

 return 0;
 }*/

#include "SCS2.hpp"
#include "Strings.hpp"
#include <cassert>
#include <map>
#include <unordered_map>
#include <iostream>
#include <algorithm>
using namespace std;

// Function to return all SCS of substrings X[0..m-1], Y[0..n-1]
vector<string> FindSCS(const string& X, const string& Y, int m, int n, int** lookup) {
	// if we have reached the end of first string, create a vector containing second substring and return
	if (m == 0) {
		vector<string> v;
		v.push_back(Y.substr(0, n));
		return v;
	}

	// if we have reached the end of second string, create a vector containing first substring and return
	else if (n == 0) {
		vector<string> v;
		v.push_back(X.substr(0, m));
		return v;
	}

	// if last character of X and Y is same, it should occur only one time in SCS
	if (X[m - 1] == Y[n - 1]) {
		// find all SCS of substring X[0..m-2], Y[0..n-2]
		vector<string> scs = FindSCS(X, Y, m - 1, n - 1, lookup);

		// append current character X[m - 1] or Y[n - 1] to all SCS of substring X[0..m-2] and Y[0..n-2]

		for (vector<string>::iterator str = scs.begin(); str != scs.end(); str++)
			str->push_back(X[m - 1]);

		return scs;
	}

	// we reach here when the last character of X and Y don't match

	// if top cell of current cell has less value than the left cell, then append current character of string X to all
	// SCS of substring X[0..m-2], Y[0..n-1]

	if (lookup[m - 1][n] < lookup[m][n - 1]) {
		vector<string> scs = FindSCS(X, Y, m - 1, n, lookup);

		for (vector<string>::iterator str = scs.begin(); str != scs.end(); str++)
			str->push_back(X[m - 1]);

		return scs;
	}

	// if left cell of current cell has less value than the top cell, then append current character of string Y to all
	// SCS of substring X[0..m-1], Y[0..n-2]

	if (lookup[m][n - 1] < lookup[m - 1][n]) {
		vector<string> scs = FindSCS(X, Y, m, n - 1, lookup);

		for (vector<string>::iterator str = scs.begin(); str != scs.end(); str++)
			str->push_back(Y[n - 1]);

		return scs;
	}

	// if top cell value is same as left cell, then go in both top and left directions

	// append current character of string X to all SCS of substring X[0..m-2], Y[0..n-1]
	vector<string> top = FindSCS(X, Y, m - 1, n, lookup);
	for (vector<string>::iterator str = top.begin(); str != top.end(); str++)
		str->push_back(X[m - 1]);

	// append current character of string Y to all SCS of substring X[0..m-1], Y[0..n-2]
	vector<string> left = FindSCS(X, Y, m, n - 1, lookup);
	for (vector<string>::iterator str = left.begin(); str != left.end(); str++)
		str->push_back(Y[n - 1]);

	// merge two vectors and return
	top.insert(top.end(), left.begin(), left.end());

	return top;
}

string FindOneSCSWildCards(const string& X, const string& Y, int m, int n, int** lookup) {
	// if we have reached the end of first string, create a vector containing second substring and return
	if (m == 0) {
		return string(n, 'Y');
	}

	// if we have reached the end of second string, create a vector containing first substring and return
	else if (n == 0) {
		return string(m, 'X');
	}

	// if last character of X and Y is same, it should occur only one time in SCS
	if (X[m - 1] == Y[n - 1]) {
		string scs = FindOneSCSWildCards(X, Y, m - 1, n - 1, lookup);
		scs.push_back(X[m - 1]);
		return scs;
	}

	// we reach here when the last character of X and Y don't match

	// if top cell of current cell has less value than the left cell, then append current character of string X to all
	// SCS of substring X[0..m-2], Y[0..n-1]

	if (lookup[m - 1][n] <= lookup[m][n - 1]) {
		string scs = FindOneSCSWildCards(X, Y, m - 1, n, lookup);
		scs.push_back('X');

		return scs;
	}

	// if left cell of current cell has less value than the top cell, then append current character of string Y to all
	// SCS of substring X[0..m-1], Y[0..n-2]

	else { //lookup[m][n - 1] < lookup[m - 1][n]
		string scs = FindOneSCSWildCards(X, Y, m, n - 1, lookup);
		scs.push_back('Y');

		return scs;
	}
}

string FindOneSCSWildCardsRandom(const string& X, const string& Y, int m, int n, int** lookup, mt19937& generator,
		uniform_int_distribution<int>& dist) {
	// if we have reached the end of first string, create a vector containing second substring and return
	if (m == 0) {
		return string(n, 'Y');
	}

	// if we have reached the end of second string, create a vector containing first substring and return
	else if (n == 0) {
		return string(m, 'X');
	}

	// if last character of X and Y is same, it should occur only one time in SCS
	if (X[m - 1] == Y[n - 1]) {
		string scs = FindOneSCSWildCardsRandom(X, Y, m - 1, n - 1, lookup, generator, dist);
		scs.push_back(X[m - 1]);
		return scs;
	}

	// we reach here when the last character of X and Y don't match

	// if top cell of current cell has less value than the left cell, then append current character of string X to all
	// SCS of substring X[0..m-2], Y[0..n-1]

	if (dist(generator)) {
		if (lookup[m - 1][n] <= lookup[m][n - 1]) {
			string scs = FindOneSCSWildCardsRandom(X, Y, m - 1, n, lookup, generator, dist);
			scs.push_back('X');

			return scs;
		}

		// if left cell of current cell has less value than the top cell, then append current character of string Y to all
		// SCS of substring X[0..m-1], Y[0..n-2]

		else { //lookup[m][n - 1] < lookup[m - 1][n]
			string scs = FindOneSCSWildCardsRandom(X, Y, m, n - 1, lookup, generator, dist);
			scs.push_back('Y');

			return scs;
		}
	}
	else {
		if (lookup[m - 1][n] < lookup[m][n - 1]) {
			string scs = FindOneSCSWildCardsRandom(X, Y, m - 1, n, lookup, generator, dist);
			scs.push_back('X');

			return scs;
		}

		// if left cell of current cell has less value than the top cell, then append current character of string Y to all
		// SCS of substring X[0..m-1], Y[0..n-2]

		else { //lookup[m][n - 1] < lookup[m - 1][n]
			string scs = FindOneSCSWildCardsRandom(X, Y, m, n - 1, lookup, generator, dist);
			scs.push_back('Y');

			return scs;
		}
	}
}

// Function to fill lookup table by finding length of SCS of
// sequences X[0..m-1] and Y[0..n-1]
int SCSLength(const string& X, const string& Y, int m, int n, int** lookup) {

	// initialize first column of the lookup table
	for (int i = 0; i <= m; i++)
		lookup[i][0] = i;

	// initialize first row of the lookup table
	for (int j = 0; j <= n; j++)
		lookup[0][j] = j;

	// fill the lookup table in bottom-up manner
	for (int i = 1; i <= m; i++) {
		for (int j = 1; j <= n; j++) {
			// if current character of X and Y matches
			if (X[i - 1] == Y[j - 1])
				lookup[i][j] = lookup[i - 1][j - 1] + 1;

			// else if current character of X and Y don't match
			else
				lookup[i][j] = min(lookup[i - 1][j] + 1, lookup[i][j - 1] + 1);
		}
	}
	return lookup[m][n];
}

// Function to find all Shortest Common Supersequence of string X and Y
vector<string> SCS2(const string& X, const string& Y) {
	int m = X.length(), n = Y.length();
	int **lookup;
	lookup = new int*[m + 1];
	for (int i = 0; i < m + 1; i++)
		lookup[i] = new int[n + 1];

	// fill lookup table
	SCSLength(X, Y, m, n, lookup);

	// find all Shortest Common Supersequence
	vector<string> v = FindSCS(X, Y, m, n, lookup);

	// since vector can contain duplicates, "copy" the vector to set
	//set<string> scs(v.begin(), v.end());

	// release memory
	for (int i = 0; i < m + 1; i++)
		delete[] lookup[i];
	delete[] lookup;
	// return set containing all SCS
	return v;
}

int EstimateSCSNum(const string& scsWithWildCards) {
	int strLen = scsWithWildCards.size();
	int score = 1;
	int xnum = 0, ynum = 0;
	for (int i = 0; i < strLen; i++) {
		if (scsWithWildCards[i] == 'X') {
			xnum++;
		}
		else if (scsWithWildCards[i] == 'Y') {
			ynum++;
		}
		else {
			vector<string> gapSCS = SCS2(string(xnum, 'X'), string(ynum, 'Y'));
			score *= gapSCS.size();
			xnum = 0;
			ynum = 0;
		}
	}
	// handle last gap (if last letter is X or Y)
	vector<string> gapSCS = SCS2(string(xnum, 'X'), string(ynum, 'Y'));
	score *= gapSCS.size();

	return score;
}

string SCS2WithWildCardsRandom(const string& X, const string& Y, mt19937& generator) {
	int m = X.length(), n = Y.length();
	int **lookup;
	lookup = new int*[m + 1];
	for (int i = 0; i < m + 1; i++)
		lookup[i] = new int[n + 1];

	// fill lookup table
	SCSLength(X, Y, m, n, lookup);
	uniform_int_distribution<int> distribution(0, 1);
	// find all Shortest Common Supersequence
	string v = FindOneSCSWildCardsRandom(X, Y, m, n, lookup,generator,distribution);

	// release memory
	for (int i = 0; i < m + 1; i++)
		delete[] lookup[i];
	delete[] lookup;
	// return scs
	return v;
}

string SCS2WithWildCards(const string& X, const string& Y) {
	int m = X.length(), n = Y.length();
	int **lookup;
	lookup = new int*[m + 1];
	for (int i = 0; i < m + 1; i++)
		lookup[i] = new int[n + 1];

	// fill lookup table
	SCSLength(X, Y, m, n, lookup);

	// find all Shortest Common Supersequence
	string v = FindOneSCSWildCards(X, Y, m, n, lookup);

	// release memory
	for (int i = 0; i < m + 1; i++)
		delete[] lookup[i];
	delete[] lookup;
	// return scs
	return v;
}

//vector<string> MostLikelySCS2(const string& X, const string& Y) {
//	int m = X.length(), n = Y.length();
//	int **lookup;
//	lookup = new int*[m + 1];
//	for (int i = 0; i < m + 1; i++)
//		lookup[i] = new int[n + 1];
//
//	// fill lookup table
//	SCSLength(X, Y, m, n, lookup);
//
//	// find all Shortest Common Supersequence
//	vector<string> v = FindMostLikelySCS2(X, Y, m, n, lookup);
//
//	// release memory
//	for (int i = 0; i < m + 1; i++)
//		delete[] lookup[i];
//	delete[] lookup;
//	return v;
//}

int SCS2Len(const string& X, const string& Y) {
	int m = X.length(), n = Y.length();
	int **lookup;
	lookup = new int*[m + 1];
	for (int i = 0; i < m + 1; i++)
		lookup[i] = new int[n + 1];

	// fill lookup table
	int len = SCSLength(X, Y, m, n, lookup);

	// release memory
	for (int i = 0; i < m + 1; i++)
		delete[] lookup[i];
	delete[] lookup;
	// return set containing all SCS
	return len;
}

int SCSLengthFast2(const string& X, const string& Y, int m, int n, int d1, int d2, vector<vector<int>>& L,
		vector<int>& mLowerBounds, vector<int>& mUpperBounds) {

	// fill the lookup table in bottom-up manner
	int maxi = min(m, n + d2), previousLowerBound = 0, previousUpperBound = 0;
	mLowerBounds = vector<int>(maxi + 1);
	mUpperBounds = vector<int>(maxi + 1);
	L = vector<vector<int>>(maxi + 1);
	for (int i = 0; i <= maxi; i++) {
		int minj = max(0, i - d2);
		int maxj = min(i + d1, n);
		mLowerBounds[i] = minj;
		mUpperBounds[i] = maxj;
		L[i] = vector<int>(maxj - minj + 1);
		for (int j = minj; j <= maxj; j++) {
			if (i == 0) {
				L[i][j - minj] = j;
			}
			else if (j == 0) {
				L[i][j - minj] = i;
			}
			// if current character of X and Y matches
			else if (X[i - 1] == Y[j - 1]) {
				//L[i][j - minj] = L[i - 1][j - 1 - mLowerBounds[i - 1]] + 1;
				L[i][j - minj] = L[i - 1][j - 1 - previousLowerBound] + 1;
			}

			// else if current character of X and Y don't match
			else {
				if (j > previousUpperBound) {	//can't go up
					L[i][j - minj] = L[i][j - 1 - minj] + 1;
				}
				else if (j - 1 < minj) {	// can't go left
					L[i][j - minj] = L[i - 1][j - previousLowerBound] + 1;
				}
				else {
					L[i][j - minj] = min(L[i][j - 1 - minj], L[i - 1][j - previousLowerBound]) + 1;
				}
			}
		}
		previousLowerBound = minj;
		previousUpperBound = maxj;
	}
	return L[maxi][n - mLowerBounds[maxi]];
}

int SCSLengthFast2(const string& X, const string& Y, int m, int n, int d1, int d2, vector<vector<int>>& L) {

	// fill the lookup table in bottom-up manner
	int maxi = min(m, n + d2), previousLowerBound = 0, previousUpperBound = 0;
	L = vector<vector<int>>(maxi + 1);
	for (int i = 0; i <= maxi; i++) {
		int minj = max(0, i - d2);
		int maxj = min(i + d1, n);
		L[i] = vector<int>(maxj - minj + 1);
		for (int j = minj; j <= maxj; j++) {
			if (i == 0) {
				L[i][j - minj] = j;
			}
			else if (j == 0) {
				L[i][j - minj] = i;
			}
			// if current character of X and Y matches
			else if (X[i - 1] == Y[j - 1]) {
				L[i][j - minj] = L[i - 1][j - 1 - previousLowerBound] + 1;
			}

			// else if current character of X and Y don't match
			else {
				if (j > previousUpperBound) {	//can't go up
					L[i][j - minj] = L[i][j - 1 - minj] + 1;
				}
				else if (j - 1 < minj) {	// can't go left
					L[i][j - minj] = L[i - 1][j - previousLowerBound] + 1;
				}
				else {
					L[i][j - minj] = min(L[i][j - 1 - minj], L[i - 1][j - previousLowerBound]) + 1;
				}
			}
		}
		previousLowerBound = minj;
		previousUpperBound = maxj;
	}
	return L[maxi][n - previousLowerBound];
}

struct Node {
	vector<string> strs;
	vector<Node*> previous;
	void Add(Node* ptr, const string& str) {
		strs.push_back(str);
		previous.push_back(ptr);
	}
	void Add(Node* ptr, const char letter) {
		strs.push_back(string(1, letter));
		previous.push_back(ptr);
	}
	int size() const {
		return strs.size();
	}
	const string& FirstString() const {
		return strs[0];
	}
	Node* FirstPointer() const {
		return previous[0];
	}
	const string& SecondString() const {
		return strs[1];
	}
	Node* SecondPointer() const {
		return previous[1];
	}
};

int SCSLengthFastest(const string& X, const string& Y, int m, int n, int d1, int d2, vector<vector<Node>>& B,
		Node** startNode) {

	// fill the lookup table in bottom-up manner
	int maxi = min(m, n + d2), previousLowerBound = 0, previousUpperBound = 0;
	vector<vector<int>> L(maxi + 1);
	B = vector<vector<Node>>(maxi + 1);
	for (int i = 0; i <= maxi; i++) {
		int minj = max(0, i - d2);
		int maxj = min(i + d1, n);
		L[i] = vector<int>(maxj - minj + 1);
		B[i] = vector<Node>(maxj - minj + 1);
		for (int j = minj; j <= maxj; j++) {
			if (i == 0) {
				L[i][j - minj] = j;
				B[i][j - minj].Add(NULL, Y.substr(0, j));
			}
			else if (j == 0) {
				L[i][j - minj] = i;
				B[i][j - minj].Add(NULL, X.substr(0, i));
			}
			// if current character of X and Y matches
			else if (X[i - 1] == Y[j - 1]) {
				L[i][j - minj] = L[i - 1][j - 1 - previousLowerBound] + 1;
				B[i][j - minj].Add(&B[i - 1][j - 1 - previousLowerBound], X[i - 1]);
			}

			// else if current character of X and Y don't match
			else {
				if (j > previousUpperBound) {	//can't go up
					L[i][j - minj] = L[i][j - 1 - minj] + 1;
					B[i][j - minj].Add(&B[i][j - 1 - minj], Y[j - 1]);
				}
				else if (j - 1 < minj) {	// can't go left
					L[i][j - minj] = L[i - 1][j - previousLowerBound] + 1;
					B[i][j - minj].Add(&B[i - 1][j - previousLowerBound], X[i - 1]);
				}
				else {	// can go in both directions
					int Lup = L[i - 1][j - previousLowerBound];
					int Lleft = L[i][j - 1 - minj];
					if (Lup < Lleft) { // go up
						L[i][j - minj] = Lup + 1;
						B[i][j - minj].Add(&B[i - 1][j - previousLowerBound], X[i - 1]);
					}
					else if (Lup > Lleft) { // go left
						L[i][j - minj] = Lleft + 1;
						B[i][j - minj].Add(&B[i][j - 1 - minj], Y[j - 1]);
					}
					else { // equal. go in both directions
						L[i][j - minj] = Lup + 1;
						B[i][j - minj].Add(&B[i - 1][j - previousLowerBound], X[i - 1]);
						B[i][j - minj].Add(&B[i][j - 1 - minj], Y[j - 1]);
					}
				}
			}
		}
		previousLowerBound = minj;
		previousUpperBound = maxj;
	}
	*startNode = &B[maxi][n - previousLowerBound];
	return L[maxi][n - previousLowerBound];
}

vector<int> BacktrackSCS(const Node* startNode, vector<string>& result, int scsLen) {
	assert(startNode!=NULL);
	Node* firstPtr = startNode->FirstPointer();
	const string& firstStr = startNode->FirstString();
	if (firstPtr == NULL) {
		result.push_back(firstStr);
		result.back().reserve(scsLen);
		return vector<int>( { (int) result.size() - 1 });
	}
	else if (startNode->size() == 1) {
		vector<int> res = BacktrackSCS(firstPtr, result, scsLen);
		for (vector<int>::iterator it = res.begin(); it != res.end(); it++) {
			result[*it] += firstStr;
		}
		return res;
	}
	else { // two paths
		vector<int> res1 = BacktrackSCS(firstPtr, result, scsLen);
		for (vector<int>::iterator it = res1.begin(); it != res1.end(); it++) {
			result[*it] += firstStr;
		}
		Node* secondPtr = startNode->SecondPointer();
		const string& secondStr = startNode->SecondString();
		vector<int> res2 = BacktrackSCS(secondPtr, result, scsLen);
		for (vector<int>::iterator it = res2.begin(); it != res2.end(); it++) {
			result[*it] += secondStr;
		}
		res1.insert(res1.end(), res2.begin(), res2.end());
		return res1;
	}
}

vector<string> FindSCSFast2(const string& X, const string& Y, int m, int n, const vector<vector<int>>& L,
		const vector<int>& mLowerBounds, const vector<int>& mUpperBounds) {
	// if we have reached the end of first string, create a vector containing second substring and return
	if (m == 0) {
		vector<string> v;
		v.push_back(Y.substr(0, n));
		return v;
	}

	// if we have reached the end of second string, create a vector containing first substring and return
	else if (n == 0) {
		vector<string> v;
		v.push_back(X.substr(0, m));
		return v;
	}

	// if last character of X and Y is same, it should occur only one time in SCS
	if (X[m - 1] == Y[n - 1]) {
		// find all SCS of substring X[0..m-2], Y[0..n-2]
		vector<string> scs = FindSCSFast2(X, Y, m - 1, n - 1, L, mLowerBounds, mUpperBounds);

		// append current character X[m - 1] or Y[n - 1] to all SCS of substring X[0..m-2] and Y[0..n-2]

		for (vector<string>::iterator str = scs.begin(); str != scs.end(); str++)
			str->push_back(X[m - 1]);

		return scs;
	}

	if (n - 1 < mLowerBounds[m]) {	// can't go left
		vector<string> scs = FindSCSFast2(X, Y, m - 1, n, L, mLowerBounds, mUpperBounds);

		for (vector<string>::iterator str = scs.begin(); str != scs.end(); str++)
			str->push_back(X[m - 1]);

		return scs;
	}

	if (n > mUpperBounds[m - 1]) {	// can't go up
		vector<string> scs = FindSCSFast2(X, Y, m, n - 1, L, mLowerBounds, mUpperBounds);

		for (vector<string>::iterator str = scs.begin(); str != scs.end(); str++)
			str->push_back(Y[n - 1]);

		return scs;
	}

	int leftL = L[m][n - 1 - mLowerBounds[m]];
	int upL = L[m - 1][n - mLowerBounds[m - 1]];

	// we reach here when the last character of X and Y don't match

	// if top cell of current cell has less value than the left cell, then append current character of string X to all
	// SCS of substring X[0..m-2], Y[0..n-1]

	if (upL < leftL) {
		vector<string> scs = FindSCSFast2(X, Y, m - 1, n, L, mLowerBounds, mUpperBounds);

		for (vector<string>::iterator str = scs.begin(); str != scs.end(); str++)
			str->push_back(X[m - 1]);

		return scs;
	}

	// if left cell of current cell has less value than the top cell, then append current character of string Y to all
	// SCS of substring X[0..m-1], Y[0..n-2]

	if (leftL < upL) {
		vector<string> scs = FindSCSFast2(X, Y, m, n - 1, L, mLowerBounds, mUpperBounds);

		for (vector<string>::iterator str = scs.begin(); str != scs.end(); str++)
			str->push_back(Y[n - 1]);

		return scs;
	}

	// if top cell value is same as left cell, then go in both top and left directions
	// append current character of string X to all SCS of substring X[0..m-2], Y[0..n-1]
	vector<string> top = FindSCSFast2(X, Y, m - 1, n, L, mLowerBounds, mUpperBounds);
	for (vector<string>::iterator str = top.begin(); str != top.end(); str++)
		str->push_back(X[m - 1]);

	// append current character of string Y to all SCS of substring X[0..m-1], Y[0..n-2]
	vector<string> left = FindSCSFast2(X, Y, m, n - 1, L, mLowerBounds, mUpperBounds);
	for (vector<string>::iterator str = left.begin(); str != left.end(); str++)
		str->push_back(Y[n - 1]);

	// merge two vectors and return
	top.insert(top.end(), left.begin(), left.end());

	return top;
}

vector<string> SCS2Fast(const string& X, const string& Y, const int originalLen) {
	int m = X.length(), n = Y.length();
	int d1 = originalLen - m, d2 = originalLen - n;
	vector<vector<int>> L;
	vector<int> mLowerBounds;
	vector<int> mUpperBounds;

// fill lookup table
	SCSLengthFast2(X, Y, m, n, d1, d2, L, mLowerBounds, mUpperBounds);

// find all Shortest Common Supersequence
	return FindSCSFast2(X, Y, m, n, L, mLowerBounds, mUpperBounds);
}

vector<string> SCS2Fastest(const string& X, const string& Y, const int originalLen) {
	int m = X.length(), n = Y.length();
	int d1 = originalLen - m, d2 = originalLen - n;
	vector<vector<Node>> B;

// fill lookup table
	Node* startNode;
	int scsLen = SCSLengthFastest(X, Y, m, n, d1, d2, B, &startNode);

// find all Shortest Common Supersequence
	vector<string> result;
	BacktrackSCS(startNode, result, scsLen);
	return result;
}

int SCS2LenFast(const string& X, const string& Y, const int originalLen) {
	int m = X.length(), n = Y.length();
	int d1 = originalLen - m, d2 = originalLen - n;
	vector<vector<int>> L;

// fill lookup table
	int result = SCSLengthFast2(X, Y, m, n, d1, d2, L);
	return result;
}

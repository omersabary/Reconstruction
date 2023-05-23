#include "SCSN.hpp"
#include <map>
#include <set>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "Strings.hpp"
using namespace std;

// map of letters and a vector of the string indexes with that letter in current position
map<char, vector<int> > LettersMap(const vector<string>& strings, const vector<int>& index) {
	map<char, vector<int> > lettersMap;
	for (int i = 0; i < (int) strings.size(); i++) {
		if (index[i] == 0) {
			lettersMap['$'].push_back(i);
		}
		else {
			lettersMap[strings[i][index[i] - 1]].push_back(i);
		}
	}
	return lettersMap;
}

// fill L array.
int SCSNLen(const vector<string>& strings, vector<int>& L, const IndexVector& firstIndex) {
	int n = strings.size();
	for (IndexVector currentIndex = firstIndex; not currentIndex.IsPastEnd(); currentIndex++) {
		map<char, vector<int> > lettersMap = LettersMap(strings, currentIndex.Vector());
		if ((int) lettersMap['$'].size() >= n - 1) { // at most 1 non zero index.
			for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
				if (it->first != '$') {
					int nonZeroDimension = it->second[0];
					vector<int> indexes = currentIndex.Vector();
					L[currentIndex.AbsInd()] = indexes[nonZeroDimension];
					break;
				}
			}
		}
		else { // at least 2 non zero indexes.
			lettersMap.erase('$');
			int minLen = INT_MAX, currentValue;
			for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
				vector<int> K = it->second;
				IndexVector kIndex = currentIndex.Beta(K);
				currentValue = L[kIndex.AbsInd()];
				if (currentValue < minLen) {
					minLen = currentValue;
				}
			}
			assert(minLen!=INT_MAX);
			L[currentIndex.AbsInd()] = minLen + 1;
		}
	}
	return L[firstIndex.ArraySize() - 1];
}

vector<string> BacktrackSCS(const vector<string>& strings, IndexVector index, const vector<int>& L) {
	int n = strings.size();
	map<char, vector<int> > lettersMap = LettersMap(strings, index.Vector());
	if ((int) lettersMap['$'].size() >= n - 1) { // at most 1 non zero index.
		vector<string> result(1);
		// find non zero dimension
		for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
			if (it->first != '$') {
				int nonZeroDimension = it->second[0];
				vector<int> indexes = index.Vector();
				result[0] = strings[nonZeroDimension].substr(0, indexes[nonZeroDimension]);
				break;
			}
		}
		return result;
	}
	else { // at least 2 non zero indexes.
		lettersMap.erase('$');
		//int minLen = INT_MAX, currentValue;
		int minValue = L[index.AbsInd()] - 1, currentValue;
		vector<char> lettersToAdd;
		vector<IndexVector> minIndexes;
		for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
			vector<int> K = it->second;
			IndexVector currentIndex = index.Beta(K);
			currentValue = L[currentIndex.AbsInd()];
			if (currentValue == minValue) {
				lettersToAdd.push_back(it->first);
				minIndexes.push_back(currentIndex);
			}
		}
		vector<string> result;
		for (unsigned int i = 0; i < minIndexes.size(); i++) {
			vector<string> currentResult = BacktrackSCS(strings, minIndexes[i], L);
			for (vector<string>::iterator it = currentResult.begin(); it != currentResult.end(); it++) {
				it->push_back(lettersToAdd[i]);
			}
			result.insert(result.end(), currentResult.begin(), currentResult.end());
		}
		return result;
	}
}

// return all SCSs of n strings
vector<string> SCSN(const vector<string>& strings) {
	int n = strings.size();
	vector<int> dimSizes(n);
	for (int i = 0; i < n; i++) {
		dimSizes[i] = strings[i].size() + 1;
	}

	IndexVector indexVector(dimSizes);
	vector<int> L(indexVector.ArraySize());
	SCSNLen(strings, L, indexVector);

	vector<string> result = BacktrackSCS(strings, indexVector.Last(), L);
	return result;
}

// return length of SCS of n strings
int SCSNLen(const vector<string>& strings) {
	int n = strings.size();
	vector<int> dimSizes(n);
	for (int i = 0; i < n; i++) {
		dimSizes[i] = strings[i].size() + 1;
	}

	IndexVector indexVector(dimSizes);
	vector<int> L(indexVector.ArraySize());

	return SCSNLen(strings, L, indexVector);
}

vector<string> BacktrackSCSFast(const vector<string>& strings, IndexVectorD index, const vector<int>& L) {
	int n = strings.size();
	map<char, vector<int> > lettersMap = LettersMap(strings, index.Vector());
	if ((int) lettersMap['$'].size() >= n - 1) { // at most 1 non zero index.
		string currentStr;
		// find non zero dimension
		for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
			if (it->first != '$') {
				int nonZeroDimension = it->second[0];
				vector<int> indexes = index.Vector();
				currentStr = strings[nonZeroDimension].substr(0, indexes[nonZeroDimension]);
				break;
			}
		}
		return vector<string>( { currentStr });
	}
	else { // at least 2 non zero indexes.
		lettersMap.erase('$');
		//int minLen = INT_MAX, currentValue;
		int minValue = L[index.AbsInd()] - 1, currentValue;
		vector<char> lettersToAdd;
		vector<IndexVectorD> minIndexes;
		for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
			vector<int> K = it->second;
			currentValue = L[index.AbsIndBeta(K)];
			if (currentValue == minValue) {
				lettersToAdd.push_back(it->first);
				minIndexes.push_back(index.BetaClass(K));
			}
		}
		vector<string> result;
		for (unsigned int i = 0; i < minIndexes.size(); i++) {
			vector<string> currentResult = BacktrackSCSFast(strings, minIndexes[i], L);
			for (vector<string>::iterator it = currentResult.begin(); it != currentResult.end(); it++) {
				it->push_back(lettersToAdd[i]);
			}
			result.insert(result.end(), currentResult.begin(), currentResult.end());
		}
		return result;
	}
}

int SCSNLenFast(const vector<string>& strings, vector<int>& L, const IndexVectorD& firstIndex) {
	int n = strings.size();
	L[0] = 0;
	for (IndexVectorD currentIndex = firstIndex; not currentIndex.IsEnd(); currentIndex++) {
		map<char, vector<int> > lettersMap = LettersMap(strings, currentIndex.Vector());
		if ((int) lettersMap['$'].size() >= n - 1) { // at most 1 non zero index.
			for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
				if (it->first != '$') {
					int nonZeroDimension = it->second[0];
					vector<int> indexes = currentIndex.Vector();
					L[currentIndex.AbsInd()] = indexes[nonZeroDimension];
					break;
				}
			}
		}
		else { // at least 2 non zero indexes.
			lettersMap.erase('$');
			int minLen = INT_MAX, currentValue;
			for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
				vector<int> K = it->second;
				currentValue = L[currentIndex.AbsIndBeta(K)];
				if (currentValue < minLen) {
					minLen = currentValue;
				}
			}

			// if minLen==INT_MAX at this stage there is no legit path to this vector, so leave it at INT_MAX
			if (minLen < INT_MAX) {
				L[currentIndex.AbsInd()] = minLen + 1;
			}
		}
	}
	assert(L[firstIndex.MaxMemorySize() - 1] != INT_MAX); // There is always a path
	return L[firstIndex.MaxMemorySize() - 1];
}

struct intVecHash {
	size_t operator()(const vector<int>& vec) const {
		size_t seed = vec.size();
		for (auto& i : vec) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

int SCSNLenFastLessMem(const vector<string>& strings, unordered_map<int, int>& L, const IndexVectorD& firstIndex) {
	int n = strings.size();
	L[0] = 0;
	for (IndexVectorD currentIndex = firstIndex; not currentIndex.IsEnd(); currentIndex++) {
		map<char, vector<int> > lettersMap = LettersMap(strings, currentIndex.Vector());
		if ((int) lettersMap['$'].size() >= n - 1) { // at most 1 non zero index.
			for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
				if (it->first != '$') {
					int nonZeroDimension = it->second[0];
					vector<int> indexes = currentIndex.Vector();
					L[currentIndex.AbsInd()] = indexes[nonZeroDimension];
					break;
				}
			}
		}
		else { // at least 2 non zero indexes.
			lettersMap.erase('$');
			int minLen = INT_MAX;
			for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
				vector<int> K = it->second;
				unordered_map<int, int>::iterator valueIt = L.find(currentIndex.AbsIndBeta(K));
				if (valueIt == L.end()) {
					continue;
				}
				if (valueIt->second < minLen) {
					minLen = valueIt->second;
				}
			}

			// if minLen==INT_MAX at this stage there is no legit path to this vector, so leave it at INT_MAX
			if (minLen < INT_MAX) {
				L[currentIndex.AbsInd()] = minLen + 1;
			}
			else {
				L[currentIndex.AbsInd()] = INT_MAX;
			}
		}
	}
	assert(L[firstIndex.MaxMemorySize() - 1] != INT_MAX); // There is always a path
	return L[firstIndex.MaxMemorySize() - 1];
}

typedef pair<int, int> mypair;
bool comparator(const mypair& l, const mypair& r) {
	return l.first < r.first;
}

int EstimateMem(const vector<int>& dimSizes, const vector<int>& D) {
	int dimNum = dimSizes.size();
	vector<mypair> DN(dimNum);
	for (int i = 0; i < dimNum; i++) {
		DN[i] = mypair(D[i], dimSizes[i]);
	}
	sort(DN.begin(), DN.end(), comparator);
	vector<int> upperLimits(dimNum), maxRanges(dimNum);
	upperLimits[0] = DN[0].second - 1;
	for (int i = 1; i < dimNum; i++) {
		upperLimits[i] = DN[i].first + DN[i].second - 1;
	}
	maxRanges[0] = *min_element(upperLimits.begin(), upperLimits.end()) + 1;
	// go over other dims. from one iteration to the next only limits for current dim and previous have to be changed
	for (int dim = 1; dim < dimNum; dim++) {
		upperLimits[dim] = DN[dim].second - 1;
		upperLimits[dim - 1] = DN[dim - 1].first;
		maxRanges[dim] = *min_element(upperLimits.begin(), upperLimits.end()) + DN[dim].first + 1;
	}
	int estimateMem = 1;
	for (int dim = 0; dim < dimNum; dim++) {
		estimateMem *= maxRanges[dim];
	}
	//cout << "Estimated memory:\t" << estimateMem << endl;
	return estimateMem;
}
// Does L[0] require fixing
vector<string> SCSNFast(const vector<string>& strings, const int originalLen) {
	int n = strings.size();
	vector<int> dimSizes(n);
	vector<int> D(n);
	for (int i = 0; i < n; i++) {
		dimSizes[i] = strings[i].size() + 1;
		D[i] = originalLen - strings[i].size();
	}
	IndexVectorD indexVector(dimSizes, D);
	vector<int> L(indexVector.MaxMemorySize(), INT_MAX);
	SCSNLenFast(strings, L, indexVector);
	return BacktrackSCSFast(strings, indexVector.Last(), L);
}

int SCSNLenFast(const vector<string>& strings, const int originalLen) {
	int n = strings.size();
	vector<int> dimSizes(n);
	vector<int> D(n);
	for (int i = 0; i < n; i++) {
		dimSizes[i] = strings[i].size() + 1;
		D[i] = originalLen - strings[i].size();
	}
	IndexVectorD indexVector(dimSizes, D);
	vector<int> L(indexVector.MaxMemorySize(), INT_MAX);
	return SCSNLenFast(strings, L, indexVector);
}

// SLOWER THAN SCSNLenFast By 20%
int SCSNLenFastLessMem(const vector<string>& strings, const int originalLen) {
	int n = strings.size();
	vector<int> dimSizes(n);
	vector<int> D(n);
	for (int i = 0; i < n; i++) {
		dimSizes[i] = strings[i].size() + 1;
		D[i] = originalLen - strings[i].size();
	}
	IndexVectorD indexVector(dimSizes, D);
	unordered_map<int, int> L;
	L.reserve(EstimateMem(dimSizes, D));
	return SCSNLenFastLessMem(strings, L, indexVector);
}

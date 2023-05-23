#include <iostream>
#include <algorithm>
#include <cassert>
#include <set>
#include "Strings.hpp"
#include "SCSN.hpp"
#include "SCS2.hpp"
#include "NChooseK.hpp"
using namespace std;

ostream& operator<<(ostream& os, const vector<string>& a) {
	for (unsigned i = 0; i < a.size(); i++) {
		os << "str" << i << ":\t" << a[i] << endl;
	}
	return os;
}

ostream& operator<<(ostream& os, const set<string>& a) {
	int i = 0;
	for (set<string>::iterator it = a.begin(); it != a.end(); it++) {
		os << "str" << i++ << ":\t" << *it << endl;
	}
	return os;
}

ostream& operator<<(ostream& os, const vector<int>& a) {
	for (unsigned i = 0; i < a.size(); i++) {
		os << a[i] << "\t";
	}
	os << endl;
	return os;
}
string MakeRandomString4(const int len, mt19937& generator) {
	string strand;
	vector<char> letters = { 'G', 'T', 'A', 'C' };
	uniform_int_distribution<int> distribution(0, 3);
	int rnd;
	for (int i = 0; i < len; i++) {
		rnd = distribution(generator);
		strand += letters[rnd];
	}
	return strand;
}

string MakeRandomString2(const int len, mt19937& generator) {
	string strand;
	vector<char> letters = { '0', '1' };
	uniform_int_distribution<int> distribution(0, 1);
	int rnd;
	for (int i = 0; i < len; i++) {
		rnd = distribution(generator);
		strand += letters[rnd];
	}
	return strand;
}

string CopyStrandWithDeletion(const string& source, mt19937& generator, double delProb) {
	uniform_real_distribution<double> distribution(0.0, 1.0);
	string result;
	for (size_t i = 0; i < source.size(); i++) {
		if (distribution(generator) > delProb) {
			result.push_back(source[i]);
		}
	}
	return result;
}

void RandomStrings(vector<string>& strings, const int stringNum, const int stringLen, mt19937& generator) {
	strings = vector<string>(stringNum);
	for (int i = 0; i < stringNum; i++) {
		strings[i] = MakeRandomString4(stringLen, generator);
	}
}

void RandomStringsRandomLen(vector<string>& strings, const int stringNum, const int minStrLen, const int maxStrLen,
		mt19937& generator) {
	strings = vector<string>(stringNum);
	uniform_int_distribution<int> strLenDistribution(minStrLen, maxStrLen);
	for (int i = 0; i < stringNum; i++) {
		strings[i] = MakeRandomString4(strLenDistribution(generator), generator);
	}
}

string RandomCopies2(vector<string>& copies, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	string original = MakeRandomString2(originalLen, generator);
	copies = vector<string>(copiesNum);
	for (int i = 0; i < copiesNum; i++) {
		copies[i] = CopyStrandWithDeletion(original, generator, delProb);
	}
	return original;
}

string RandomCopies4(vector<string>& copies, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	string original = MakeRandomString4(originalLen, generator);
	copies = vector<string>(copiesNum);
	for (int i = 0; i < copiesNum; i++) {
		copies[i] = CopyStrandWithDeletion(original, generator, delProb);
	}
	return original;
}

// A function to find the number of times string b occurs in string a as a subsequence

int CountEmbeddings(const string& a, const string& b) {
	int m = a.length();
	int n = b.length();

	// Create a table to store results of sub-problems
	int lookup[m + 1][n + 1] = { { 0 } };

	// If first string is empty
	for (int i = 0; i <= n; ++i)
		lookup[0][i] = 0;

	// If second string is empty
	for (int i = 0; i <= m; ++i)
		lookup[i][0] = 1;

	// Fill lookup[][] in bottom up manner
	for (int i = 1; i <= m; i++) {
		for (int j = 1; j <= n; j++) {
			// If last characters are same, we have two options:
			// 1. consider last characters of both strings in solution
			// 2. ignore last character of first string
			if (a[i - 1] == b[j - 1])
				lookup[i][j] = lookup[i - 1][j - 1] + lookup[i - 1][j];

			else
				// If last character are different, ignore last character of first string
				lookup[i][j] = lookup[i - 1][j];
		}
	}

	return lookup[m][n];
}

int min(const int a, const int b, const int c) {
	int min;
	min = (a < b) ? a : b;
	min = (c < min) ? c : min;
	return min;
}

int EditDistance(const string& s, const string& t) {
	int m = s.size();
	int n = t.size();
	int deletion, insertion, substitution;
	// For all i and j, d[i,j] will hold the Levenshtein distance between the first i characters of s and the first j
	// characters of t.
	// Note that d has (m+1) x (n+1) values.
	//let d be a 2-d array of int with dimensions [0..m, 0..n]
	int d[m + 1][n + 1];
	for (int i = 0; i <= m; i++) {
		// the distance of any first string to an empty second string
		// (transforming the string of the first i characters of s into the empty string requires i deletions)
		d[i][0] = i;
	}
	for (int j = 0; j <= n; j++) {
		d[0][j] = j;	// the distance of any second string to an empty first string
	}
	for (int j = 1; j <= n; j++) {
		for (int i = 1; i <= m; i++) {
			if (s[i - 1] == t[j - 1]) {	// no operation required
				d[i][j] = d[i - 1][j - 1];
			}
			else {	//minimum of a deletion, an insertion, a substitution.
				deletion = d[i - 1][j] + 1;
				insertion = d[i][j - 1] + 1;
				substitution = d[i - 1][j - 1] + 1;
				d[i][j] = min(deletion, insertion, substitution);
			}
		}
	}
	return d[m][n];
}

int EditDistanceDI(const string& s, const string& t) {
	int lcslen = LCS2Len(s, t);
	int m = s.size();
	int n = t.size();
	return m + n - 2 * lcslen;
}

vector<long> SCSScores(const vector<string>& inputStrings, const vector<string>& SCSs) {
	int SCSNum = SCSs.size();
	int inputStringNum = inputStrings.size();
	vector<long> scsScores(SCSNum, 1);
	for (int i = 0; i < SCSNum; i++) {
		for (int inputIndex = 0; inputIndex < inputStringNum; inputIndex++) {
			scsScores[i] *= CountEmbeddings(SCSs[i], inputStrings[inputIndex]);
		}
	}
	return scsScores;
}

vector<string> AllMostLikely(const vector<string>& inputStrings, const vector<string>& SCSs) {
	int SCSNum = SCSs.size();
	int inputStringNum = inputStrings.size();
	vector<long> scsScores(SCSNum, 1);
	for (int i = 0; i < SCSNum; i++) {
		for (int inputIndex = 0; inputIndex < inputStringNum; inputIndex++) {
			scsScores[i] *= CountEmbeddings(SCSs[i], inputStrings[inputIndex]);
		}
	}
	long maxScore = *max_element(scsScores.begin(), scsScores.end());
	vector<string> result;
	for (int i = 0; i < SCSNum; i++) {
		if (scsScores[i] == maxScore) {
			result.push_back(SCSs[i]);
		}
	}
	return result;
}

string FirstMostLikely(const vector<string>& inputStrings, const vector<string>& SCSs) {
	int SCSNum = SCSs.size();
	int inputStringNum = inputStrings.size();
	vector<long> scsScores(SCSNum, 1);
	for (int i = 0; i < SCSNum; i++) {
		for (int inputIndex = 0; inputIndex < inputStringNum; inputIndex++) {
			scsScores[i] *= CountEmbeddings(SCSs[i], inputStrings[inputIndex]);
		}
	}
	int maxScoreIndex = max_element(scsScores.begin(), scsScores.end()) - scsScores.begin();
	return SCSs[maxScoreIndex];
}

string MinEditDistanceSum(const vector<string>& inputStrings, const vector<string>& SCSs) {
	int SCSNum = SCSs.size();
	int inputStringNum = inputStrings.size();
	vector<long> scsScores(SCSNum);
	for (int i = 0; i < SCSNum; i++) {
		for (int inputIndex = 0; inputIndex < inputStringNum; inputIndex++) {
			scsScores[i] += EditDistanceDI(SCSs[i], inputStrings[inputIndex]);
		}
	}
	int minScoreIndex = min_element(scsScores.begin(), scsScores.end()) - scsScores.begin();
	return SCSs[minScoreIndex];
}

string FirstMostLikelyAndCount(const vector<string>& inputStrings, const vector<string>& SCSs, int& mostLikelyNum) {
	int SCSNum = SCSs.size();
	int inputStringNum = inputStrings.size();
	vector<long> scsScores(SCSNum, 1);
	for (int i = 0; i < SCSNum; i++) {
		for (int inputIndex = 0; inputIndex < inputStringNum; inputIndex++) {
			scsScores[i] *= CountEmbeddings(SCSs[i], inputStrings[inputIndex]);
		}
	}
	int maxScoreIndex = max_element(scsScores.begin(), scsScores.end()) - scsScores.begin();
	long maxScore = scsScores[maxScoreIndex];
	int countMostLikely = 0;
	for (int i = 0; i < SCSNum; i++) {
		if (scsScores[i] == maxScore) {
			countMostLikely++;
		}
	}
	mostLikelyNum = countMostLikely;
	return SCSs[maxScoreIndex];
}

vector<string> AllMostLikely(const string& X, const string& Y, const vector<string>& SCSs) {
	int SCSNum = SCSs.size();
	vector<long> scsScores(SCSNum, 1);
	for (int i = 0; i < SCSNum; i++) {
		scsScores[i] = CountEmbeddings(SCSs[i], X) * CountEmbeddings(SCSs[i], Y);
	}
	long maxScore = *max_element(scsScores.begin(), scsScores.end());
	vector<string> result;
	for (int i = 0; i < SCSNum; i++) {
		if (scsScores[i] == maxScore) {
			result.push_back(SCSs[i]);
		}
	}
	return result;
}

// Find all SCSs of inputStrings. Return most likely SCS.
string MostLikelySCS(const vector<string>& inputStrings, const long double p, const int originalLen,
		mt19937& generator) {
	vector<string> SCS = SCSNFast(inputStrings, originalLen);
	int SCSNum = SCS.size();
	assert(SCSNum != 0);
	if (SCSNum == 1) {
		return SCS[0];
	}
	return FirstMostLikely(inputStrings, SCS);;
}

//string MostLikelyKCS(const vector<string>& inputStrings, const long double p, const int originalLen,
//		mt19937& generator) {
//	vector<string> SCS = SCSNFast(inputStrings, originalLen);
//	int SCSNum = SCS.size();
//	assert(SCSNum != 0);
//	if ((int) SCS[0].size() < originalLen) {
//		SCS = KCSDP(inputStrings[0], inputStrings[1], originalLen);
//		SCSNum = SCS.size();
//		assert(SCSNum != 0);
//	}
//	if (SCSNum == 1) {
//		return SCS[0];
//	}
//	return FirstMostLikely(inputStrings, SCS);;
//}

bool IsSupersequenceOfAll(const string& candidate, const vector<string>& inputStrings) {
	for (unsigned i = 0; i < inputStrings.size(); i++) {
		if (not IsSubSequence(inputStrings[i], candidate)) {
			return false;
		}
	}
	return true;
}

//	compute scs len for all k tuples
//	find longest scs len -> maxLen
//	find all SCSs of a k tuple with max scs len
//	filter scss which are not a supersequence of any of the input strings. if filtered is empty, return first SCS.
//	return most likely among remaining scss

string MostLikelyFilteredSCS(const vector<string>& inputStrings, const int k, const long double p,
		const int originalLen) {
	vector<vector<string> > ktuples = StringKtuples(inputStrings, k);
	int ktuplesNum = ktuples.size();
	vector<int> scsLen(ktuplesNum);
	for (int i = 0; i < ktuplesNum; i++) {
		scsLen[i] = SCSNLenFast(ktuples[i], originalLen);
	}
	// find longest scs (first if many).
	int maxLenSCSIndex = max_element(scsLen.begin(), scsLen.end()) - scsLen.begin();
	vector<string> SCSs = SCSNFast(ktuples[maxLenSCSIndex], originalLen);

	// filter SCSs which are not supersequences of all input strings
	vector<string> filteredSCSs;
	for (unsigned i = 0; i < SCSs.size(); i++) {
		if (IsSupersequenceOfAll(SCSs[i], inputStrings)) {
			filteredSCSs.push_back(SCSs[i]);
		}
	}

	// return most likely among filtered
	int SCSNum = filteredSCSs.size();
	if (SCSNum == 0) { // all SCSs were filtered
		return SCSs[0];
	}
	if (SCSNum == 1) {
		return filteredSCSs[0];
	}
	return FirstMostLikely(inputStrings, filteredSCSs);;
}

//	compute scs len for all k tuples
//	find longest scs len -> maxLen
//	find all SCSs of a k tuple with max scs len
//	filter scss which are not a supersequence of any of the input strings.
//  if filtered is empty, return scs with min sum of Levinstein distance from all copies.
//	Otherwise, return most likely among remaining scss. first, if many.

string MostLikelyMaxLenSCS(const vector<string>& inputStrings, const int k, const long double p, const int originalLen,
		bool& foundCorrectSize, int& foundKtupleIndex) {
	vector<vector<string> > ktuples = SortedStringKtuplesBySumLen(inputStrings, k);
	int ktuplesNum = ktuples.size();
	int maxLenSCSIndex;
	bool foundCorrect = false;
	foundKtupleIndex = -1; // index for not found
	vector<int> scsLen(ktuplesNum);
	// try to find ktuple whose scs len is correct
	for (int i = 0; i < ktuplesNum; i++) {
		scsLen[i] = SCSNLenFast(ktuples[i], originalLen);
		if (scsLen[i] == originalLen) {
			foundCorrect = true;
			foundKtupleIndex = i;
			maxLenSCSIndex = i;
			break;
		}
	}
	foundCorrectSize = foundCorrect;
	// find longest scs (first if many).
	if (not foundCorrect) {
		maxLenSCSIndex = max_element(scsLen.begin(), scsLen.end()) - scsLen.begin();
	}
	vector<string> SCSs = SCSNFast(ktuples[maxLenSCSIndex], originalLen);

	// filter SCSs which are not supersequences of all input strings
	vector<string> filteredSCSs;
	for (unsigned i = 0; i < SCSs.size(); i++) {
		if (IsSupersequenceOfAll(SCSs[i], inputStrings)) {
			filteredSCSs.push_back(SCSs[i]);
		}
	}

	// return most likely among filtered
	int SCSNum = filteredSCSs.size();
	if (SCSNum == 0) { // all SCSs were filtered. return scs with min sum of Levinstein distance.
		return MinEditDistanceSum(inputStrings, SCSs);
	}
	if (SCSNum == 1) {
		return filteredSCSs[0];
	}
	return FirstMostLikely(inputStrings, filteredSCSs);;
}

string MostLikelyMaxLenSCS(const vector<string>& inputStrings, const int k, const long double p,
		const int originalLen) {
	vector<vector<string> > ktuples = StringKtuples(inputStrings, k);
	int ktuplesNum = ktuples.size();
	int maxLenSCSIndex;
	bool foundCorrect = false;
	vector<int> scsLen(ktuplesNum);
	// try to find ktuple whose scs len is correct
	for (int i = 0; i < ktuplesNum; i++) {
		scsLen[i] = SCSNLenFast(ktuples[i], originalLen);
		if (scsLen[i] == originalLen) {
			foundCorrect = true;
			maxLenSCSIndex = i;
			break;
		}
	}
	// find longest scs (first if many).
	if (not foundCorrect) {
		maxLenSCSIndex = max_element(scsLen.begin(), scsLen.end()) - scsLen.begin();
	}
	vector<string> SCSs = SCSNFast(ktuples[maxLenSCSIndex], originalLen);

	// filter SCSs which are not supersequences of all input strings
	vector<string> filteredSCSs;
	for (unsigned i = 0; i < SCSs.size(); i++) {
		if (IsSupersequenceOfAll(SCSs[i], inputStrings)) {
			filteredSCSs.push_back(SCSs[i]);
		}
	}

	// return most likely among filtered
	int SCSNum = filteredSCSs.size();
	if (SCSNum == 0) { // all SCSs were filtered. return scs with min sum of Levinstein distance.
		return MinEditDistanceSum(inputStrings, SCSs);
	}
	if (SCSNum == 1) {
		return filteredSCSs[0];
	}
	return FirstMostLikely(inputStrings, filteredSCSs);;
}

string MostLikelyCorrectSizeSCS(const vector<string>& inputStrings, const int k, const long double p,
		const int originalLen) {
	vector<vector<string> > ktuples = SortedStringKtuplesBySumLen(inputStrings, k);
	int ktuplesNum = ktuples.size();
	bool foundCorrect = false;
	int correctLenSCSIndex;
	vector<int> scsLen(ktuplesNum);
	for (int i = 0; i < ktuplesNum; i++) {
		scsLen[i] = SCSNLenFast(ktuples[i], originalLen);
		if (scsLen[i] == originalLen) {
			foundCorrect = true;
			correctLenSCSIndex = i;
			break;
		}
	}

	if (not foundCorrect) { // no ktuple with correct size scs
		return "";
	}

	vector<string> SCSs = SCSNFast(ktuples[correctLenSCSIndex], originalLen);

	// filter SCSs which are not supersequences of all input strings
	vector<string> filteredSCSs;
	for (unsigned i = 0; i < SCSs.size(); i++) {
		if (IsSupersequenceOfAll(SCSs[i], inputStrings)) {
			filteredSCSs.push_back(SCSs[i]);
		}
	}

	// return most likely among filtered
	int SCSNum = filteredSCSs.size();
	assert(SCSNum > 0);
	if (SCSNum == 1) {
		return filteredSCSs[0];
	}
	return FirstMostLikely(inputStrings, filteredSCSs);;
}

//
int CorrectSizeKtupleNumSortByLetterMax(const vector<string>& inputStrings, const int k, const int originalLen) {
	vector<vector<string> > ktuples = SortedStringKtuplesBySumLetterMax(inputStrings, k);
	int ktuplesNum = ktuples.size();
	for (int i = 0; i < ktuplesNum; i++) {
		int tempLen = SCSNLenFast(ktuples[i], originalLen);
		if (tempLen == originalLen) {
			return i;
		}
	}
	return -1; // not found
}

int CorrectSizeKtupleNumSortBySumLen(const vector<string>& inputStrings, const int k, const int originalLen) {
	vector<vector<string> > ktuples = SortedStringKtuplesBySumLen(inputStrings, k);
	int ktuplesNum = ktuples.size();
	for (int i = 0; i < ktuplesNum; i++) {
		int tempLen = SCSNLenFast(ktuples[i], originalLen);
		if (tempLen == originalLen) {
			return i;
		}
	}
	return -1; // not found
}

//	compute scs len for all k tuples
//	find longest scs len -> maxLen
//	if maxLen < originalLen return empty string. Otherwise:
//	find all SCSs of a ktuple with correct scs len
//	filter scss which are not a supersequence of any of the input strings
//	return most likely among remaining scss. first, if more than one.

string MostLikelyCorrectSizeSCS(const vector<string>& inputStrings, const int k, const long double p,
		const int originalLen, int& foundKtupleIndex) {
	vector<vector<string> > ktuples = SortedStringKtuplesBySumLen(inputStrings, k);
	int ktuplesNum = ktuples.size();
	bool foundCorrect = false;
	int correctLenSCSIndex;
	foundKtupleIndex = -1; //-1 for not found
	vector<int> scsLen(ktuplesNum);
	for (int i = 0; i < ktuplesNum; i++) {
		scsLen[i] = SCSNLenFast(ktuples[i], originalLen);
		if (scsLen[i] == originalLen) {
			foundCorrect = true;
			foundKtupleIndex = i;
			correctLenSCSIndex = i;
			break;
		}
	}

	if (not foundCorrect) { // no ktuple with correct size scs
		return "";
	}

	vector<string> SCSs = SCSNFast(ktuples[correctLenSCSIndex], originalLen);

	// filter SCSs which are not supersequences of all input strings
	vector<string> filteredSCSs;
	for (unsigned i = 0; i < SCSs.size(); i++) {
		if (IsSupersequenceOfAll(SCSs[i], inputStrings)) {
			filteredSCSs.push_back(SCSs[i]);
		}
	}

	// return most likely among filtered
	int SCSNum = filteredSCSs.size();
	assert(SCSNum > 0);
	if (SCSNum == 1) {
		return filteredSCSs[0];
	}
	return FirstMostLikely(inputStrings, filteredSCSs);;
}

bool HasPairWithCorrectSizeSCS(const vector<string>& inputStrings, const int originalLen) {
	vector<vector<string> > ktuples = StringKtuples(inputStrings, 2);
	int ktuplesNum = ktuples.size();
	int currentLen;
	for (int i = 0; i < ktuplesNum; i++) {
		currentLen = SCSNLenFast(ktuples[i], originalLen);
		if (currentLen == originalLen) {
			return true;
		}
	}
	return false;
}

string MinEditDistanceFilteredSCS(const vector<string>& inputStrings, const int k, const long double p,
		const int originalLen) {
	vector<vector<string> > ktuples = StringKtuples(inputStrings, k);
	int ktuplesNum = ktuples.size();
	vector<int> scsLen(ktuplesNum);
	for (int i = 0; i < ktuplesNum; i++) {
		scsLen[i] = SCSNLenFast(ktuples[i], originalLen);
	}
	// find longest scs (first if many).
	int maxLenSCSIndex = max_element(scsLen.begin(), scsLen.end()) - scsLen.begin();
	vector<string> SCSs = SCSNFast(ktuples[maxLenSCSIndex], originalLen);

	// filter SCSs which are not supersequences of all input strings
	vector<string> filteredSCSs;
	for (unsigned i = 0; i < SCSs.size(); i++) {
		if (IsSupersequenceOfAll(SCSs[i], inputStrings)) {
			filteredSCSs.push_back(SCSs[i]);
		}
	}

	// return most likely among filtered
	int SCSNum = filteredSCSs.size();
	if (SCSNum == 0) { // all SCSs were filtered
		return SCSs[0];
	}
	if (SCSNum == 1) {
		return filteredSCSs[0];
	}
	return MinEditDistanceSum(inputStrings, filteredSCSs);;
}

// k=2 for now
//string MostLikelyFilteredKCS(const vector<string>& inputStrings, const long double p, const int originalLen) {
//	vector<vector<string> > ktuples = StringKtuples(inputStrings, 2);
//	int ktuplesNum = ktuples.size();
//	vector<int> scsLen(ktuplesNum);
//	for (int i = 0; i < ktuplesNum; i++) {
//		scsLen[i] = SCSNLenFast(ktuples[i], originalLen);
//	}
//	// find longest scs (first if many).
//	int maxLenSCSIndex = max_element(scsLen.begin(), scsLen.end()) - scsLen.begin();
//	int maxSCSLen = scsLen[maxLenSCSIndex];
//	vector<string> SCSs;
//	if (maxSCSLen == originalLen) {
//		SCSs = SCSNFast(ktuples[maxLenSCSIndex], originalLen);
//	}
//	else {
//		SCSs = KCSDP(ktuples[maxLenSCSIndex][0], ktuples[maxLenSCSIndex][1], originalLen);
//	}
//
//	// filter SCSs which are not supersequences of all input strings
//	vector<string> filteredSCSs;
//	for (unsigned i = 0; i < SCSs.size(); i++) {
//		if (IsSupersequenceOfAll(SCSs[i], inputStrings)) {
//			filteredSCSs.push_back(SCSs[i]);
//		}
//	}
//
//	// return most likely among filtered
//	int SCSNum = filteredSCSs.size();
//	if (SCSNum == 0) { // all SCSs were filtered
//		return SCSs[0];
//	}
//	if (SCSNum == 1) {
//		return filteredSCSs[0];
//	}
//
//	return FirstMostLikely(inputStrings, filteredSCSs);
//
//}

void MostLikelyEitan(const vector<string>& inputStrings, const int originalLen, int& scsNum, int& mostLikelyNum) {
	vector<string> SCSs = SCS2Fastest(inputStrings[0], inputStrings[1], originalLen);
	scsNum = SCSs.size();
	FirstMostLikelyAndCount(inputStrings, SCSs, mostLikelyNum);
}

string MostLikelyEitanFull(const vector<string>& inputStrings, const int originalLen, int& scsNum, int& mostLikelyNum) {
	vector<string> SCSs = SCSNFast(inputStrings, originalLen);
	scsNum = SCSs.size();
	return FirstMostLikelyAndCount(inputStrings, SCSs, mostLikelyNum);
}

vector<string> AllMostLikelyEitan(const vector<string>& inputStrings, const int originalLen) {
	vector<string> SCSs = SCSNFast(inputStrings, originalLen);
	return AllMostLikely(inputStrings, SCSs);
}

//void MostLikelyNewAlgorithmVsOld(const vector<string>& inputStrings) {
//	vector<string> SCSs = SCS2(inputStrings[0], inputStrings[1]);
//	vector<string> oldSCSs = AllMostLikely(inputStrings, SCSs);
//	set<string> oldSet(oldSCSs.begin(), oldSCSs.end());
//
//	vector<string> newSCSs = MostLikelySCS2(inputStrings[0], inputStrings[1]);
//	set<string> newSet(newSCSs.begin(), newSCSs.end());
//	set<string> intersection;
//	set_intersection(oldSet.begin(), oldSet.end(), newSet.begin(), newSet.end(),
//			inserter(intersection, intersection.begin()));
//	if (oldSet == newSet) {
//		//cout<<"Sets Equal"<<endl;
//	}
//	else {
//		if (intersection.empty()) {
//			cout << "Sets Not Equal" << endl;
//			cout << "Old Set: " << endl << oldSet;
//			cout << "New Set: " << endl << newSet;
//			cout << "Intersection Empty!" << endl;
//		}
//	}
//}

// for SCS index from 0 to originalLen-1: max frequency of {A,G,T,C} up to index-delj for each string j.
vector<map<char, int>> MaxLetterFrequency(const vector<string>& inputStrings, const int originalLen) {
	int stringNum = inputStrings.size();
	vector<map<char, int>> maxFreq(originalLen);
	map<char, int> currentMax;
	vector<map<char, int>> frequencies(stringNum);
	vector<int> lastIndex(stringNum, -1);
	for (int index = 0; index < originalLen; index++) {
		//int halfIndex = index / 2;
		//int maxHalfIndex = max(5, halfIndex);
		for (int j = 0; j < stringNum; j++) {
			int delj = originalLen - (int) inputStrings[j].size();
			//delj = min(delj, maxHalfIndex);
			if (index - delj >= 0) {
				if (lastIndex[j] != index - delj) {
					lastIndex[j] = index - delj;
					char currentLetter = inputStrings[j][index - delj];
					frequencies[j][currentLetter]++;
					if (frequencies[j][currentLetter] > currentMax[currentLetter]) {
						currentMax[currentLetter] = frequencies[j][currentLetter];
					}
				}
			}

		}
		maxFreq[index] = currentMax;
	}
	return maxFreq;
}

// Returns true if str1 is a subsequence of str2. m is length of str1 and n is length of str2

bool IsSubSequence(const string& str1, const string& str2) {
	int j = 0; // For index of str1 (or subsequence
	int m = str1.size();
	int n = str2.size();
// Traverse str2 and str1, and compare current character of str2 with first unmatched char of str1,
// if matched then move ahead in str1
	for (int i = 0; i < n && j < m; i++)
		if (str1[j] == str2[i])
			j++;

// If all characters of str1 were found in str2
	return (j == m);
}

/* Returns length of LCS for X[0..m-1], Y[0..n-1] */
int LCSLength(const string& X, const string& Y, int m, int n, int** L) {
// Build L[m+1][n+1] in bottom up fashion
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			if (i == 0 || j == 0)
				L[i][j] = 0;
			else if (X[i - 1] == Y[j - 1])
				L[i][j] = L[i - 1][j - 1] + 1;
			else
				L[i][j] = max(L[i - 1][j], L[i][j - 1]);
		}
	}
	return L[m][n];
}

vector<string> LCSDP(const string& X, const string& Y, int m, int n, vector<vector<vector<string>>>& L) {
// Build L[m+1][n+1] in bottom up fashion
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			vector<string>& current = L[i][j];
			if (i == 0 || j == 0)
			current.push_back("");
			else if (X[i - 1] == Y[j - 1]) {
				current = L[i - 1][j - 1];
				for (unsigned k = 0; k < current.size(); k++) {
					current[k] += X[i - 1];
				}
			}
			else {
				int upLen = L[i - 1][j][0].size();
				int leftLen = L[i][j - 1][0].size();
				if (upLen > leftLen) {
					current = L[i - 1][j];
					for (unsigned k = 0; k < current.size(); k++) {
						current[k] += X[i - 1];
					}
				}
				else if (upLen < leftLen) {
					current = L[i][j - 1];
					for (unsigned k = 0; k < current.size(); k++) {
						current[k] += Y[j - 1];
					}
				}
				else {	//upLen==leftLen
					current = L[i - 1][j];
					for (unsigned k = 0; k < current.size(); k++) {
						current[k] += X[i - 1];
					}
					vector<string> tmp = L[i][j - 1];
					for (unsigned k = 0; k < tmp.size(); k++) {
						tmp[k] += Y[j - 1];
					}
					current.insert(current.end(),tmp.begin(),tmp.end());
				}
			}
		}
	}
	return L[m][n];
}

vector<string> LCS2DP(const string& X, const string& Y) {
	//allocate memory
	int m = X.length(), n = Y.length();
	vector<vector<vector<string>>> L(m+1, vector<vector<string>>(n+1, vector<string>()));

//	L = new vector<string>*[m + 1];
//	for (int i = 0; i < m + 1; i++)
//		L[i] = new vector<string> [n + 1];

	vector<string> result = LCSDP(X, Y, m, n, L);
	// release memory
//	for (int i = 0; i < m + 1; i++)
//		delete[] L[i];
//	delete[] L;
	return result;
}

/* Returns set containing all LCS for X[0..m-1], Y[0..n-1] */
vector<string> findLCS(const string& X, const string& Y, int m, int n, int** L) {
	// construct a set to store possible LCS
	vector<string> s;

	// If we reach end of either string, return a set with empty string
	if (m == 0 || n == 0) {
		s.push_back("");
		return s;
	}

	// If the last characters of X and Y are the same
	if (X[m - 1] == Y[n - 1]) {
		// recurse for X[0..m-2] and Y[0..n-2] in
		// the matrix
		s = findLCS(X, Y, m - 1, n - 1, L);

		// append current character to all possible LCS
		// of substring X[0..m-2] and Y[0..n-2].
		for (vector<string>::iterator str = s.begin(); str != s.end(); str++)
			str->push_back(X[m - 1]);
	}

	// If the last characters of X and Y are not same
	else {
		// If LCS can be constructed from top side of
		// the matrix, recurse for X[0..m-2] and Y[0..n-1]
		if (L[m - 1][n] >= L[m][n - 1])
			s = findLCS(X, Y, m - 1, n, L);

		// If LCS can be constructed from left side of
		// the matrix, recurse for X[0..m-1] and Y[0..n-2]
		if (L[m][n - 1] >= L[m - 1][n]) {
			vector<string> tmp = findLCS(X, Y, m, n - 1, L);

			// merge two sets if L[m-1][n] == L[m][n-1]
			// Note s will be empty if L[m-1][n] != L[m][n-1]
			s.insert(s.end(), tmp.begin(), tmp.end());
		}
	}
	return s;
}

vector<vector<pair<int, int>>> findLCSIndexes(const string& X, const string& Y, int m, int n, int** L) {
	// construct a set to store possible LCS
	vector<vector<pair<int,int>>> s;

	// If we reach end of either string, return a set with empty string
	if (m == 0 || n == 0) {
		s.push_back(vector<pair<int,int>>());
		return s;
	}

	// If the last characters of X and Y are the same
	if (X[m - 1] == Y[n - 1]) {
		// recurse for X[0..m-2] and Y[0..n-2] in
		// the matrix
		s = findLCSIndexes(X, Y, m - 1, n - 1, L);

		// append current character to all possible LCS
		// of substring X[0..m-2] and Y[0..n-2].
		for (vector<vector<pair<int,int>>>::iterator vec = s.begin(); vec != s.end(); vec++)
		vec->push_back(make_pair(m-1,n-1));
	}

	// If the last characters of X and Y are not same
	else {
		// If LCS can be constructed from top side of
		// the matrix, recurse for X[0..m-2] and Y[0..n-1]
		if (L[m - 1][n] >= L[m][n - 1])
		s = findLCSIndexes(X, Y, m - 1, n, L);

		// If LCS can be constructed from left side of
		// the matrix, recurse for X[0..m-1] and Y[0..n-2]
		if (L[m][n - 1] >= L[m - 1][n]) {
			vector<vector<pair<int,int>>> tmp = findLCSIndexes(X, Y, m, n - 1, L);

			// merge two sets if L[m-1][n] == L[m][n-1]
			// Note s will be empty if L[m-1][n] != L[m][n-1]
			s.insert(s.end(),tmp.begin(), tmp.end());
		}
	}
	return s;
}

vector<string*> findLCSPointers(const string& X, const string& Y, int m, int n, int** L, const int lcsSize) {
	// construct a set to store possible LCS
	vector<string*> s;

	// If we reach end of either string, return a set with empty string
	if (m == 0 || n == 0) {
		string* newString = new string();
		newString->reserve(lcsSize);
		s.push_back(newString);
		return s;
	}

	// If the last characters of X and Y are the same
	if (X[m - 1] == Y[n - 1]) {
		// recurse for X[0..m-2] and Y[0..n-2] in the matrix
		s = findLCSPointers(X, Y, m - 1, n - 1, L, lcsSize);

		// append current character to all possible LCS of substring X[0..m-2] and Y[0..n-2].
		for (vector<string*>::iterator vec = s.begin(); vec != s.end(); vec++)
			*(*vec) += X[m - 1];
	}

	// If the last characters of X and Y are not same
	else {
		// If LCS can be constructed from top side of
		// the matrix, recurse for X[0..m-2] and Y[0..n-1]
		if (L[m - 1][n] >= L[m][n - 1])
			s = findLCSPointers(X, Y, m - 1, n, L, lcsSize);

		// If LCS can be constructed from left side of
		// the matrix, recurse for X[0..m-1] and Y[0..n-2]
		if (L[m][n - 1] >= L[m - 1][n]) {
			vector<string*> tmp = findLCSPointers(X, Y, m, n - 1, L, lcsSize);

			// merge two sets if L[m-1][n] == L[m][n-1]
			// Note s will be empty if L[m-1][n] != L[m][n-1]
			s.insert(s.end(), tmp.begin(), tmp.end());
		}
	}
	return s;
}

int LCS2Len(const string& X, const string& Y) {
	//allocate memory
	int m = X.length(), n = Y.length();
	int **L;
	L = new int*[m + 1];
	for (int i = 0; i < m + 1; i++)
		L[i] = new int[n + 1];

	int lcsLen = LCSLength(X, Y, m, n, L);

	// release memory
	for (int i = 0; i < m + 1; i++)
		delete[] L[i];
	delete[] L;
	return lcsLen;
}

vector<string> LCS2(const string& X, const string& Y) {
	//allocate memory
	int m = X.length(), n = Y.length();
	int **L;
	L = new int*[m + 1];
	for (int i = 0; i < m + 1; i++)
		L[i] = new int[n + 1];

	LCSLength(X, Y, m, n, L);
	vector<string> result = findLCS(X, Y, m, n, L);
	// release memory
	for (int i = 0; i < m + 1; i++)
		delete[] L[i];
	delete[] L;
	return result;
}

vector<vector<pair<int, int>>> LCS2Indexes(const string& X, const string& Y) {
	//allocate memory
	int m = X.length(), n = Y.length();
	int **L;
	L = new int*[m + 1];
	for (int i = 0; i < m + 1; i++)
	L[i] = new int[n + 1];

	LCSLength(X, Y, m, n, L);
	vector<vector<pair<int,int>>> result = findLCSIndexes(X, Y, m, n, L);

	// release memory
	for (int i = 0; i < m + 1; i++)
	delete[] L[i];
	delete[] L;
	return result;
}

vector<string> LCS2Pointers(const string& X, const string& Y) {
	//allocate memory
	int m = X.length(), n = Y.length();
	int **L;
	L = new int*[m + 1];
	for (int i = 0; i < m + 1; i++)
		L[i] = new int[n + 1];

	LCSLength(X, Y, m, n, L);
	vector<string*> temp = findLCSPointers(X, Y, m, n, L, L[m][n]);
	vector<string> result(temp.size());
	for (unsigned int i = 0; i < temp.size(); i++) {
		result[i] = *temp[i];
		delete temp[i];
	}

	// release memory
	for (int i = 0; i < m + 1; i++)
		delete[] L[i];
	delete[] L;
	return result;
}

// create Word by alphabet letter indexes
string Word(const vector<int>& indexes, const vector<char>& alphabet) {
	string result;
	for (unsigned i = 0; i < indexes.size(); i++) {
		result += alphabet[indexes[i]];
	}
	return result;
}

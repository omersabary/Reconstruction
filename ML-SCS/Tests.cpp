#include <iostream>
#include <algorithm>
#include <cassert>
#include <set>
#include <map>
#include <iomanip>
#include <fstream>
#include "SCSN.hpp"
#include "SCS2.hpp"
#include "Tests.hpp"
#include "Strings.hpp"
#include "NChooseK.hpp"

using namespace std;

//void TestKCSNewAndOld(const int reps, const int minStrLen, const int maxStrLen, mt19937& generator) {
//	uniform_int_distribution<int> lenDistribution(minStrLen, maxStrLen);
//	for (int i = 0; i < reps; i++) {
//		int len1 = lenDistribution(generator), len2 = lenDistribution(generator);
//		string X = MakeRandomString(len1, generator), Y = MakeRandomString(len2, generator);
//		int lcsLen = LCS2Len(X, Y);
//		vector<string> kcsdp = KCSDP(X, Y, lcsLen + 1);
//		set<string> setKCSDP(kcsdp.begin(), kcsdp.end());
//		vector<string> kcs = KCS(X, Y, lcsLen + 1);
//		set<string> setKCS(kcs.begin(), kcs.end());
//		assert(setKCSDP == setKCS); // KCS with csubLen more than LCS length should fail;
//		for (int csubLen = 0; csubLen <= lcsLen; csubLen++) {
//			kcsdp = KCSDP(X, Y, csubLen);
//			kcs = KCS(X, Y, csubLen);
//			set<string> setKCSDP(kcsdp.begin(), kcsdp.end());
//			set<string> setKCS(kcs.begin(), kcs.end());
//
//			if (setKCSDP != setKCS) {
//				cout << "Test KCSDP KCS Failure!" << endl;
//				cout << "First String:\t" << X << endl;
//				cout << "Second String:\t" << Y << endl;
//				cout << "csubLen:\t" << csubLen << endl;
//				cout << "KCS:" << endl << kcs;
//				cout << "KCSDP:" << endl << kcsdp;
//				return;
//			}
//
//		}
//	}
//	cout << "Test KCSDP KCS Success!" << endl;
//}

// In what ratio of the cases where regular SCS of 2 strings is too short,
// does the set produced by KCS contain original?
//void TestKCSOriginalFoundFraction(const int reps, const int originalLen, const double delProb, mt19937& generator) {
//	vector<string> copies;
//	double count = 0;
//	double countSuccess = 0;
//	for (int i = 0; i < reps; i++) {
//		string original = RandomCopies(copies, originalLen, 2, delProb, generator);
//		int regSCSNlen = SCSNLen(copies);
//		if (regSCSNlen == originalLen) {
//			continue;
//		}
//		//cout << "was here!" << endl;
//		count++;
//		int csubLen = copies[0].size() + copies[1].size() - originalLen;
//		vector<string> kcs = KCSDP(copies[0], copies[1], csubLen);
//		vector<string>::iterator it = find(kcs.begin(), kcs.end(), original);
//		if (it != kcs.end()) {
//			countSuccess++;
//		}
//	}
//	cout << "Original Found Fraction:\t" << countSuccess / count;
//}

// Test all KCS strings are supersequences of input strings
//void TestKCSIsSuperseq(const int reps, mt19937& generator, const int minStrLen, const int maxStrLen) {
//	uniform_int_distribution<int> lenDistribution(minStrLen, maxStrLen);
//	for (int i = 0; i < reps; i++) {
//		int len1 = lenDistribution(generator), len2 = lenDistribution(generator);
//		string X = MakeRandomString(len1, generator), Y = MakeRandomString(len2, generator);
//		int m = X.size(), n = Y.size();
//		int lcsLen = LCS2Len(X, Y);
//		vector<string> kcs = KCSDP(X, Y, lcsLen + 1);
//		assert(kcs.empty()); // KCS with csubLen more than LCS length should fail;
//		for (int csubLen = 0; csubLen <= lcsLen; csubLen++) {
//			int targetCSLen = m + n - csubLen;
//			kcs = KCSDP(X, Y, csubLen);
//			for (vector<string>::iterator it = kcs.begin(); it != kcs.end(); it++) {
//				assert((int )it->size() == targetCSLen);
//				if (not IsSubSequence(X, *it) or not IsSubSequence(Y, *it)) {
//					cout << "Test KCS Is Subsequence Failure!" << endl;
//					cout << "First String:\t" << X << endl;
//					cout << "Second String:\t" << Y << endl;
//					cout << "csubLen:\t" << csubLen << endl;
//					cout << "KCS:" << endl << kcs;
//					return;
//				}
//			}
//		}
//	}
//	cout << "Test KCS Is Subsequence Success!" << endl;
//}

// test if for csubLen = length of LCS, SCS is the same as normal SCS.
//void TestKCSLCS(const int reps, mt19937& generator, const int minStrLen, const int maxStrLen) {
//	// Create Two Random Strings
//	uniform_int_distribution<int> lenDistribution(minStrLen, maxStrLen);
//	vector<string> strs(2);
//	for (int i = 0; i < reps; i++) {
//		int len1 = lenDistribution(generator), len2 = lenDistribution(generator);
//		string X = MakeRandomString(len1, generator), Y = MakeRandomString(len2, generator);
//		strs[0] = X;
//		strs[1] = Y;
//
//		// Find Common Shortest Supersequences
//		vector<string> SCS = SCSN(strs);
//		set<string> setSCS(SCS.begin(), SCS.end());
//		int scsLen = SCS[0].size();
//		// Find KCS with same length as SCS
//		int csubLen = X.size() + Y.size() - scsLen;
//		vector<string> CS = KCSDP(X, Y, csubLen);
//		set<string> setCS(CS.begin(), CS.end());
//
//		if (setSCS != setCS) {
//			cout << "Test KCS-LCS Failure!" << endl;
//			cout << "First String:\t" << X << endl;
//			cout << "Second String:\t" << Y << endl;
//			cout << "SCS:" << endl << SCS;
//			cout << "KCS:" << endl << CS;
//			return;
//		}
//	}
//	cout << "Test KCS-LCS Success!" << endl;
//}

// Is candidate a result of KCS algorithm
//bool CandidateIsSuperX(const string& X, const string& Y, const string& cand, const int m, const int n, const int p) {
//	if (m == 0) {
//		return Y.substr(0, n) == cand.substr(0, p);
//	}
//	if (n == 0) {
//		return X.substr(0, m) == cand.substr(0, p);
//	}
//	if (p == 0) {
//		return false;
//	}
//	if (cand[p - 1] != X[m - 1] and cand[p - 1] != Y[n - 1]) {
//		return false;
//	}
//	if (X[m - 1] == Y[n - 1]) { //==candidate[p-1]
//		return CandidateIsSuperX(X, Y, cand, m - 1, n - 1, p - 1) or CandidateIsSuperX(X, Y, cand, m, n - 1, p - 1)
//				or CandidateIsSuperX(X, Y, cand, m - 1, n, p - 1);
//	}
//	else { //X[m - 1] != Y[n - 1] and candidate is equal to one of them
//		if (X[m - 1] == cand[p - 1]) {
//			return CandidateIsSuperX(X, Y, cand, m - 1, n, p - 1);
//		}
//		else { //Y[n-1]==cand[p-1];
//			return CandidateIsSuperX(X, Y, cand, m, n - 1, p - 1);
//		}
//	}
//}

//const int NA = -1;
//
//int CandidateIsSuperXArray(const string& X, const string& Y, const string& cand, const int m, const int n, const int p,
//		int*** lookup) {
//	if (lookup[m][n][p] != NA) {
//		return lookup[m][n][p];
//	}
//	if (m == 0) {
//		lookup[m][n][p] = (Y.substr(0, n) == cand.substr(0, p));
//		return lookup[m][n][p];
//	}
//	if (n == 0) {
//		lookup[m][n][p] = (X.substr(0, m) == cand.substr(0, p));
//		return lookup[m][n][p];
//	}
//	if (p == 0) {
//		lookup[m][n][p] = 0;
//		return 0;
//	}
//	if (cand[p - 1] != X[m - 1] and cand[p - 1] != Y[n - 1]) {
//		lookup[m][n][p] = 0;
//		return 0;
//	}
//	if (X[m - 1] == Y[n - 1]) { //==candidate[p-1]
//		lookup[m][n][p] = CandidateIsSuperXArray(X, Y, cand, m - 1, n - 1, p - 1, lookup);
//		if (lookup[m][n][p] == 1) {
//			return 1;
//		}
//		lookup[m][n][p] = CandidateIsSuperXArray(X, Y, cand, m, n - 1, p - 1, lookup);
//		if (lookup[m][n][p] == 1) {
//			return 1;
//		}
//		lookup[m][n][p] = CandidateIsSuperXArray(X, Y, cand, m - 1, n, p - 1, lookup);
//		return lookup[m][n][p];
//	}
//	else { //X[m - 1] != Y[n - 1] and candidate is equal to one of them
//		if (X[m - 1] == cand[p - 1]) {
//			lookup[m][n][p] = CandidateIsSuperXArray(X, Y, cand, m - 1, n, p - 1, lookup);
//			return lookup[m][n][p];
//		}
//		else { //Y[n-1]==cand[p-1];
//			lookup[m][n][p] = CandidateIsSuperXArray(X, Y, cand, m, n - 1, p - 1, lookup);
//			return lookup[m][n][p];
//		}
//	}
//}

//bool CandidateIsSupX(const string& X, const string& Y, const string& cand) {
//	int ***lookup;
//	int m = X.size(), n = Y.size(), p = cand.size();
//	lookup = new int**[m + 1];
//	for (int i = 0; i < m + 1; i++) {
//		lookup[i] = new int*[n + 1];
//		for (int j = 0; j < n + 1; j++) {
//			lookup[i][j] = new int[p + 1];
//			for (int k = 0; k < p + 1; k++) {
//				lookup[i][j][k] = NA;
//			}
//		}
//	}
//
//	int result = CandidateIsSuperXArray(X, Y, cand, m, n, p, lookup);
//
//	// Release memory
//	for (int i = 0; i < m + 1; i++) {
//		for (int j = 0; j < n + 1; j++) {
//			delete[] lookup[i][j];
//		}
//		delete[] lookup[i];
//	}
//	delete[] lookup;
//	return result == 1;
//}

// Test KCS has all common sequences it is supposed to have
//void TestAllKCS(const int reps, mt19937& generator, const int minStrLen, const int maxStrLen) {
//	// Create Two Random Strings
//	uniform_int_distribution<int> lenDistribution(minStrLen, maxStrLen);
//	vector<string> strs(2);
//	vector<char> alphabet = { 'G', 'T', 'A', 'C' };
//	for (int i = 0; i < reps; i++) {
//		int len1 = lenDistribution(generator), len2 = lenDistribution(generator);
//		string X = MakeRandomString(len1, generator), Y = MakeRandomString(len2, generator);
//		strs[0] = X;
//		strs[1] = Y;
//
//		// Find Common Shortest Supersequences
//		int lcsLen = LCS2Len(X, Y);
//		if (lcsLen == 0) {
//			continue;
//		}
//		//uniform_int_distribution<int> cSubLenDistribution(0, lcsLen);
//		//int csubLen = cSubLenDistribution(generator);
//		for (int csubLen = 0; csubLen <= lcsLen; csubLen++) {
//			int csupLen = X.size() + Y.size() - csubLen;
//			vector<string> CS = KCSDP(X, Y, csubLen);
//			set<string> setCS(CS.begin(), CS.end());
//
//			vector<int> dimSizes(csupLen);
//			for (int i = 0; i < csupLen; i++) {
//				dimSizes[i] = 4;
//			}
//			set<string> candidates;
//			IndexVector indexVector(dimSizes);
//			for (IndexVector currentIndex = indexVector; not currentIndex.IsPastEnd(); currentIndex++) {
//				string word = Word(currentIndex.Vector(), alphabet);
//				if (CandidateIsSuperX(X, Y, word, X.size(), Y.size(), csupLen)) {
//					candidates.insert(word);
//				}
//			}
//
//			if (candidates != setCS) {
//				cout << "---------------" << endl;
//				cout << "Test All KCS Failure!" << endl;
//				cout << "X:\t" << X << endl;
//				cout << "Y:\t" << Y << endl;
//				cout << "KCS:\t" << endl << setCS;
//				cout << "Candidates:" << endl << candidates;
//				return;
//			}
//		}
//	}
//	cout << "Test All KCS Success!" << endl;
//}

void TestSCSNum(const int reps, const int originalLen, const int K, const double delProb, mt19937& generator) {
	vector<string> copies;
	double countCorrect = 0, correctNum = 0, countShort = 0, shortNum = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, K, delProb, generator);
		vector<string> scs = SCSNFast(copies, originalLen);
		int scsSize = scs[0].size();
		int scsNum = scs.size();
		if (scsSize == originalLen) {
			countCorrect += scsNum;
			correctNum++;
		}
		else {
			countShort += scsNum;
			shortNum++;
		}

	}
	cout << "Fraction of correct Size SCS:\t\t" << correctNum / reps << endl;
	cout << "Average number of correct size SCS:\t" << countCorrect / correctNum << endl;
	cout << "Average number of too short SCS:\t" << countShort / shortNum << endl;
}

void TestKLongestSCSNum(const int reps, const int originalLen, const int copiesNum, const int K, const double delProb,
		mt19937& generator) {
	vector<string> copies, BestK;
	double countCorrect = 0, correctNum = 0, countShort = 0, shortNum = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		BestK = ChooseKLongest(copies, K);
		vector<string> scs = SCSNFast(BestK, originalLen);
		int scsSize = scs[0].size();
		int scsNum = scs.size();
		if (scsSize == originalLen) {
			countCorrect += scsNum;
			correctNum++;
		}
		else {
			countShort += scsNum;
			shortNum++;
		}
	}
	cout << "Fraction of correct Size SCS:\t\t" << correctNum / reps << endl;
	cout << "Average number of correct size SCS:\t" << countCorrect / correctNum << endl;
	cout << "Average number of too short SCS:\t" << countShort / shortNum << endl;
}

void TestMinSCSNum(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	vector<string> copies;
	double countCorrect = 0, correctNum = 0; // countShort = 0, shortNum = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		vector<pair<string, string> > pairs = StringPairings(copies);
		int pairNum = pairs.size();
		vector<pair<int, int>> scsNumLen(pairNum);
		for (int i = 0; i < pairNum; i++) {
			vector<string> pairVec = { pairs[i].first, pairs[i].second };
			vector<string> scs = SCSNFast(pairVec, originalLen);
			int scsSize = scs[0].size();
			int scsNum = scs.size();
			scsNumLen[i] = make_pair(scsNum, scsSize);
		}
		int minCorrect = INT_MAX;
		//int minShort = INT_MAX;
		for (int i = 0; i < pairNum; i++) {
			int scsSize = scsNumLen[i].second;
			int scsNum = scsNumLen[i].first;
			if (scsSize == originalLen) {
				if (scsNum < minCorrect) {
					minCorrect = scsNum;
				}
			}
//			else {
//				if (scsNum < minShort) {
//					minShort = scsNum;
//				}
//			}
		}
		if (minCorrect != INT_MAX) {
			countCorrect += minCorrect;
			correctNum++;
		}
//		else {
//			countShort += minShort;
//			shortNum++;
//		}

	}
	//cout << "Fraction of correct Size SCS:\t\t" << correctNum / reps << endl;
	cout << "Average number of min correct size SCS:\t" << countCorrect / correctNum << endl;
	//cout << "Average number of too short SCS:\t" << countShort / shortNum << endl;
}

void TestCorrectSizeFraction(const int reps, const int originalLen, const int copiesNum, const int k,
		const double delProb, mt19937& generator) {
	vector<string> copies;
	double countCorrect = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		vector<vector<string>> ktuples = StringKtuples(copies, k);
		int ktuplesNum = ktuples.size();
		vector<int> scsLen(ktuplesNum);
		for (int i = 0; i < ktuplesNum; i++) {
			scsLen[i] = SCSNLenFast(ktuples[i], originalLen);
		}
		int maxSCSLen = *max_element(scsLen.begin(), scsLen.end());
		if (maxSCSLen == originalLen) {
			countCorrect++;
		}
	}
	cout << "Fraction of correct Size SCS:\t\t" << countCorrect / reps << endl;
}

int ChoseFirst(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int chosenIndex = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			chosenIndex = i;
			break;
		}
	}
	return chosenIndex;
}

int ChoseMaxLenDiff(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int maxLenDiff = 0;
	int indexOfMaxPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			if (abs(firstLen - secondLen) > maxLenDiff) {
				maxLenDiff = abs(firstLen - secondLen);
				indexOfMaxPair = i;
			}
		}
	}
	return indexOfMaxPair;
}

int ChoseMinLenDiff(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int minLenDiff = INT_MAX;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			if (abs(firstLen - secondLen) < minLenDiff) {
				minLenDiff = abs(firstLen - secondLen);
				indexOfMinPair = i;
			}
		}
	}
	return indexOfMinPair;
}

int ChoseMaxSumLen(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int maxSum = 0;
	int indexOfMaxPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			if (firstLen + secondLen > maxSum) {
				maxSum = firstLen + secondLen;
				indexOfMaxPair = i;
			}
		}
	}
	return indexOfMaxPair;
}

//int ChoseMaxLen(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
//	int pairNum = pairs.size();
//	int maxLen = 0;
//	int indexOfMaxPair = -1;
//	for (int i = 0; i < pairNum; i++) {
//		if (scsLen[i] == originalLen) {
//			int firstLen = pairs[i].first.size();
//			int secondLen = pairs[i].second.size();
//			if (max(firstLen,secondLen) > maxLen) {
//				maxLen = max(firstLen,secondLen);
//				indexOfMaxPair = i;
//			}
//		}
//	}
//	return indexOfMaxPair;
//}

int ChoseMaxSumLenOnEqualMinStr(const vector<pair<string, string> >& pairs, const vector<int>& scsLen,
		int originalLen) {
	int pairNum = pairs.size();
	int maxSum = 0;
	int minStrLen = INT_MAX;
	int indexOfMaxPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			if (firstLen + secondLen > maxSum) {
				maxSum = firstLen + secondLen;
				minStrLen = min(firstLen, secondLen);
				indexOfMaxPair = i;
			}
			else if (firstLen + secondLen == maxSum) {
				int minLen = min(firstLen, secondLen);
				if (minLen < minStrLen) {
					minStrLen = minLen;
					indexOfMaxPair = i;
				}
			}
		}
	}
	return indexOfMaxPair;
}

int ChoseMaxSumLenOnEqualMinED(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int maxSum = 0;
	int minEditDistance = INT_MAX;
	int indexOfMaxPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			if (firstLen + secondLen > maxSum) {
				maxSum = firstLen + secondLen;
				minEditDistance = EditDistance(pairs[i].first, pairs[i].second);
				indexOfMaxPair = i;
			}
			else if (firstLen + secondLen == maxSum) {
				int currentEditDistance = EditDistance(pairs[i].first, pairs[i].second);
				if (currentEditDistance < minEditDistance) {
					minEditDistance = currentEditDistance;
					indexOfMaxPair = i;
				}
			}
		}
	}
	return indexOfMaxPair;
}

int ChoseMaxSumLenOnEqualMaxED(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int maxSum = 0;
	int maxEditDistance = 0;
	int indexOfMaxPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			if (firstLen + secondLen > maxSum) {
				maxSum = firstLen + secondLen;
				maxEditDistance = EditDistance(pairs[i].first, pairs[i].second);
				indexOfMaxPair = i;
			}
			else if (firstLen + secondLen == maxSum) {
				int currentEditDistance = EditDistance(pairs[i].first, pairs[i].second);
				if (currentEditDistance > maxEditDistance) {
					maxEditDistance = currentEditDistance;
					indexOfMaxPair = i;
				}
			}
		}
	}
	return indexOfMaxPair;
}

int ChoseMinSumLen(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int minSum = INT_MAX;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			if (firstLen + secondLen < minSum) {
				minSum = firstLen + secondLen;
				indexOfMinPair = i;
			}
		}
	}
	return indexOfMinPair;
}

int ChoseMinEditDistance(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int minEditDistance = INT_MAX;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int currentED = EditDistance(pairs[i].first, pairs[i].second);
			if (currentED < minEditDistance) {
				minEditDistance = currentED;
				indexOfMinPair = i;
			}
		}
	}
	return indexOfMinPair;
}

int ChoseMaxEditDistance(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int maxEditDistance = 0;
	int indexOfMaxPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int currentED = EditDistance(pairs[i].first, pairs[i].second);
			if (currentED > maxEditDistance) {
				maxEditDistance = currentED;
				indexOfMaxPair = i;
			}
		}
	}
	return indexOfMaxPair;
}

// chose pair with minimum difference between Insertions deletions edit distance and insertion deletions
// substitutions edit distance
int ChoseMinEDDiff(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int minEditDistanceDiff = INT_MAX;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			int lcslen = firstLen + secondLen - scsLen[i];
			int currentEDID = firstLen + secondLen - 2 * lcslen;
			int currentED = EditDistance(pairs[i].first, pairs[i].second);
			int currentDiff = currentEDID - currentED;
			assert(currentDiff >= 0);
			if (currentDiff < minEditDistanceDiff) {
				minEditDistanceDiff = currentDiff;
				indexOfMinPair = i;
			}
		}
	}
	return indexOfMinPair;
}

int ChoseMinEDDiffOnEqualMaxSumLen(const vector<pair<string, string> >& pairs, const vector<int>& scsLen,
		int originalLen) {
	int pairNum = pairs.size();
	int minEditDistanceDiff = INT_MAX;
	int maxSum = 0;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			int lcslen = firstLen + secondLen - scsLen[i];
			int currentEDID = firstLen + secondLen - 2 * lcslen;
			int currentED = EditDistance(pairs[i].first, pairs[i].second);
			int currentDiff = currentEDID - currentED;
			assert(currentDiff >= 0);
			if (currentDiff < minEditDistanceDiff) {
				minEditDistanceDiff = currentDiff;
				maxSum = firstLen + secondLen;
				indexOfMinPair = i;
			}
			else if (currentDiff == minEditDistanceDiff) {
				if (firstLen + secondLen > maxSum) {
					maxSum = firstLen + secondLen;
					indexOfMinPair = i;
				}
			}
		}
	}
	return indexOfMinPair;
}

int ChoseMinEDDiffOnEqualMinSumLen(const vector<pair<string, string> >& pairs, const vector<int>& scsLen,
		int originalLen) {
	int pairNum = pairs.size();
	int minEditDistanceDiff = INT_MAX;
	int minSum = INT_MAX;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			int lcslen = firstLen + secondLen - scsLen[i];
			int currentEDID = firstLen + secondLen - 2 * lcslen;
			int currentED = EditDistance(pairs[i].first, pairs[i].second);
			int currentDiff = currentEDID - currentED;
			assert(currentDiff >= 0);
			if (currentDiff < minEditDistanceDiff) {
				minEditDistanceDiff = currentDiff;
				minSum = firstLen + secondLen;
				indexOfMinPair = i;
			}
			else if (currentDiff == minEditDistanceDiff) {
				if (firstLen + secondLen < minSum) {
					minSum = firstLen + secondLen;
					indexOfMinPair = i;
				}
			}
		}
	}
	return indexOfMinPair;
}

int ChoseMinEDDiffOnEqualMaxLen(const vector<pair<string, string> >& pairs, const vector<int>& scsLen,
		int originalLen) {
	int pairNum = pairs.size();
	int minEditDistanceDiff = INT_MAX;
	int maxLen = 0;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			int lcslen = firstLen + secondLen - scsLen[i];
			int currentEDID = firstLen + secondLen - 2 * lcslen;
			int currentED = EditDistance(pairs[i].first, pairs[i].second);
			int currentDiff = currentEDID - currentED;
			assert(currentDiff >= 0);
			if (currentDiff < minEditDistanceDiff) {
				minEditDistanceDiff = currentDiff;
				maxLen = max(firstLen, secondLen);
				indexOfMinPair = i;
			}
			else if (currentDiff == minEditDistanceDiff) {
				if (max(firstLen, secondLen) > maxLen) {
					maxLen = max(firstLen, secondLen);
					indexOfMinPair = i;
				}
			}
		}
	}
	return indexOfMinPair;
}

int ChoseMinEDDiffOnEqualMinLen(const vector<pair<string, string> >& pairs, const vector<int>& scsLen,
		int originalLen) {
	int pairNum = pairs.size();
	int minEditDistanceDiff = INT_MAX;
	int minLen = INT_MAX;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			int lcslen = firstLen + secondLen - scsLen[i];
			int currentEDID = firstLen + secondLen - 2 * lcslen;
			int currentED = EditDistance(pairs[i].first, pairs[i].second);
			int currentDiff = currentEDID - currentED;
			assert(currentDiff >= 0);
			if (currentDiff < minEditDistanceDiff) {
				minEditDistanceDiff = currentDiff;
				minLen = min(firstLen, secondLen);
				indexOfMinPair = i;
			}
			else if (currentDiff == minEditDistanceDiff) {
				if (min(firstLen, secondLen) < minLen) {
					minLen = min(firstLen, secondLen);
					indexOfMinPair = i;
				}
			}
		}
	}
	return indexOfMinPair;
}

int ChoseMinEDDiffOnEqualMinED(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int minEditDistanceDiff = INT_MAX;
	int minED = INT_MAX;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			int lcslen = firstLen + secondLen - scsLen[i];
			int currentEDID = firstLen + secondLen - 2 * lcslen;
			int currentED = EditDistance(pairs[i].first, pairs[i].second);
			int currentDiff = currentEDID - currentED;
			assert(currentDiff >= 0);
			if (currentDiff < minEditDistanceDiff) {
				minEditDistanceDiff = currentDiff;
				minED = currentED;
				indexOfMinPair = i;
			}
			else if (currentDiff == minEditDistanceDiff) {
				if (currentED < minED) {
					minED = currentED;
					indexOfMinPair = i;
				}
			}
		}
	}
	return indexOfMinPair;
}

int ChoseMinEDDiffOnEqualMaxED(const vector<pair<string, string> >& pairs, const vector<int>& scsLen, int originalLen) {
	int pairNum = pairs.size();
	int minEditDistanceDiff = INT_MAX;
	int maxED = 0;
	int indexOfMinPair = -1;
	for (int i = 0; i < pairNum; i++) {
		if (scsLen[i] == originalLen) {
			int firstLen = pairs[i].first.size();
			int secondLen = pairs[i].second.size();
			int lcslen = firstLen + secondLen - scsLen[i];
			int currentEDID = firstLen + secondLen - 2 * lcslen;
			int currentED = EditDistance(pairs[i].first, pairs[i].second);
			int currentDiff = currentEDID - currentED;
			assert(currentDiff >= 0);
			if (currentDiff < minEditDistanceDiff) {
				minEditDistanceDiff = currentDiff;
				maxED = currentED;
				indexOfMinPair = i;
			}
			else if (currentDiff == minEditDistanceDiff) {
				if (currentED > maxED) {
					maxED = currentED;
					indexOfMinPair = i;
				}
			}
		}
	}
	return indexOfMinPair;
}

// Tests: copies num 5, len 100, delProb 0.10.
// First pair with correct size scs: reps 10000: 26.32.
// pair with max sum of copies len: reps 10000: 18.07 10000: 15.4

// find pairs whose scs is of correct size. if more than one - choose longest pair.
void TestChoseFromCorrectSizeSCSs(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	vector<string> copies;
	double countCorrect = 0, correctNum = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		vector<pair<string, string> > pairs = StringPairings(copies);
		int pairNum = pairs.size();
		vector<int> scsLen(pairNum);
		// create vector of scs length for each string pair
		for (int i = 0; i < pairNum; i++) {
			vector<string> pairVec = { pairs[i].first, pairs[i].second };
			int currentScsLen = SCSNLenFast(pairVec, originalLen);
			scsLen[i] = currentScsLen;
		}

//		int indexOfChosenPair = ChoseFirst(pairs, scsLen, originalLen);
		int indexOfChosenPair = ChoseMaxSumLen(pairs, scsLen, originalLen);
//		int indexOfChosenPair = ChoseFirst(pairs,scsLen, originalLen);
//		int indexOfChosenPair = ChoseMaxLenDiff(pairs, scsLen, originalLen);
//		int indexOfChosenPair = ChoseMinLenDiff(pairs, scsLen, originalLen);
//		int indexOfChosenPair = ChoseMaxSumLenOnEqualMaxED(pairs, scsLen, originalLen);
//		int indexOfChosenPair = ChoseMinEDDiff(pairs, scsLen, originalLen);
//		int indexOfChosenPair = ChoseMinEDDiffOnEqualMaxSumLen(pairs, scsLen, originalLen);
//		int indexOfChosenPair = ChoseMinEDDiffOnEqualMaxSumLen(pairs, scsLen, originalLen);

		if (indexOfChosenPair >= 0) {		// there was at least one pair with correct size SCS
			vector<string> pairVec = { pairs[indexOfChosenPair].first, pairs[indexOfChosenPair].second };
			vector<string> scs = SCSNFast(pairVec, originalLen);
			countCorrect += scs.size(); // accumulate number of scss for chosen pair
			correctNum++;
		}
	}

	cout << "Fraction of correct Size SCS:\t\t" << correctNum / reps << endl;
	cout << "Average number of correct size SCS:\t" << countCorrect / correctNum << endl;
}

// Find Letter Error Rate of most likely SCS.

void TestMostLikelySCS(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	vector<string> copies;
	double cumEditDistance = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		string mostLikelySCS = MostLikelySCS(copies, delProb, originalLen, generator);
		cumEditDistance += EditDistanceDI(original, mostLikelySCS);
	}
	cout << "Letter Error Rate:\t" << cumEditDistance / (reps * originalLen) << endl;
}

// only for 2 copies for now
//void TestMostLikelyKCS(const int reps, const int originalLen, const double delProb, mt19937& generator) {
//	vector<string> copies;
//	double cumEditDistance = 0;
//	for (int i = 0; i < reps; i++) {
//		string original = RandomCopies(copies, originalLen, 2, delProb, generator);
//		string mostLikelySCS = MostLikelyKCS(copies, delProb, originalLen, generator);
//		cumEditDistance += EditDistanceDI(original, mostLikelySCS);
//	}
//	cout << "KCS Letter Error Rate:\t" << cumEditDistance / (reps * originalLen) << endl;
//}

void TestMostLikelyFilteredSCS(const int reps, const int originalLen, const int copiesNum, const int k,
		const double delProb, mt19937& generator) {
	vector<string> copies;
	map<int, int> SCSofLenCount, cumSCSNum, cumMostLikelyNum, cumEditDistance;
	int currentSCSNum = 0, currentMostLikelyNum = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		string mostLikelySCS = MostLikelyFilteredSCS(copies, k, delProb, originalLen);
		int scsLen = mostLikelySCS.size();
		SCSofLenCount[scsLen]++;
		cumSCSNum[scsLen] += currentSCSNum;
		cumMostLikelyNum[scsLen] += currentMostLikelyNum;
		cumEditDistance[scsLen] += EditDistanceDI(original, mostLikelySCS);
	}
	cout << "Copies Num:\t\t" << copiesNum << endl;
	cout << "Original Len:\t\t" << originalLen << endl;
	cout << "Deletion Probability:\t\t" << delProb << endl;
	cout << "Repetitions:\t\t" << reps << endl;
	cout << "SCS Len\tOccurences\tSCS Total\tML SCS Total\tTotal ED" << endl;
	for (map<int, int>::reverse_iterator rit = SCSofLenCount.rbegin(); rit != SCSofLenCount.rend(); rit++) {
		int scsLen = rit->first;
		cout << scsLen << "\t" << SCSofLenCount[scsLen] << "\t" << cumSCSNum[scsLen] << "\t" << cumMostLikelyNum[scsLen]
				<< "\t" << cumEditDistance[scsLen] << endl;
	}
}

void TestAlgorithm2(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	vector<string> copies;
	double cumEditDistance = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		string guess = MostLikelyMaxLenSCS(copies, 2, delProb, originalLen);
		cumEditDistance += EditDistanceDI(original, guess);
	}
	cout << fixed << setprecision(8) << "Algorithm 2 Letter Error Rate:\t" << cumEditDistance / (reps * originalLen)
			<< setprecision(2) << endl;
}

void TestAlgorithm3(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	vector<string> copies;
	double cumEditDistance = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		string guess = MostLikelyCorrectSizeSCS(copies, 2, delProb, originalLen);
		if (guess.empty()) {
			guess = MostLikelyMaxLenSCS(copies, 3, delProb, originalLen);
		}
		cumEditDistance += EditDistanceDI(original, guess);
	}
	cout << fixed << setprecision(8) << "Algorithm 3 Letter Error Rate:\t" << cumEditDistance / (reps * originalLen)
			<< setprecision(2) << endl;
}

void TestBigClusterAlgorithm(const int reps, const int originalLen, const int copiesNum, const double delProb, const int k,
		mt19937& generator) {
	vector<string> copies;
	vector<int> twoHist(copiesNum * (copiesNum - 1) / 2);
	int nothing = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		// which pair number in sorted pairs gives original length
		int ktupleNum = CorrectSizeKtupleNumSortBySumLen(copies,k,originalLen);
		//int ktupleNum = CorrectSizeKtupleNumSortByLetterMax(copies,k,originalLen);
		if (ktupleNum==-1){
			nothing++;
		}
		else{
			twoHist[ktupleNum]++;
		}
	}
	cout << "Original Len:\t" << originalLen << endl;
	cout << "Copies Num:\t" << copiesNum << endl;
	cout << "Deletion Prob:\t" << delProb << endl;
	cout << "Reps:\t" << reps << endl;
	int cumHist = 0;
	for (unsigned i = 0; i < twoHist.size(); i++) {
		cumHist+=twoHist[i];
		cout << i << "\t" << cumHist << endl;
	}
	cout << "None:\t" << nothing << endl;
}

void TestAlgorithmNStepsOld(const int reps, const int maxStageNum, const int originalLen, const int copiesNum,
		const double delProb, mt19937& generator) {
	vector<string> copies;
	string guess;
	int stageNum;
	map<int, int> editDistanceDist, stoppingStepDist;

	double cumEditDistance = 0, successNum = 0;
	for (int i = 0; i < reps; i++) {
		stageNum = 2;
		guess.clear();
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		while (guess.empty() and stageNum < maxStageNum) {
			guess = MostLikelyCorrectSizeSCS(copies, stageNum, delProb, originalLen);
			stageNum++;
		}
		if (guess.empty()) { // stageNum == maxStageNum
			guess = MostLikelyMaxLenSCS(copies, maxStageNum, delProb, originalLen);
			stoppingStepDist[maxStageNum]++;
		}
		else { // guess not empty
			stoppingStepDist[stageNum - 1]++;
		}
		int currentEditDistance = EditDistanceDI(original, guess);
		editDistanceDist[currentEditDistance]++;
		cumEditDistance += currentEditDistance;
		if (guess == original) {
			successNum++;
		}
	}
	cout << "Original Len:\t" << originalLen << endl;
	cout << "Copies Num:\t" << copiesNum << endl;
	cout << "Deletion Prob:\t" << setprecision(2) << delProb << endl;
	cout << "Max Stage Num:\t" << maxStageNum << endl;
	cout << "Reps:\t" << reps << endl;
	cout << fixed << setprecision(8) << "Letter Error Rate:\t" << cumEditDistance / (reps * originalLen)
			<< setprecision(6) << endl;
	cout << "Success Rate:\t" << successNum / reps << endl;
	cout << "Guess ED Distribution:" << endl;
	for (map<int, int>::iterator it = editDistanceDist.begin(); it != editDistanceDist.end(); it++) {
		cout << it->first << "\t" << it->second << endl;
	}
	cout << "Stopping Stage Distribution:" << endl;
	for (map<int, int>::iterator it = stoppingStepDist.begin(); it != stoppingStepDist.end(); it++) {
		cout << it->first << "\t" << it->second << endl;
	}
	cout << "********************" << endl;
}

void TestAlgorithmNSteps(const int reps, const int maxStageNum, const int originalLen, const int copiesNum,
		const double delProb, mt19937& generator) {
	string clusterSizeStr = to_string(copiesNum);
	string deletionProbStr = to_string((int) ((delProb + 0.001) * 100));
	string testNumStr = to_string(reps);
	string originalLenStr = to_string(originalLen);
	string fileName = clusterSizeStr + "_" + deletionProbStr + "_" + originalLenStr + "_" + testNumStr;
	cout << fileName << endl;
	ofstream file(fileName + ".txt");
	if (not file.is_open()) {
		cout << "Error opening output file!" << endl;
		return;
	}
	vector<string> copies;
	string guess;
	int nothingCount = 0; // how many time 2,3,4 didn't find correct size;
	int stageNum, foundKtupleIndex;
	map<int, int> editDistanceHist;
	vector<int> stoppingStageHist(maxStageNum + 1);
	vector<int> twoHist(copiesNum * (copiesNum - 1) / 2);
	vector<int> threeHist(copiesNum * (copiesNum - 1) * (copiesNum - 2) / 6);
	vector<int> fourHist(copiesNum * (copiesNum - 1) * (copiesNum - 2) * (copiesNum - 3) / 24);

	double cumEditDistance = 0, successNum = 0;
	for (int i = 0; i < reps; i++) {
		stageNum = 2;
		guess.clear();
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		while (guess.empty() and stageNum < maxStageNum) {
			guess = MostLikelyCorrectSizeSCS(copies, stageNum, delProb, originalLen, foundKtupleIndex);
			stageNum++;
		}
		bool foundCorrectSize;
		if (guess.empty()) { // stageNum == maxStageNum
			guess = MostLikelyMaxLenSCS(copies, maxStageNum, delProb, originalLen, foundCorrectSize, foundKtupleIndex);
			if (foundCorrectSize) {
				stoppingStageHist[stageNum]++;
				fourHist[foundKtupleIndex]++;
			}
			else {
				nothingCount++;
			}
		}
		else { // guess not empty
			stoppingStageHist[stageNum - 1]++;
			if (stageNum - 1 == 2) {
				twoHist[foundKtupleIndex]++;
			}
			if (stageNum - 1 == 3) {
				threeHist[foundKtupleIndex]++;
			}
		}
		int currentEditDistance = EditDistanceDI(original, guess);
		editDistanceHist[currentEditDistance]++;
		cumEditDistance += currentEditDistance;
		if (guess == original) {
			successNum++;
		}
	}
	cout << "Cluster size:\t" << copiesNum << endl;
	cout << "Deletion Prob:\t" << setprecision(2) << delProb << endl;
	cout << "Number of tests:\t" << reps << endl;
	cout << "str Len:\t" << originalLen << endl;
	cout << "q:\t" << maxStageNum << endl;
	cout << fixed << setprecision(8) << "Error rate:\t" << cumEditDistance / (reps * originalLen) << setprecision(6)
			<< endl;
	cout << "Success rate:\t" << successNum / reps << endl;

	file << "Cluster size:\t" << copiesNum << endl;
	file << "Deletion Prob:\t" << setprecision(2) << delProb << endl;
	file << "Number of tests:\t" << reps << endl;
	file << "str Len:\t" << originalLen << endl;
	file << "q:\t" << maxStageNum << endl;
	file << fixed << setprecision(8) << "Error rate:\t" << cumEditDistance / (reps * originalLen) << setprecision(6)
			<< endl;
	file << "Success rate:\t" << successNum / reps << endl;

	// Stopping Stage Hist
	cout << "Two: " << stoppingStageHist[2] << endl;
	cout << "Three: " << stoppingStageHist[3] << endl;
	cout << "Four: " << stoppingStageHist[4] << endl;
	cout << "Nothing: " << nothingCount << endl;

	file << "Two: " << stoppingStageHist[2] << endl;
	file << "Three: " << stoppingStageHist[3] << endl;
	file << "Four: " << stoppingStageHist[4] << endl;
	file << "Nothing: " << nothingCount << endl;

	cout << "SR histo:" << endl;
	file << "SR histo:" << endl;
	map<int, int>::reverse_iterator rit = editDistanceHist.rbegin(); // points to last element in map
	int highestED = rit->first;
	int cumDist = 0;
	for (int i = 0; i <= highestED; i++) {
		cumDist += editDistanceHist[i];
		cout << i << "\t" << cumDist << endl;
		file << i << "\t" << cumDist << endl;
	}

	cout << "Two histo:" << endl;
	file << "Two histo:" << endl;
	for (unsigned i = 0; i < twoHist.size(); i++) {
		cout << i << "\t" << twoHist[i] << endl;
		file << i << "\t" << twoHist[i] << endl;
	}
	cout << "Three histo:" << endl;
	file << "Three histo:" << endl;
	for (unsigned i = 0; i < threeHist.size(); i++) {
		cout << i << "\t" << threeHist[i] << endl;
		file << i << "\t" << threeHist[i] << endl;
	}
	cout << "Four histo:" << endl;
	file << "Four histo:" << endl;
	for (unsigned i = 0; i < fourHist.size(); i++) {
		cout << i << "\t" << fourHist[i] << endl;
		file << i << "\t" << fourHist[i] << endl;
	}
	cout << "********************" << endl;
	file.close();
}

// Test only cases in which no 2 pairs have scs of correct size
// hardNum - number of hard cases
void TestAlgorithmNStepsHard(const int hardNum, const int maxStageNum, const int originalLen, const int copiesNum,
		const double delProb, mt19937& generator) {
	vector<string> copies;
	string guess;
	int stageNum;
	double cumEditDistance = 0, successNum = 0;
	int caseNum = 0;
	while (caseNum < hardNum) {
		guess.clear();
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		if (HasPairWithCorrectSizeSCS(copies, originalLen)) {
			continue;
		}
		stageNum = 3;
		while (guess.empty() and stageNum < maxStageNum) {
			guess = MostLikelyCorrectSizeSCS(copies, stageNum, delProb, originalLen);
			stageNum++;
		}
		if (guess.empty()) {
			guess = MostLikelyMaxLenSCS(copies, maxStageNum, delProb, originalLen);
		}
		cumEditDistance += EditDistanceDI(original, guess);
		if (guess == original) {
			successNum++;
		}
		caseNum++;
	}
	cout << "Original Len:\t" << originalLen << endl;
	cout << "Copies Num:\t" << copiesNum << endl;
	cout << "Deletion Prob:\t" << setprecision(3) << delProb << endl;
	cout << "Stage Num:\t" << maxStageNum << endl;
	cout << "Hard Cases Num:\t" << hardNum << endl;
	cout << fixed << setprecision(8) << "Letter Error Rate:\t" << cumEditDistance / (hardNum * originalLen)
			<< setprecision(6) << endl;
	cout << "Success Rate:\t" << successNum / hardNum << endl;
	cout << "********************" << endl;
}

void TestSCSEstimate(const int reps, const int originalLen, const double delProb, mt19937& generator) {
	vector<string> copies;
	double cumRealNum = 0, cumEstimatedNum = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, 2, delProb, generator);
		vector<string> scss = SCS2(copies[0], copies[1]);
		int realSCSNum = scss.size();
		cumRealNum += realSCSNum;
		//string scsWithWildcards = SCS2WithWildCardsRandom(copies[0], copies[1],generator);
		string scsWithWildcards = SCS2WithWildCards(copies[0], copies[1]);
		int estimatedSCSNum = EstimateSCSNum(scsWithWildcards);
		cumEstimatedNum += estimatedSCSNum;
//		if (realSCSNum != estimatedSCSNum){
//			cout<<copies[0]<<endl;
//			cout<<copies[1]<<endl;
//			cout<<scsWithWildcards<<endl;
//			cout<<estimatedSCSNum<<endl;
//			cout<<scss;
//			break;
//		}
	}
	cout << "Real SCS Num Avg:\t" << cumRealNum / reps << endl;
	cout << "Estimated SCS Num Avg:\t" << cumEstimatedNum / reps << endl;
}

void TestMinEditDistanceFilteredSCS(const int reps, const int originalLen, const int copiesNum, const int k,
		const double delProb, mt19937& generator) {
	vector<string> copies;
	double cumEditDistance = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		string mostLikelySCS = MinEditDistanceFilteredSCS(copies, k, delProb, originalLen);
		cumEditDistance += EditDistanceDI(original, mostLikelySCS);
	}
	cout << "Min Edit Distance Filtered SCS Letter Error Rate:\t" << cumEditDistance / (reps * originalLen) << endl;
}

void TestEitan(const int reps, const int originalLen, const double delProb, mt19937& generator) {
	vector<string> copies;
	double cumSCSNum = 0, cumMostLikelyNum = 0;
	int currentSCSNum = 0, currentMostLikelyNum = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, 2, delProb, generator);
		MostLikelyEitan(copies, originalLen, currentSCSNum, currentMostLikelyNum);
		cumSCSNum += currentSCSNum;
		cumMostLikelyNum += currentMostLikelyNum;
	}
	cout << "Original Len:\t" << originalLen << endl;
	cout << "Deletion Probability:\t" << delProb << endl;
	cout << "Average Number of SCSs:\t\t\t" << cumSCSNum / reps << endl;
	cout << "Average Number of Most Likely SCSs:\t" << cumMostLikelyNum / reps << endl;
}

// Error Vector to transform X to Y
string ErrorVector(const string& X, const string& Y) {
	vector<vector<pair<int, int>>> lcs2IndexesAll = LCS2Indexes(X, Y);
	vector<pair<int, int>> lcs2Indexes = lcs2IndexesAll[0];
	string errorVector;
	int lcsLen = lcs2Indexes.size();
	int xIndex = 0, yIndex = 0, currentLCSXindex, currentLCSYindex;
	for (int i = 0; i < lcsLen; i++) {
		currentLCSXindex = lcs2Indexes[i].first;
		currentLCSYindex = lcs2Indexes[i].second;
		while (xIndex < currentLCSXindex) {
			errorVector += '1';
			xIndex++;
		}
		while (yIndex < currentLCSYindex) {
			errorVector += '2';
			yIndex++;
		}
		assert(xIndex == currentLCSXindex && yIndex == currentLCSYindex);
		errorVector += '0';
		xIndex++;
		yIndex++;
	}
	// after last lcs letter
	currentLCSXindex = X.size();
	while (xIndex < currentLCSXindex) {
		errorVector += '1';
		xIndex++;
	}
	currentLCSYindex = Y.size();
	while (yIndex < currentLCSYindex) {
		errorVector += '2';
		yIndex++;
	}
	return errorVector;
}

void TestEitanOutput(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator, const int maxOutputNum) {
	vector<string> copies;
	int outputCount = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		vector<string> allMostLikelySCS = AllMostLikelyEitan(copies, originalLen);
		string mostLikelySCS = allMostLikelySCS[0];
		if (mostLikelySCS != original) {
			cout << "original:\t" << original << endl;
			cout << endl;
			for (int i = 0; i < copiesNum; i++) {
				cout << "copy" << i << "\t\t" << copies[i] << endl;
			}
			cout << endl;
			for (unsigned i = 0; i < allMostLikelySCS.size(); i++) {
				cout << "ML" << i << "\t\t" << allMostLikelySCS[i] << endl;
			}
			cout << endl;
			cout << "chosen ML:\t" << mostLikelySCS << endl;
			string errorVector = ErrorVector(mostLikelySCS, original);
			cout << "errors:\t\t" << errorVector << endl;
			cout << "original:\t" << original << endl;
			cout << endl;
			cout << "deletions('1'):\t\t" << count(errorVector.begin(), errorVector.end(), '1') << endl;
			cout << "insertions('2'):\t" << count(errorVector.begin(), errorVector.end(), '2') << endl;
			outputCount++;
			if (outputCount == maxOutputNum) {
				return;
			}
			cout << endl;
			cout << "*********************************************************************" << endl;
		}

	}
}

void TestEitanFull(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	vector<string> copies;
	map<int, int> SCSofLenCount, cumSCSNum, cumMostLikelyNum, cumEditDistance;
	int currentSCSNum = 0, currentMostLikelyNum = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		//string original = RandomCopies2(copies, originalLen, copiesNum, delProb, generator);
		//cout<<"original:\t"<<original<<endl;
		string mostLikelySCS = MostLikelyEitanFull(copies, originalLen, currentSCSNum, currentMostLikelyNum);
		//cout<<"most likely:\t"<<mostLikelySCS<<endl;
		int scsLen = mostLikelySCS.size();
		SCSofLenCount[scsLen]++;
		cumSCSNum[scsLen] += currentSCSNum;
		cumMostLikelyNum[scsLen] += currentMostLikelyNum;
		cumEditDistance[scsLen] += EditDistanceDI(original, mostLikelySCS);
	}
	cout << "Copies Num:\t\t" << copiesNum << endl;
	cout << "Original Len:\t\t" << originalLen << endl;
	cout << "Deletion Probability:\t\t" << delProb << endl;
	cout << "Repetitions:\t\t" << reps << endl;
	cout << "SCS Len\tOccurences\tSCS Total\tML SCS Total\tTotal ED" << endl;
	for (map<int, int>::reverse_iterator rit = SCSofLenCount.rbegin(); rit != SCSofLenCount.rend(); rit++) {
		int scsLen = rit->first;
		cout << scsLen << "\t" << SCSofLenCount[scsLen] << "\t" << cumSCSNum[scsLen] << "\t" << cumMostLikelyNum[scsLen]
				<< "\t" << cumEditDistance[scsLen] << endl;
	}
}

void TestFrequencies(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	vector<string> copies;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		//cout << copies;
		//cout << SCSNLenFast(copies,originalLen)<<endl;
		vector<map<char, int>> result = MaxLetterFrequency(copies, originalLen);
		for (int i = 0; i < originalLen; i++) {
			cout << "A " << result[i]['A'] << "\t";
			cout << "C " << result[i]['C'] << "\t";
			cout << "G " << result[i]['G'] << "\t";
			cout << "T " << result[i]['T'] << "\t";
			cout << "index: " << i + 1 << "\t";
			int sum = result[i]['A'] + result[i]['C'] + result[i]['G'] + result[i]['T'];
			cout << "Sum: " << sum << "\tDiff: " << i + 1 - sum << endl;
		}
		cout << "**********************************************************" << endl;
		vector<string> reverseCopies(copies);
		for (int i = 0; i < copiesNum; i++) {
			reverse(reverseCopies[i].begin(), reverseCopies[i].end());
		}
		result = MaxLetterFrequency(reverseCopies, originalLen);
		for (int i = 0; i < originalLen; i++) {
			cout << "A " << result[i]['A'] << "\t";
			cout << "C " << result[i]['C'] << "\t";
			cout << "G " << result[i]['G'] << "\t";
			cout << "T " << result[i]['T'] << "\t";
			cout << "index: " << i + 1 << "\t";
			int sum = result[i]['A'] + result[i]['C'] + result[i]['G'] + result[i]['T'];
			cout << "Sum: " << sum << "\tDiff: " << i + 1 - sum << endl;
		}
	}
}

bool ThereIsCorrrectSizePair(const vector<string>& copies, const int originalLen) {
	int copiesNum = copies.size();
	int tempSCSLen;
	for (int i = 0; i < copiesNum; i++) {
		for (int j = i + 1; j < copiesNum; j++) {
			tempSCSLen = SCS2LenFast(copies[i], copies[j], originalLen);
			if (tempSCSLen == originalLen) {
				return true;
			}
		}
	}
	return false;
}

void TestCountPairwiseCorrectSize(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	vector<string> copies;
	int countShort = 0;
	for (int counter = 0; counter < reps; counter++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		if (not ThereIsCorrrectSizePair(copies, originalLen)) {
			countShort++;
		}
	}
	cout << "Too short count:\t" << countShort << endl;
}

struct LetterStats {
	int aNum;
	int aRunNum;
	int cNum;
	int cRunNum;
	int gNum;
	int gRunNum;
	int tNum;
	int tRunNum;
	int SumRun() {
		return aRunNum + cRunNum + gRunNum + tRunNum;
	}
};

LetterStats CountLetters(const string& str) {
	int strLen = str.size();
	LetterStats stats;
	map<char, int> count;
	map<char, int> runCount;
	char lastLetter = 'X', currentLetter;
	for (int i = 0; i < strLen; i++) {
		currentLetter = str[i];
		count[currentLetter]++;
		if (currentLetter != lastLetter) {
			runCount[currentLetter]++;
		}
		lastLetter = currentLetter;
	}
	stats.aNum = count['A'];
	stats.aRunNum = runCount['A'];
	stats.cNum = count['C'];
	stats.cRunNum = runCount['C'];
	stats.gNum = count['G'];
	stats.gRunNum = runCount['G'];
	stats.tNum = count['T'];
	stats.tRunNum = runCount['T'];
	return stats;
}

int ChooseValue(int compA, int compB, int valueA, int valueB) {
	if (compA > compB) {
		return valueA - valueB;
	}
	else if (compA == compB) {
		return abs(valueA - valueB);
	}
	else {
		return valueB - valueA;
	}
}

// sum max:									28801
// sum max break ties by sum len:			20616
// sum max break ties by max total run num:	32197
// sum max break ties by min total run num:	30649
// sum max break ties by min run num diff:	28803
// sum max pairwise run:					373830
// sum max break ties first: by sum len. break ties second: by min sum run 18874

// TODO: count number of runs for each letter/string - the more the better? or do i need many in one and few in other?
void TestSCSLen(const int reps, const int originalLen, const int copiesNum, const double delProb, mt19937& generator) {
	vector<string> copies;
	int countShort = 0;
	for (int counter = 0; counter < reps; counter++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		vector<LetterStats> stats(copiesNum);
		for (int i = 0; i < copiesNum; i++) {
			stats[i] = CountLetters(copies[i]);
		}

		//int firstIndex=0,secondIndex=0,sumMax;
		int maxA, maxC, maxG, maxT, sumMax, firstIndex = 0, secondIndex = 0; //tieBreaker1 = 0, tieBreaker2 = 0;
		int maxSumMax = 0;

		// find pair with max sum
		for (int i = 0; i < copiesNum; i++) {
			for (int j = i + 1; j < copiesNum; j++) {
				maxA = max(stats[i].aNum, stats[j].aNum);
				maxC = max(stats[i].cNum, stats[j].cNum);
				maxG = max(stats[i].gNum, stats[j].gNum);
				maxT = max(stats[i].tNum, stats[j].tNum);
				sumMax = maxA + maxC + maxG + maxT;

				if (sumMax > maxSumMax) {
					maxSumMax = sumMax;
					firstIndex = i;
					secondIndex = j;
//					tieBreaker1 = copies[i].size() + copies[j].size();
//					tieBreaker2 = stats[i].SumRun() + stats[j].SumRun();
				}
//				else if (sumMax == maxSumMax) {
//					int temp1 = copies[i].size() + copies[j].size();
//					int temp2 = stats[i].SumRun() + stats[j].SumRun();
//					if (temp1 > tieBreaker1) {
//						firstIndex = i;
//						secondIndex = j;
//						tieBreaker1 = temp1;
//						tieBreaker2 = temp2;
//					}
//					else if (temp1 == tieBreaker1) {
//						if (temp2 < tieBreaker2) {
//							firstIndex = i;
//							secondIndex = j;
//							tieBreaker2 = temp2;
//						}
//					}
//				}
			}
		}
		int chosenSCSLen = SCS2LenFast(copies[firstIndex], copies[secondIndex], originalLen);

		if (chosenSCSLen < originalLen) {
			countShort++;
		}
	}
	cout << "Max Short Num: " << countShort << endl;
}

// create all words of len K from alphabet GTAC
void AllKLenWords(int K) {
	vector<char> alphabet = { 'G', 'T', 'A', 'C' };
	vector<int> dimSizes(K);
	for (int i = 0; i < K; i++) {
		dimSizes[i] = 4;
	}
	IndexVector indexVector(dimSizes);
	int j = 0;
	for (IndexVector currentIndex = indexVector; not currentIndex.IsPastEnd(); currentIndex++) {
		cout << Word(currentIndex.Vector(), alphabet) << endl;
		j++;
	}
	cout << "string num: " << j << endl;
}

// Brute force find a common supersequence of len K

string KLenOneCommonSupersequence(const vector<string>& strings, const int K) {
	if (K < 0) {
		cout << "Invalid length!" << endl;
		return "";
	}
	if (K == 0) {	//0 len string (the empty string) can be supersequence only of empty strings.
		for (unsigned i = 0; i < strings.size(); i++) {
			if (not strings[i].empty()) {
				return "X";
			}
		}
		return "";
	}

	vector<char> alphabet = { 'G', 'T', 'A', 'C' };
	vector<int> dimSizes(K);
	for (int i = 0; i < K; i++) {
		dimSizes[i] = 4;
	}
	IndexVector indexVector(dimSizes);
	for (IndexVector currentIndex = indexVector; not currentIndex.IsPastEnd(); currentIndex++) {
		string word = Word(currentIndex.Vector(), alphabet);
		bool isSupersequence = true;
		for (unsigned j = 0; j < strings.size(); j++) {
			if (not IsSubSequence(strings[j], word)) {
				isSupersequence = false;
				break;
			}
		}
		if (isSupersequence) {
			return word;
		}
	}
	return "X";
}

// Brute force find all common supersequences of len K

vector<string> KLenAllCommonSupersequences(const vector<string>& strings, const int K) {
	vector<string> result;
	if (K < 0) {
		cout << "Invalid length!" << endl;
		return result;
	}
	if (K == 0) {	//0 len string (the empty string) can be supersequence only of empty strings.
		for (unsigned i = 0; i < strings.size(); i++) {
			if (not strings[i].empty()) {
				return result; // empty string is not common supersequence, return empty vector.
			}
		}
		result.push_back(""); // empty string is common supersequence.
		return result;
	}

	vector<char> alphabet = { 'G', 'T', 'A', 'C' };
	vector<int> dimSizes(K);
	for (int i = 0; i < K; i++) {
		dimSizes[i] = 4;
	}
	IndexVector indexVector(dimSizes);
	for (IndexVector currentIndex = indexVector; not currentIndex.IsPastEnd(); currentIndex++) {
		string word = Word(currentIndex.Vector(), alphabet);
		bool isSupersequence = true;
		for (unsigned j = 0; j < strings.size(); j++) {
			if (not IsSubSequence(strings[j], word)) {
				isSupersequence = false;
				break;
			}
		}
		if (isSupersequence) {
			result.push_back(word);
		}
	}
	return result;
}

// check SCSN and SCS2 is the same for two strings
// generate two random strings X,Y, which will be input strings for SCS2
// for SCSN input strings are X,Y and 4 empty strings
// check equality of SCSs

void Test1(int reps, mt19937& generator, const int minStrLen, const int maxStrLen) {
	uniform_int_distribution<int> lenDistribution(minStrLen, maxStrLen);
	vector<int> indexes = { 0, 1, 2, 3, 4, 5 };
	for (int i = 0; i < reps; i++) {
		vector<string> scsnInput(6);
		int len1 = lenDistribution(generator);
		int len2 = lenDistribution(generator);
		string X = MakeRandomString4(len1, generator);
		string Y = MakeRandomString4(len2, generator);

		shuffle(indexes.begin(), indexes.end(), generator);
		scsnInput[indexes[0]] = X;
		scsnInput[indexes[1]] = Y;

		vector<string> SCS = SCSN(scsnInput);
		vector<string> scs2 = SCS2(X, Y);

		set<string> SCS2Set(scs2.begin(), scs2.end());
		set<string> SCSNSet(SCS.begin(), SCS.end());

		assert(SCS2Set.size() == SCSNSet.size());

		if (SCS2Set != SCSNSet) {
			cout << "Sets NOT Equal!" << endl;
			cout << "SCSN input strings:" << endl << scsnInput;
			cout << "SCS2 input strings:" << endl;
			cout << X << endl;
			cout << Y << endl;
			cout << "SCSN/SCS2 Test Failure" << endl;
			return;
		}
	}
	cout << "SCSN/SCS2 Test Success!" << endl;
}

// check SCSN and SCS2 is the same for two strings
// generate two random strings X,Y, which will be input strings for SCS2
// for SCSN input strings are X,Y 2 random subsequences of X,2 random subsequences of Y
// check equality of SCSs

void Test2(int reps, mt19937& generator, const int minStrLen, const int maxStrLen) {
	uniform_int_distribution<int> lenDistribution(minStrLen, maxStrLen);
	uniform_real_distribution<double> delDistribution(0.0, 1.0);
	vector<int> indexes = { 0, 1, 2, 3, 4, 5 };
	for (int i = 0; i < reps; i++) {
		vector<string> scsnInput(6);
		int len1 = lenDistribution(generator);
		int len2 = lenDistribution(generator);
		string X = MakeRandomString4(len1, generator);
		string Y = MakeRandomString4(len2, generator);
		shuffle(indexes.begin(), indexes.end(), generator);
		scsnInput[indexes[0]] = X;
		scsnInput[indexes[1]] = Y;
		scsnInput[indexes[2]] = CopyStrandWithDeletion(X, generator, delDistribution(generator));
		scsnInput[indexes[3]] = CopyStrandWithDeletion(X, generator, delDistribution(generator));
		scsnInput[indexes[4]] = CopyStrandWithDeletion(Y, generator, delDistribution(generator));
		scsnInput[indexes[5]] = CopyStrandWithDeletion(Y, generator, delDistribution(generator));

		vector<string> SCS = SCSN(scsnInput);
		vector<string> scs2 = SCS2(X, Y);

		set<string> SCS2Set(scs2.begin(), scs2.end());
		set<string> SCSNSet(SCS.begin(), SCS.end());

		assert(SCS2Set.size() == SCSNSet.size());

		if (SCS2Set != SCSNSet) {
			cout << "Sets NOT Equal!" << endl;
			cout << "SCSN input strings:" << endl << scsnInput;
			cout << "SCS2 input strings:" << endl;
			cout << X << endl;
			cout << Y << endl;
			cout << "SCSN/SCS2 Test Failure" << endl;
			return;
		}
	}
	cout << "SCSN/SCS2 Test Success!" << endl;
}

void RandomStrings(vector<string>& strings, uniform_int_distribution<int>& strNumDistribution,
		uniform_int_distribution<int>& strLenDistribution, mt19937& generator) {
	int stringNum = strNumDistribution(generator);
	strings = vector<string>(stringNum);
	for (int i = 0; i < stringNum; i++) {
		strings[i] = MakeRandomString4(strLenDistribution(generator), generator);
	}
}

// TEST SCSN IS SUPERSEQUENCE
// generate random number (2 to 6) of random strings of random len(0 to 10)
// test if all SCS are supersequences of all input strings

void TestSCSNIsSuperseq(int reps, mt19937& generator, const int minStrNum, const int maxStrNum, const int minStrLen,
		const int maxStrLen) {
	uniform_int_distribution<int> strNumDist(minStrNum, maxStrNum);
	uniform_int_distribution<int> strLenDist(minStrLen, maxStrLen);
	vector<string> scsnInput;
	for (int k = 0; k < reps; k++) {
		RandomStrings(scsnInput, strNumDist, strLenDist, generator);

		vector<string> SCS = SCSN(scsnInput);
		for (unsigned i = 0; i < scsnInput.size(); i++) {
			for (unsigned j = 0; j < SCS.size(); j++) {
				if (not IsSubSequence(scsnInput[i], SCS[j])) {
					cout << "Supersequence Test Failure" << endl;
					cout << "SCSN input strings:" << endl << scsnInput;
					return;
				}

			}
		}
	}
	cout << "Supersequence Test Success!" << endl;
}

// TEST SCSN INPUT ORDER
// generate random number (2 to 6) of random strings of random len(0 to 10)
// compute SCSN
// shuffle input strings.
// Test SCSN hasn't changed

void TestSCSNInputOrder(int reps, mt19937& generator, const int minStrNum, const int maxStrNum, const int minStrLen,
		const int maxStrLen) {
	uniform_int_distribution<int> strNumDist(minStrNum, maxStrNum);
	uniform_int_distribution<int> strLenDist(minStrLen, maxStrLen);
	vector<string> scsnInput;
	for (int k = 0; k < reps; k++) {
		RandomStrings(scsnInput, strNumDist, strLenDist, generator);

		vector<string> SCS1 = SCSN(scsnInput);
		set<string> SCS1Set(SCS1.begin(), SCS1.end());
		vector<string> suffledscsnInput = scsnInput;
		shuffle(suffledscsnInput.begin(), suffledscsnInput.end(), generator);
		vector<string> scs2 = SCSN(suffledscsnInput);
		set<string> SCS2Set(scs2.begin(), scs2.end());

		if (SCS1Set != SCS2Set) {
			cout << "SCSN Input Order Test Failure" << endl;
			cout << "SCSN input strings:" << endl << scsnInput;
			cout << "Shuffled SCSN input strings:" << endl << suffledscsnInput;
			return;
		}
	}
	cout << "SCSN Input Order Test Success!" << endl;
}

// TEST SCSN HAS ALL SCSS
// generate random number (2 to 6) of random strings of random len(0 to 8)
// compute SCSN
// brute force to find all SCSs of same len as SCSN
// compare if the two sets are equal

void TestSCSNHasAll(int reps, mt19937& generator, const int minStrNum, const int maxStrNum, const int minStrLen,
		const int maxStrLen, const int maxSCSLenToProcess) {

	uniform_int_distribution<int> strNumDist(minStrNum, maxStrNum);
	uniform_int_distribution<int> strLenDist(minStrLen, maxStrLen);
	vector<string> scsnInput;
	for (int k = 0; k < reps; k++) {
		RandomStrings(scsnInput, strNumDist, strLenDist, generator);

		int scsLen = SCSNLen(scsnInput);

		if (scsLen > maxSCSLenToProcess) {
			continue;
		}

		vector<string> SCS = SCSN(scsnInput);
		set<string> SCSset(SCS.begin(), SCS.end());

		vector<string> brute = KLenAllCommonSupersequences(scsnInput, scsLen);
		set<string> BruteSet(brute.begin(), brute.end());

		if (SCSset != BruteSet) {
			cout << "SCSN Has All SCSs Test Failure" << endl;
			cout << "SCSN input strings:" << endl << scsnInput;
			return;
		}
	}
	cout << "SCSN Has All SCSs Test Success!" << endl;
}

// SCSN IS SHORTEST TEST
// generate random number (2 to 6) of random strings of random len(0 to 8)
// compute SCSN len
// brute force to find SCS of len-1, i.e. a shorter SCS

void TestSCSNIsShortest(int reps, mt19937& generator, const int minStrNum, const int maxStrNum, const int minStrLen,
		const int maxStrLen, const int maxSCSLenToProcess) {
	uniform_int_distribution<int> strNumDist(minStrNum, maxStrNum);
	uniform_int_distribution<int> strLenDist(minStrLen, maxStrLen);
	vector<string> scsnInput;
	for (int k = 0; k < reps; k++) {
		RandomStrings(scsnInput, strNumDist, strLenDist, generator);

		int scsLen = SCSNLen(scsnInput);
		if (scsLen == 0) {
			continue;
		}
		if (scsLen > maxSCSLenToProcess) {
			continue;
		}

		string brute = KLenOneCommonSupersequence(scsnInput, scsLen - 1);
		if (brute != "X") {
			cout << "SCSN Is Shortest Test Failure" << endl;
			cout << "SCSN input strings:" << endl << scsnInput;
		}
	}
	cout << "SCSN Is Shortest Test Success" << endl;
}

bool CompareLenDes(const string& a, const string& b) {
	return b.size() < a.size();
}

bool CompareEditDistance(const pair<string, string>& a, const pair<string, string>& b) {
	return EditDistanceDI(a.first, a.second) < EditDistanceDI(b.first, b.second);
}

pair<string, string> MaxEditDistancePair(const vector<pair<string, string> >& pairs) {
	int pairNum = pairs.size();
	vector<pair<string, string> > temp(pairs);
	sort(temp.begin(), temp.end(), CompareEditDistance);
	return temp[pairNum - 1];
}

// return vector with K longest strings from a
vector<string> ChooseKLongest(const vector<string>& a, const int K) {
	vector<string> temp(a), result;
	sort(temp.begin(), temp.end(), CompareLenDes);
	for (int i = 0; i < K; i++) {
		result.push_back(temp[i]);
	}
	return result;
}

// test most likely scs of max edit distance pair
void TestMostLikelySCSMaxEditDistance(const int reps, const int originalLen, const int copiesNum, const double delProb,
		mt19937& generator) {
	vector<string> copies;
	double cumEditDistance = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		vector<pair<string, string> > pairs = StringPairings(copies);
		pair<string, string> bestPair = MaxEditDistancePair(pairs);
		vector<string> pairVec = { bestPair.first, bestPair.second };
		string mostLikelySCS = MostLikelySCS(pairVec, delProb, originalLen, generator);
		cumEditDistance += EditDistanceDI(original, mostLikelySCS);
	}
	cout << "Letter Error Rate:\t" << cumEditDistance / (reps * originalLen) << endl;
}

// test most likely scs of k longest input strings
void TestMostLikelySCSKLongest(const int reps, const int originalLen, const int copiesNum, const int K,
		const double delProb, mt19937& generator) {
	vector<string> copies, BestK;
	double cumEditDistance = 0;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		BestK = ChooseKLongest(copies, K);
		string mostLikelySCS = MostLikelySCS(BestK, delProb, originalLen, generator);
		cumEditDistance += EditDistanceDI(original, mostLikelySCS);
	}
	cout << "Letter Error Rate:\t" << cumEditDistance / (reps * originalLen) << endl;
}

void TimeSCS2Len(const int reps, const int originalLen, const double delProb, mt19937& generator) {
	string original, X, Y;
	for (int i = 0; i < reps; i++) {
		original = MakeRandomString4(originalLen, generator);
		X = CopyStrandWithDeletion(original, generator, delProb);
		Y = CopyStrandWithDeletion(original, generator, delProb);
		SCS2Len(X, Y);
	}
}
void TimeSCS2(const int reps, const int originalLen, const double delProb, mt19937& generator) {
	string original, X, Y;
	for (int i = 0; i < reps; i++) {
		original = MakeRandomString4(originalLen, generator);
		X = CopyStrandWithDeletion(original, generator, delProb);
		Y = CopyStrandWithDeletion(original, generator, delProb);
		vector<string> v = SCS2(X, Y);
	}
}


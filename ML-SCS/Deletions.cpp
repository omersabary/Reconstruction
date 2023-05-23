#include <iostream>
#include <set>
#include <random>
#include <chrono>
#include <cassert>
#include <algorithm>
#include "Tests.hpp"
#include "Strings.hpp"
#include "SCS2.hpp"
#include "SCSN.hpp"
#include "NChooseK.hpp"
#include "MultiString.hpp"

using namespace std;

void TestLCS(const int reps, const int originalLen, const double delProb, mt19937& generator) {
	double totalLCSNum = 0;
	for (int i = 0; i < reps; i++) {
		string original, X, Y;
		original = MakeRandomString4(originalLen, generator);
		X = CopyStrandWithDeletion(original, generator, delProb);
		Y = CopyStrandWithDeletion(original, generator, delProb);
		vector<vector<pair<int, int>>> test = LCS2Indexes(X,Y);
		set<vector<pair<int, int>>> st(test.begin(),test.end());
		totalLCSNum += st.size();
	}
	cout << "Avg. LCS Num:\t" << totalLCSNum / reps << endl;
}

void TestSCSN(int reps, const int originalLen, const int copiesNum, const double delProb, mt19937& generator) {
	vector<string> copies;
	for (int i = 0; i < reps; i++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
		//SCSNLen(copies);
		SCS2Fastest(copies[0], copies[1], originalLen);
		//SCSNFast(copies,originalLen);
	}
}

int max(const int a, const int b, const int c) {
	int max;
	max = (a > b) ? a : b;
	max = (c > max) ? c : max;
	return max;
}

int min(const int a, const int b, const int c, const int d) {
	return min(min(a, b, c), d);
}

int max(const int a, const int b, const int c, const int d) {
	return max(max(a, b, c), d);
}

void TestFastSCS2(const int reps, const int originalLen, const double delProb, mt19937& generator) {
	string original, X, Y;
	for (int i = 0; i < reps; i++) {
		original = MakeRandomString4(originalLen, generator);
		X = CopyStrandWithDeletion(original, generator, delProb);
		Y = CopyStrandWithDeletion(original, generator, delProb);

		vector<string> scs2 = SCS2(X, Y);
		vector<string> scs2fast = SCS2Fast(X, Y, originalLen);
		vector<string> scs2fastest = SCS2Fastest(X, Y, originalLen);
		bool successStrings = (scs2 == scs2fast and scs2 == scs2fastest);

		int scs2len = SCS2Len(X, Y);
		int scs2lenfast = SCS2LenFast(X, Y, originalLen);
		bool successLengths = (scs2len == scs2lenfast and scs2len == (int) scs2[0].size());

		if ((not successStrings) or (not successLengths)) {
			cout << "Fast SCS2 Test Failed!" << endl;
			cout << X << endl;
			cout << Y << endl;
			return;
		}
	}
	cout << "Fast SCS2 Test Success!" << endl;
}

void TestSCSNFast(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator) {
	vector<string> copies;
	for (int k = 0; k < reps; k++) {
		string original = RandomCopies4(copies, originalLen, copiesNum, delProb, generator);
//		SCSNLen(copies);
//		SCSNLenFast(copies, original.size());
//		SCSNLenFastLessMem(copies, original.size());

//		vector<string> scsn = SCSN(copies);
		SCSNFast(copies, original.size());
//		int scsLen = SCSNLen(copies);
//		int scsnLenFast = SCSNLenFast(copies, original.size());
//
//		if (scsn != scsnfast) {
//			cout << "SCSNFast same as SCSN Test Failure!" << endl;
//			cout << "Input strings:" << endl << copies;
//			return;
//		}
	}
//	cout << "SCSNFast same as SCSN Test Success" << endl;
}

// TODO: check SCS len as a function of number of substitutions in Edit distance
// 		 if number of substitutions is zero can we lose a letter other than those lost in both strings?
//       consider also edit distance with transpositions. does it help?

// TODO: best choice of a pair among pairs with correct SCS len
// TODO: find a clever way to transform a too short SCS to supersequence of all
// TODO: tests of len 200 of algorithm. 5 copies.
// Algorithm: for all possible pairs find a pair whose SCS is of correct len.
//				if found:	1. filter SCS which are not supersequence of all copies
//							2. return most likely SCS. first, if more than 1.
//				if not found:	1. return most likely SCS, if score is not zero.
//								2. if score is zero, choose scs with min total Levinshtein distance
int main() {
	clock_t begin = clock();

	//unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	unsigned sd = 3865506088;
	mt19937 generator(sd);

	int reps = 100;
	int originalLen = 100;
	double delProb = 0.10;
	int copiesNum = 10;
	int N = 4;

//	TestBigClusterAlgorithm(reps, originalLen, copiesNum, delProb, k, generator);

// For OMER ************************************************************************************
//	void TestAlgorithmNSteps(const int reps, const int maxStageNum, const int originalLen, const int copiesNum,
//	const double delProb, mt19937& generator);

//  reps - number of iterations
// 	maxStageNum - the highest the function gets with ktupling - should be 4 - i.e. reaches at most SCS of 4.
// 	originalLen - length of original strand
//	copiesNum - number of copies
//  delProb - deletion probability
//	generator - for the random function. NOTE: currently seeds each time with same number. see line 120

	TestAlgorithmNSteps(reps, N, originalLen, copiesNum, delProb, generator);

//	for (int copiesNum = 4; copiesNum <= 10; copiesNum++) {
//		mt19937 generator(sd);
//		TestAlgorithmNSteps(reps, N, originalLen, copiesNum, delProb, generator);
//	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << endl << "Time elapsed: " << elapsed_secs << "\tseconds" << endl;
	return 0;
}

//vector<vector<vector<int>>> test(x, vector<vector<int>>(y, vector<int>(z, value)));

//vector<string> input= {"A","BB","CCC","D","EEEEE"};
//	int k=3;
//	vector<vector<string>> stringktuples= SortedStringKtuples(input,k);
//	for (unsigned i=0;i<stringktuples.size();i++){
//		for (int j=0;j<k;j++){
//			cout<< stringktuples[i][j]<<"\t";
//		}
//		cout << endl;
//	}

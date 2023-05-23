#include "EditDistance.hpp"
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <algorithm>
using namespace std;

int EditDistanceArray(const string& X, const string& Y, int m, int n, vector<vector<int> >& dp) {
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			if (i == 0) {
				dp[i][j] = j;
			}
			else if (j == 0) {
				dp[i][j] = i;
			}
			else if (X[i - 1] == Y[j - 1]) {
				dp[i][j] = dp[i - 1][j - 1];
			}
			else {
				dp[i][j] = 1 + min(min(dp[i][j - 1], dp[i - 1][j]), dp[i - 1][j - 1]);
			}
		}
	}
	return dp[m][n];
}

//vector<LetterOps> BacktrackEditDistanceOperations(const string& X, const string& Y, int m, int n,
//		const vector<vector<int> >& dp) {
//	vector<LetterOps> s(X.size() + 1);
//	if (m == 0) {
//		s[0].insert = Y.substr(0, n);
//		return s;
//	}
//	else if (n == 0) {
//		for (int index = 0; index < m; index++) {
//			s[index].CDR = "D";
//		}
//		return s;
//	}
//	if (X[m - 1] == Y[n - 1]) {
//		vector<LetterOps> s = BacktrackEditDistanceOperations(X, Y, m - 1, n - 1, dp);
//		s[m - 1].CDR = "C";
//		s[m - 1].opParam = string(1, X[m - 1]);
//		return s;
//	}
//	else {
//
//		if (dp[m][n] == dp[m - 1][n - 1] + 1) {  	// replace
//			vector<LetterOps> s1 = BacktrackEditDistanceOperations(X, Y, m - 1, n - 1, dp);
//			s1[m - 1].CDR = "R";
//			s1[m - 1].opParam = string(1, Y[n - 1]);
//			return s1;
//		}
//
//		if (dp[m][n] == dp[m][n - 1] + 1) { 		// insert
//			vector<LetterOps> s2 = BacktrackEditDistanceOperations(X, Y, m, n - 1, dp);
//			s2[m].insert += Y[n - 1];
//			return s2;
//		}
//		if (dp[m][n] == dp[m - 1][n] + 1) {			// delete
//			vector<LetterOps> s3 = BacktrackEditDistanceOperations(X, Y, m - 1, n, dp);
//			s3[m - 1].CDR = "D";
//			return s3;
//		}
//		assert(0);
//	}
//}

// subPriority, delPriority, insPriority - permutation of 1,2,3
vector<LetterOps> BacktrackEditDistancePriority(const string& X, const string& Y, int m, int n,
		const vector<vector<int> >& dp, const int subPriority, const int insPriority, const int delPriority) {
	vector<LetterOps> s(X.size() + 1);
	if (m == 0) {
		s[0].insert = Y.substr(0, n);
		return s;
	}
	else if (n == 0) {
		for (int index = 0; index < m; index++) {
			s[index].CDR = "D";
		}
		return s;
	}
	if (X[m - 1] == Y[n - 1]) {
		vector<LetterOps> s = BacktrackEditDistancePriority(X, Y, m - 1, n - 1, dp, subPriority, insPriority,
				delPriority);
		s[m - 1].CDR = "C";
		s[m - 1].opParam = string(1, X[m - 1]);
		return s;
	}
	else {
		int subEqual = 0, delEqual = 0, insEqual = 0;
		if (dp[m][n] == dp[m - 1][n - 1] + 1) {
			subEqual = 1;
		}
		if (dp[m][n] == dp[m][n - 1] + 1) {
			insEqual = 1;
		}
		if (dp[m][n] == dp[m - 1][n] + 1) {
			delEqual = 1;
		}
		int subScore = subEqual * subPriority, delScore = delEqual * delPriority, insScore = insEqual * insPriority;
		bool chooseSub = false, chooseDel = false, chooseIns = false;
		if (subScore > delScore and subScore > insScore) {
			chooseSub = true;
		}
		else if (delScore > subScore and delScore > insScore) {
			chooseDel = true;
		}
		else {
			assert(insScore > subScore and insScore > delScore);
			chooseIns = true;
		}
		if (chooseSub) {  	// replace
			vector<LetterOps> s1 = BacktrackEditDistancePriority(X, Y, m - 1, n - 1, dp, subPriority, insPriority,
					delPriority);
			s1[m - 1].CDR = "R";
			s1[m - 1].opParam = string(1, Y[n - 1]);
			return s1;
		}

		else if (chooseIns) { 		// insert
			vector<LetterOps> s2 = BacktrackEditDistancePriority(X, Y, m, n - 1, dp, subPriority, insPriority,
					delPriority);
			s2[m].insert += Y[n - 1];
			return s2;
		}
		else {			// delete
			assert(chooseDel);
			vector<LetterOps> s3 = BacktrackEditDistancePriority(X, Y, m - 1, n, dp, subPriority, insPriority,
					delPriority);
			s3[m - 1].CDR = "D";
			return s3;
		}
	}
}

vector<LetterOps> BacktrackEditDistanceRandom(const string& X, const string& Y, int m, int n,
		const vector<vector<int> >& dp, vector<int>& priorities, mt19937& generator) {
	vector<LetterOps> s(X.size() + 1);
	if (m == 0) {
		s[0].insert = Y.substr(0, n);
		return s;
	}
	else if (n == 0) {
		for (int index = 0; index < m; index++) {
			s[index].CDR = "D";
		}
		return s;
	}
	if (X[m - 1] == Y[n - 1]) {
		vector<LetterOps> s = BacktrackEditDistanceRandom(X, Y, m - 1, n - 1, dp, priorities, generator);
		s[m - 1].CDR = "C";
		s[m - 1].opParam = string(1, X[m - 1]);
		return s;
	}
	else {
		//random_shuffle(priorities.begin(),priorities.end());
		shuffle(priorities.begin(), priorities.end(), generator);
		int subEqual = 0, delEqual = 0, insEqual = 0;
		if (dp[m][n] == dp[m - 1][n - 1] + 1) {
			subEqual = 1;
		}
		if (dp[m][n] == dp[m][n - 1] + 1) {
			insEqual = 1;
		}
		if (dp[m][n] == dp[m - 1][n] + 1) {
			delEqual = 1;
		}
		int subScore = subEqual * priorities[0], delScore = delEqual * priorities[2], insScore = insEqual
				* priorities[1];
		bool chooseSub = false, chooseDel = false, chooseIns = false;
		if (subScore > delScore and subScore > insScore) {
			chooseSub = true;
		}
		else if (delScore > subScore and delScore > insScore) {
			chooseDel = true;
		}
		else {
			assert(insScore > subScore and insScore > delScore);
			chooseIns = true;
		}
		if (chooseSub) {  	// replace
			vector<LetterOps> s1 = BacktrackEditDistanceRandom(X, Y, m - 1, n - 1, dp, priorities, generator);
			s1[m - 1].CDR = "R";
			s1[m - 1].opParam = string(1, Y[n - 1]);
			return s1;
		}

		else if (chooseIns) { 		// insert
			vector<LetterOps> s2 = BacktrackEditDistanceRandom(X, Y, m, n - 1, dp, priorities, generator);
			s2[m].insert += Y[n - 1];
			return s2;
		}
		else {			// delete
			assert(chooseDel);
			vector<LetterOps> s3 = BacktrackEditDistanceRandom(X, Y, m - 1, n, dp, priorities, generator);
			s3[m - 1].CDR = "D";
			return s3;
		}
	}
}

// priority: R - replace, I - insert, D - delete
// 0 - RID = [3,2,1]
// 1 - RDI = [3,1,2]
// 2 - IRD = [2,3,1]
// 3 - IDR = [1,3,2]
// 4 - DRI = [2,1,3]
// 5 - DIR = [1,2,3]
// 6 - Random branch each time
vector<LetterOps> ComputeEditDistancePriority(const string& X, const string& Y, const int priority,
		mt19937& generator) {

	int m = X.length(), n = Y.length();
	vector<vector<int> > dp(m + 1, vector<int>(n + 1));
	vector<LetterOps> result;
	EditDistanceArray(X, Y, m, n, dp);
	if (priority == 0) {
		result = BacktrackEditDistancePriority(X, Y, m, n, dp, 3, 2, 1);
	}
	else if (priority == 1) {
		result = BacktrackEditDistancePriority(X, Y, m, n, dp, 3, 1, 2);
	}
	else if (priority == 2) {
		result = BacktrackEditDistancePriority(X, Y, m, n, dp, 2, 3, 1);
	}
	else if (priority == 3) {
		result = BacktrackEditDistancePriority(X, Y, m, n, dp, 1, 3, 2);
	}
	else if (priority == 4) {
		result = BacktrackEditDistancePriority(X, Y, m, n, dp, 2, 1, 3);
	}
	else if (priority == 5) {
		result = BacktrackEditDistancePriority(X, Y, m, n, dp, 1, 2, 3);
	}
	else {
		assert(priority == 6);
		vector<int> priorities { 1, 2, 3 };
		result = BacktrackEditDistanceRandom(X, Y, m, n, dp, priorities, generator);
	}
	return result;
}


//vector<LetterOps> ComputeEditDistanceOperations(const string& X, const string& Y) {
//
//	int m = X.length(), n = Y.length();
//	vector<vector<int> > dp(m + 1, vector<int>(n + 1));
//
//	EditDistanceArray(X, Y, m, n, dp);
//	return BacktrackEditDistanceOperations(X, Y, m, n, dp);
//}

int ComputeEditDistanceNum(const string& X, const string& Y) {

	int m = X.length(), n = Y.length();
	vector<vector<int> > dp(m + 1, vector<int>(n + 1));

	return EditDistanceArray(X, Y, m, n, dp);
}

map<string, double> CountOperations(const vector<LetterOps>& opList) {
	map<string, double> count;
	count["C"] = 0;
	count["I"] = 0;
	count["D"] = 0;
	count["R"] = 0;

	for (vector<LetterOps>::const_iterator it = opList.begin(); it != opList.end(); it++) {
		if (not it->CDR.empty()) {
			count[it->CDR]++;
		}
		count["I"] += it->insert.size();
	}
	return count;
}

vector<LetterOps> RemoveSubstitutions(const vector<LetterOps>& ops) {
	vector<LetterOps> newOps = ops;
	for (unsigned index = 0; index < ops.size(); index++) {
		if (newOps[index].CDR == "R") {
			newOps[index].insert += newOps[index].opParam;
			newOps[index].CDR = "D";
		}
	}
	return newOps;
}

string BuildStringByOperations(const vector<LetterOps>& ops) {
	string newString;
	for (unsigned index = 0; index < ops.size(); index++) {
		newString += ops[index].insert;
		if (ops[index].CDR != "D") {
			newString += ops[index].opParam;
		}
	}
	return newString;
}

vector<LetterOps> ReverseOps(const vector<LetterOps>& ops) {
	int opNum = ops.size();
	vector<LetterOps> reversed(opNum);
	for (int index = opNum - 1; index >= 0; index--) {
		reversed[opNum - 1 - index] = ops[index];
	}
	reversed.push_back(LetterOps());			//add 1 operation for possible insert
	string prevInsString, currInsString;
	prevInsString = reversed[0].insert;
	for (int index = 1; index < opNum + 1; index++) {
		currInsString = reversed[index].insert;
		reverse(prevInsString.begin(), prevInsString.end());
		reversed[index].insert = prevInsString;
		prevInsString = currInsString;
	}
	//delete first LetterOps which is now empty
	reversed.erase(reversed.begin());
	return reversed;
}

vector<LetterOps> ComputeEditDistancePriorityReverse(const string& X, const string& Y, const int priority,
		mt19937& generator){
	string XReverse = X, YReverse=Y;
	reverse(XReverse.begin(),XReverse.end());
	reverse(YReverse.begin(),YReverse.end());
	vector<LetterOps> revOps = ComputeEditDistancePriority(XReverse,YReverse,priority,generator);
	return ReverseOps(revOps);
}

//void TestEdit() {
//	string str1 = "xxxx";
//	string str2 = "yyyy";
//	vector<LetterOps> s = ComputeEditDistancePriority(str1, str2, 0);
//	for (unsigned index = 0; index < s.size(); index++) {
//		cout << s[index].insert << '\t' << s[index].CDR << '\t' << s[index].opParam << endl;
//	}
//	cout << BuildStringByOperations(s) << endl;
//}

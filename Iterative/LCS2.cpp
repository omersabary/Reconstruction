#include "LCS2.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
using namespace std;

const int NA = -1;

LCS2::LCS2(const string& string1, const string& string2) :
		string1(string1), string2(string2), commonSubstrings(), commonSubstringsIndex(-1), commonActive(false) {

}

void LCS2::AddLetter(const char& letter, const int firstIndex, const int secondIndex) {
	if (letter == '_') {
		commonActive = false;
	}
	else { // letter is from GTAC
		if (commonActive) { // There is an active common substring. just add to it.
			commonSubstrings[commonSubstringsIndex].AddLetter(letter);
		}
		else { // no active common substrings. Open a new one.
			commonSubstrings.push_back(CommonSubstring2(letter, firstIndex, secondIndex));
			commonSubstringsIndex++;
			commonActive = true;
		}
	}
}

int LCS2::Len() const {
	int totalLen = 0;
	for (vector<CommonSubstring2>::const_iterator it = commonSubstrings.begin(); it != commonSubstrings.end(); it++) {
		totalLen += it->Len();
	}
	return totalLen;
}

vector<vector<int> > LCS2::Edges() const {
	vector<vector<int> > edges;
	for (vector<CommonSubstring2>::const_iterator it = commonSubstrings.begin(); it != commonSubstrings.end(); it++) {
		vector<int> current(3);
		int start1, end1, start2;
		it->StartEnd(&start1, &end1, &start2);
		current[0] = start1;
		current[1] = end1;
		current[2] = start2;
		edges.push_back(current);
	}
	return edges;
}

vector<int> LCS2::String1LCSIndexes() const {
	vector<int> string1Indexes(string1.size());
	for (vector<CommonSubstring2>::const_iterator it = commonSubstrings.begin(); it != commonSubstrings.end(); it++) {
		pair<int, int> range = it->Range1();
		for (int index = range.first; index < range.second; index++) {
			string1Indexes[index] = 1;
		}
	}
	return string1Indexes;
}

vector<int> LCS2::String1Mirrors() const {
	int start1, end1, start2;
	int vectorSize = string1.size();
	vector<vector<int> > edges = Edges();
	vector<int> mirrors(vectorSize, NA);
	for (vector<vector<int> >::iterator it = edges.begin(); it != edges.end(); it++) {
		start1 = (*it)[0];
		end1 = (*it)[1];
		start2 = (*it)[2];
		for (int index = start1; index < end1; index++) {
			mirrors[index] = start2 + index - start1;
		}
	}
	return mirrors;
}

vector<string> LCS2::DelStringGaps(const vector<int>& string1Mirrors) const {
	int vectorSize = string1.size() + 1;
	vector<string> string1MirrorsGaps(vectorSize);
	int lastIndex = -1;
	for (unsigned int index = 0; index < string1Mirrors.size(); index++) {
		if (string1Mirrors[index] != NA) {
			int startPos = lastIndex + 1;
			int substrLen = string1Mirrors[index] - startPos;
			string1MirrorsGaps[index] = string2.substr(startPos, substrLen);
			lastIndex = string1Mirrors[index];
		}
	}

	// End Gap
	int startPos = lastIndex + 1;
	string1MirrorsGaps[vectorSize - 1] = string2.substr(startPos);

	return string1MirrorsGaps;
}

string LCS2::RepStringWithLowercaseGap(const int patternLen, const int index, const vector<int>& string1Mirrors,
		const vector<string>& string1MirrorsGaps) const {
	assert(patternLen > 1);
	string posStr;
	// handle first letter in string1
	if (index == 0) {
		posStr += 'S';
		for (int j = 0; j < patternLen - 1; j++) {
			if (string1Mirrors[index + j] != NA) {
				string gap2 = string1MirrorsGaps[index + j];
				transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
				posStr += gap2;
				posStr += string1[index + j];
			}
			else {
				posStr += NOT_IN_LCS;
			}
		}
		return posStr;
	}

	int lettersToEnd = string1.size() - index;
	int maxPatternLen = ((patternLen - 1) < lettersToEnd) ? patternLen - 1 : lettersToEnd;
	for (int j = -1; j < maxPatternLen; j++) {
		if (string1Mirrors[index + j] != NA) {
			if (j > -1) {
				string gap2 = string1MirrorsGaps[index + j];
				transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
				posStr += gap2;
			}
			posStr += string1[index + j];
		}
		else {
			posStr += NOT_IN_LCS;
		}
	}

	if (patternLen - 1 > lettersToEnd) {
		string gap2 = string1MirrorsGaps[string1.size()];
		transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
		posStr += gap2;
		posStr += 'S';
	}
	return posStr;
}

string LCS2::RepStringWithNumberGap(const int patternLen, const int index, const vector<int>& string1Mirrors,
		const vector<string>& string1MirrorsGaps) const {
	assert(patternLen > 1);
	string posStr;
	// handle first letter in string1
	if (index == 0) {
		posStr += 'S';
		for (int j = 0; j < patternLen - 1; j++) {
			if (string1Mirrors[index + j] != NA) {
				string gap2str = string1MirrorsGaps[index + j];
				string gap2 = to_string(gap2str.size());
				posStr += gap2;
				posStr += string1[index + j];
			}
			else {
				posStr += NOT_IN_LCS;
			}
		}
		return posStr;
	}

	int lettersToEnd = string1.size() - index;
	int maxPatternLen = ((patternLen - 1) < lettersToEnd) ? patternLen - 1 : lettersToEnd;
	for (int j = -1; j < maxPatternLen; j++) {
		if (string1Mirrors[index + j] != NA) {
			if (j > -1) {
				string gap2str = string1MirrorsGaps[index + j];
				string gap2 = to_string(gap2str.size());
				posStr += gap2;
			}
			posStr += string1[index + j];
		}
		else {
			posStr += NOT_IN_LCS;
		}
	}

	if (patternLen - 1 > lettersToEnd) {
		string gap2str = string1MirrorsGaps[string1.size()];
		string gap2 = to_string(gap2str.size());
		posStr += gap2;
		posStr += 'S';
	}
	return posStr;
}

string LCS2::RepStringNoGap(const int patternLen, const int index, const vector<int>& string1Mirrors) const {
	assert(patternLen > 1);
	string posStr;
	// handle first letter in string1
	if (index == 0) {
		posStr += 'S';
		for (int j = 0; j < patternLen - 1; j++) {
			if (string1Mirrors[index + j] != NA) {
				posStr += string1[index + j];
			}
			else {
				posStr += NOT_IN_LCS;
			}
		}
		return posStr;
	}

	int lettersToEnd = string1.size() - index;
	int maxPatternLen = ((patternLen - 1) < lettersToEnd) ? patternLen - 1 : lettersToEnd;
	for (int j = -1; j < maxPatternLen; j++) {
		if (string1Mirrors[index + j] != NA) {
			posStr += string1[index + j];
		}
		else {
			posStr += NOT_IN_LCS;
		}
	}

	if (patternLen - 1 > lettersToEnd) {
		posStr += 'S';
	}
	return posStr;
}

vector<string> LCS2::RepStringWithLowercaseGapArray(const int patternLen) const {
	int vectorLen = string1.size() + 1;
	vector<int> string1Mirrors = String1Mirrors();
	vector<string> string1MirrorsGaps = DelStringGaps(string1Mirrors);
	vector<string> repArray;
	for (int index = 0; index < vectorLen; index++) {
		string posStr = RepStringWithLowercaseGap(patternLen, index, string1Mirrors, string1MirrorsGaps);
		repArray.push_back(posStr);
	}
	return repArray;
}

vector<string> LCS2::RepStringWithNumberGapArray(const int patternLen) const {
	int vectorLen = string1.size() + 1;
	vector<int> string1Mirrors = String1Mirrors();
	vector<string> string1MirrorsGaps = DelStringGaps(string1Mirrors);
	vector<string> repArray;
	for (int index = 0; index < vectorLen; index++) {
		string posStr = RepStringWithNumberGap(patternLen, index, string1Mirrors, string1MirrorsGaps);
		repArray.push_back(posStr);
	}
	return repArray;
}

vector<string> LCS2::RepStringNoGapArray(const int patternLen) const {
	int vectorLen = string1.size() + 1;
	vector<int> string1Mirrors = String1Mirrors();
	vector<string> repArray;
	for (int index = 0; index < vectorLen; index++) {
		string posStr = RepStringNoGap(patternLen, index, string1Mirrors);
		repArray.push_back(posStr);
	}
	return repArray;
}

std::ostream& operator<<(std::ostream& os, const LCS2& a) {
	if (a.commonSubstrings.empty()) {
		return os;
	}
	for (unsigned commonIndex = 0; commonIndex < a.commonSubstrings.size() - 1; commonIndex++) {
		os << a.commonSubstrings[commonIndex];
		os << "_";
	}
	os << a.commonSubstrings[a.commonSubstrings.size() - 1];
	return os;
}

/*vector<LCS2> FindLCS(const string& X, const string& Y, int m, int n,
 int** L) {

 // If we reach end of either string, return a vector with empty LCS2
 if (m == 0 || n == 0) {
 return vector<LCS2>(1, LCS2(X, Y));
 }

 // If the last characters of X and Y are the same
 if (X[m - 1] == Y[n - 1]) {
 // recurse for X[0..m-2] and Y[0..n-2] in the matrix
 vector<LCS2> s = FindLCS(X, Y, m - 1, n - 1, L);

 // append current character to all possible LCS of substring X[0..m-2] and Y[0..n-2].
 for (vector<LCS2>::iterator lcs2 = s.begin(); lcs2 != s.end(); lcs2++)
 lcs2->AddLetter(X[m - 1], m - 1, n - 1);
 return s;
 }

 // If the last characters of X and Y are not same
 else {
 vector<LCS2> s;
 // If LCS can be constructed from top side of the matrix, recurse for X[0..m-2] and Y[0..n-1]
 if (L[m - 1][n] >= L[m][n - 1])
 s = FindLCS(X, Y, m - 1, n, L);

 // If LCS can be constructed from left side of the matrix, recurse for X[0..m-1] and Y[0..n-2]
 if (L[m][n - 1] >= L[m - 1][n]) {
 vector<LCS2> tmp = FindLCS(X, Y, m, n - 1, L);

 // merge two vectors if L[m-1][n] == L[m][n-1] Note s will be empty if L[m-1][n] != L[m][n-1]
 s.insert(s.end(), tmp.begin(), tmp.end());
 }
 for (vector<LCS2>::iterator lcs2 = s.begin(); lcs2 != s.end(); lcs2++)
 lcs2->AddLetter('_', m - 1, n - 1);
 return s;
 }
 }*/

// find one LCS
LCS2 FindLCS(const string& X, const string& Y, int m, int n, const vector<vector<int> >& L) {

// If we reach end of either string, return a vector with empty LCS2
	if (m == 0 || n == 0) {
		return LCS2(X, Y);
	}

// If the last characters of X and Y are the same
	if (X[m - 1] == Y[n - 1]) {
		// recurse for X[0..m-2] and Y[0..n-2] in the matrix
		LCS2 lcs2 = FindLCS(X, Y, m - 1, n - 1, L);

		// append current character to LCS of substring X[0..m-2] and Y[0..n-2].
		lcs2.AddLetter(X[m - 1], m - 1, n - 1);
		return lcs2;
	}

// If the last characters of X and Y are not same
	else {

		// If LCS can be constructed from top side of the matrix, recurse for X[0..m-2] and Y[0..n-1]
		if (L[m - 1][n] >= L[m][n - 1]) {
			LCS2 lcs2 = FindLCS(X, Y, m - 1, n, L);
			lcs2.AddLetter('_', m - 1, n);
			return lcs2;
		}
		// If LCS can be constructed from left side of the matrix, recurse for X[0..m-1] and Y[0..n-2]
		else {
			LCS2 lcs2 = FindLCS(X, Y, m, n - 1, L);
			lcs2.AddLetter('_', m, n - 1);
			return lcs2;
		}
	}
}

/* Returns length of LCS for X[0..m-1], Y[0..n-1] */
int LCSLength(const string& X, const string& Y, int m, int n, vector<vector<int> >& L) {
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

LCS2 ComputeLCS2(const string& X, const string& Y) {

	int m = X.length(), n = Y.length();
	vector<vector<int> > L(m + 1, vector<int>(n + 1));

	LCSLength(X, Y, m, n, L);
	return FindLCS(X, Y, m, n, L);
}

int ComputeEditDistWithoutReplace(const string& X, const string& Y) {

	int m = X.length(), n = Y.length();
	vector<vector<int> > L(m + 1, vector<int>(n + 1));

	int lcsLen = LCSLength(X, Y, m, n, L);
	return X.size() - lcsLen + Y.size() - lcsLen;
}


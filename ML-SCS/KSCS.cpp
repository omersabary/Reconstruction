#include "KSCS.hpp"
#include <set>

using namespace std;

// Common supersequence with "common subsequence" of length csubLen. if possible, return csubLen.
// Otherwise, return INT_MAX

int CSupByCSubLen(const string& X, const string&Y, int m, int n, int csubLen) {

	if (m == 0) {
		if (csubLen == 0) {
			return n;
		}
		else {
			return INT_MAX;
		}
	}
	if (n == 0) {
		if (csubLen == 0) {
			return m;
		}
		else {
			return INT_MAX;
		}
	}

	int targetSCSLen = m + n - csubLen;
	if (X[m - 1] == Y[n - 1]) {
		if (csubLen > 0) {
			if (CSupByCSubLen(X, Y, m - 1, n - 1, csubLen - 1) + 1 == targetSCSLen) {
				return targetSCSLen;
			}
		}
		if (CSupByCSubLen(X, Y, m - 1, n, csubLen) + 1 == targetSCSLen) {
			return targetSCSLen;
		}
		if (CSupByCSubLen(X, Y, m, n - 1, csubLen) + 1 == targetSCSLen) {
			return targetSCSLen;
		}
		return INT_MAX;
	}
	else {
		if (CSupByCSubLen(X, Y, m - 1, n, csubLen) + 1 == targetSCSLen) {
			return targetSCSLen;
		}
		if (CSupByCSubLen(X, Y, m, n - 1, csubLen) + 1 == targetSCSLen) {
			return targetSCSLen;
		}
		return INT_MAX;
	}

}

// One Common supersequence with "common subsequence" of length csubLen. if possible, return common supersequence.
// Otherwise, return "!"
string CSupByCSubLenString(const string& X, const string&Y, int m, int n, int csubLen) {
	if (m == 0) {
		if (csubLen == 0) {
			return Y.substr(0, n);
		}
		else {
			return "!";
		}
	}
	if (n == 0) {
		if (csubLen == 0) {
			return X.substr(0, m);
		}
		else {
			return "!";
		}
	}

	int targetSCSLen = m + n - csubLen;
	string current;
	int currentSize;
	if (X[m - 1] == Y[n - 1]) {
		if (csubLen > 0) {
			current = CSupByCSubLenString(X, Y, m - 1, n - 1, csubLen - 1);
			currentSize = current.size();
			if (current != "!" and currentSize + 1 == targetSCSLen) {
				return current + X[m - 1];
			}
		}
		current = CSupByCSubLenString(X, Y, m - 1, n, csubLen);
		currentSize = current.size();
		if (current != "!" and currentSize + 1 == targetSCSLen) {
			return current + X[m - 1];
		}
		current = CSupByCSubLenString(X, Y, m, n - 1, csubLen);
		currentSize = current.size();
		if (current != "!" and currentSize + 1 == targetSCSLen) {
			return current + X[m - 1];
		}
		return "!";
	}
	else {
		current = CSupByCSubLenString(X, Y, m - 1, n, csubLen);
		currentSize = current.size();
		if (current != "!" and currentSize + 1 == targetSCSLen) {
			return current + X[m - 1];
		}
		current = CSupByCSubLenString(X, Y, m, n - 1, csubLen);
		currentSize = current.size();
		if (current != "!" and currentSize + 1 == targetSCSLen) {
			return current + Y[n - 1];
		}
		return "!";
	}

}
// All Common supersequences with "common subsequence" of length csubLen. if possible, return common subsequences.
// Otherwise, return empty vector.
vector<string> CSupByCSubLenAll(const string& X, const string&Y, int m, int n, int csubLen) {
	vector<string> result;
	if (m == 0) {
		if (csubLen == 0) {
			result.push_back(Y.substr(0, n));
		}
		return result;
	}
	if (n == 0) {
		if (csubLen == 0) {
			result.push_back(X.substr(0, m));
		}
		return result;
	}

	int targetSCSLen = m + n - csubLen;
	int currentLen;
	if (X[m - 1] == Y[n - 1]) {
		vector<string> v1, v2, v3;
		if (csubLen > 0) {
			v1 = CSupByCSubLenAll(X, Y, m - 1, n - 1, csubLen - 1);
			if (not v1.empty()) {
				currentLen = v1[0].size();
				if (currentLen + 1 != targetSCSLen) {
					v1.clear();
				}
			}
		}
		v2 = CSupByCSubLenAll(X, Y, m - 1, n, csubLen);
		if (not v2.empty()) {
			currentLen = v2[0].size();
			if (currentLen + 1 != targetSCSLen) {
				v2.clear();
			}
		}
		v3 = CSupByCSubLenAll(X, Y, m, n - 1, csubLen);
		if (not v3.empty()) {
			currentLen = v3[0].size();
			if (currentLen + 1 != targetSCSLen) {
				v3.clear();
			}
		}

		result.insert(result.end(), v1.begin(), v1.end());
		if (v2 != v1) {
			result.insert(result.end(), v2.begin(), v2.end());
		}
		if (v3 != v1 and v3 != v2) {
			result.insert(result.end(), v3.begin(), v3.end());
		}
		// add last letter
		for (vector<string>::iterator it = result.begin(); it != result.end(); it++) {
			it->push_back(X[m - 1]);
		}
		return result;
	}
	else {
		vector<string> v2, v3;
		v2 = CSupByCSubLenAll(X, Y, m - 1, n, csubLen);
		if (not v2.empty()) {
			currentLen = v2[0].size();
			if (currentLen + 1 != targetSCSLen) {
				v2.clear();
			}
			else {
				for (vector<string>::iterator it = v2.begin(); it != v2.end(); it++) {
					it->push_back(X[m - 1]);
				}
			}
		}
		v3 = CSupByCSubLenAll(X, Y, m, n - 1, csubLen);
		if (not v3.empty()) {
			currentLen = v3[0].size();
			if (currentLen + 1 != targetSCSLen) {
				v3.clear();
			}
			else {
				for (vector<string>::iterator it = v3.begin(); it != v3.end(); it++) {
					it->push_back(Y[n - 1]);
				}
			}
		}
		result.insert(result.end(), v2.begin(), v2.end());
		result.insert(result.end(), v3.begin(), v3.end());
		return result;
	}

}

vector<string> KCS(const string& X, const string&Y, int csubLen) {
	int m = X.size();
	int n = Y.size();
	vector<string> CS = CSupByCSubLenAll(X, Y, m, n, csubLen);
	set<string> setCS(CS.begin(), CS.end());
	CS.assign(setCS.begin(), setCS.end());
	return CS;
}

int KCSDPLen(const string& X, const string&Y, const int csubLen, int*** lookup) {
	int m = X.size(), n = Y.size();
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int csub = 0; csub <= csubLen; csub++) {
				lookup[i][j][csub] = INT_MAX;
				int targetCSLen = i + j - csub;
				if (i == 0) {
					if (csub == 0) {
						lookup[i][j][csub] = j;
					}
				}
				else if (j == 0) {
					if (csub == 0) {
						lookup[i][j][csub] = i;
					}
				}
				else if (X[i - 1] == Y[j - 1]) {
					if (csub > 0 and lookup[i - 1][j - 1][csub - 1] + 1 == targetCSLen) {
						lookup[i][j][csub] = targetCSLen;
					}
					else if (lookup[i - 1][j][csub] + 1 == targetCSLen) {
						lookup[i][j][csub] = targetCSLen;
					}
					else if (lookup[i][j - 1][csub] + 1 == targetCSLen) {
						lookup[i][j][csub] = targetCSLen;
					}
				}
				else {
					if (lookup[i - 1][j][csub] + 1 == targetCSLen) {
						lookup[i][j][csub] = targetCSLen;
					}
					else if (lookup[i][j - 1][csub] + 1 == targetCSLen) {
						lookup[i][j][csub] = targetCSLen;
					}
				}
			}
		}
	}
	return lookup[m][n][csubLen];
}

set<string> KCSBacktrack(const string& X, const string&Y, const int m, const int n, const int csubLen, int*** lookup) {
	set<string> result;
	if (m == 0) {
		if (csubLen == 0) {
			result.insert(Y.substr(0, n));
		}
		return result;
	}
	if (n == 0) {
		if (csubLen == 0) {
			result.insert(X.substr(0, m));
		}
		return result;
	}

	int targetSCSLen = m + n - csubLen;
	if (X[m - 1] == Y[n - 1]) {
		set<string> v1, v2, v3, all;
		if (csubLen > 0) {
			if (lookup[m - 1][n - 1][csubLen - 1] + 1 == targetSCSLen) {
				v1 = KCSBacktrack(X, Y, m - 1, n - 1, csubLen - 1, lookup);
			}
		}
		if (lookup[m - 1][n][csubLen] + 1 == targetSCSLen) {
			v2 = KCSBacktrack(X, Y, m - 1, n, csubLen, lookup);
		}
		if (lookup[m][n - 1][csubLen] + 1 == targetSCSLen) {
			v3 = KCSBacktrack(X, Y, m, n - 1, csubLen, lookup);
		}

		all = v1;
		all.insert(v2.begin(), v2.end());
		all.insert(v3.begin(), v3.end());
		vector<string> iterable(all.begin(), all.end());

		// add last letter
		for (vector<string>::iterator it = iterable.begin(); it != iterable.end(); it++) {
			it->push_back(X[m - 1]);
		}
		result.insert(iterable.begin(), iterable.end());
		return result;
	}
	else {
		set<string> v2, v3;
		if (lookup[m - 1][n][csubLen] + 1 == targetSCSLen) {
			v2 = KCSBacktrack(X, Y, m - 1, n, csubLen, lookup);
			vector<string> iterable(v2.begin(), v2.end());
			for (vector<string>::iterator it = iterable.begin(); it != iterable.end(); it++) {
				it->push_back(X[m - 1]);
			}
			result.insert(iterable.begin(), iterable.end());
		}

		if (lookup[m][n - 1][csubLen] + 1 == targetSCSLen) {
			v3 = KCSBacktrack(X, Y, m, n - 1, csubLen, lookup);
			vector<string> iterable(v3.begin(), v3.end());
			for (vector<string>::iterator it = iterable.begin(); it != iterable.end(); it++) {
				it->push_back(Y[n - 1]);
			}
			result.insert(iterable.begin(), iterable.end());
		}

		return result;
	}

}

vector<string> KCSDP(const string& X, const string&Y, const int scslen) {
	vector<string> result;
	int m = X.size(), n = Y.size();
	int csubLen = m+n - scslen;
	int ***lookup;
	lookup = new int**[m + 1];
	for (int i = 0; i < m + 1; i++) {
		lookup[i] = new int*[n + 1];
		for (int j = 0; j < n + 1; j++) {
			lookup[i][j] = new int[csubLen + 1];
		}
	}
	// Fill lookup table.
	int kcsLen = KCSDPLen(X, Y, csubLen, lookup);

	if (kcsLen != INT_MAX) {
		set<string> setRes = KCSBacktrack(X, Y, m, n, csubLen, lookup);
		result.assign(setRes.begin(), setRes.end());
	}

	// Release memory
	for (int i = 0; i < m + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			delete[] lookup[i][j];
		}
		delete[] lookup[i];
	}
	delete[] lookup;
	return result;
}

#include <iostream>
#include <set>
#include "MultiString.hpp"

using namespace std;

const bool LCS_STR = true;
const bool NON_LCS_STR = false;

void MultiString::AddLetterLCS(const char& letter) {
	int vecNum = sequence.size();
	if (vecNum == 0 or isLCS[vecNum - 1] == NON_LCS_STR) { // new vector for LCS letters
		vector<string> newVec = { string(1, letter) };
		sequence.push_back(newVec);
		isLCS.push_back(LCS_STR);
	}
	else { // add to current LCS string (only string in vector)
		sequence[vecNum - 1][0] += letter;
	}

}
void MultiString::AddLetterNONLCS(const char& letter) {
	int vecNum = sequence.size();
	if (vecNum == 0 or isLCS[vecNum - 1] == LCS_STR) { // new vector for NON LCS letters
		vector<string> newVec = { string(1, letter) };
		sequence.push_back(newVec);
		isLCS.push_back(NON_LCS_STR);
	}
	else { // add letter to all strings in current vector
		int strNum = sequence[vecNum - 1].size();
		for (int i = 0; i < strNum; i++) {
			sequence[vecNum - 1][i] += letter;
		}
	}
}

void MultiString::AddStringNONLCS(const std::string& str) {
	int vecNum = sequence.size();
	if (vecNum == 0 or isLCS[vecNum - 1] == LCS_STR) { // new vector for NON LCS letters
		vector<string> newVec = { str };
		sequence.push_back(newVec);
		isLCS.push_back(NON_LCS_STR);
	}
	else { // add letter to all strings in current vector
		int strNum = sequence[vecNum - 1].size();
		for (int i = 0; i < strNum; i++) {
			sequence[vecNum - 1][i] += str;
		}
	}
}

bool MultiString::empty() const {
	return sequence.empty();
}
void MultiString::Print() const {
	if (empty()) {
		cout << "<empty>" << endl;
		return;
	}
	for (unsigned i = 0; i < sequence.size(); i++) {
		if (isLCS[i] == LCS_STR) {
			cout << sequence[i][0];
		}
		else {
			cout << "{";
			for (unsigned j = 0; j < sequence[i].size() - 1; j++) {
				cout << sequence[i][j] << ", ";
			}
			cout << sequence[i][sequence[i].size() - 1] << "}";
		}
		if (i != sequence.size() - 1) {
			cout << " ";
		}
	}
	cout << endl;
}

// If multistrings are the same return one of them.
// If last string vector is NON LCS for both and the rest is the same, merge unique strings of the last string vector
// otherwise, return empty multistring
MultiString Merge(const MultiString& a, const MultiString& b) {
	int vecNum = a.sequence.size();
	if (vecNum == 0) {
		return MultiString();
	}
	if (a.isLCS != b.isLCS) {
		return MultiString();
	}
	else { // isLCS vector is the same -> same sequence length.
		   // test if all vectors but last are the same. if not, return empty multistring
		for (int i = 0; i < vecNum - 1; i++) {
			if (a.sequence[i] != b.sequence[i]) {
				return MultiString();
			}
		}
		// at this point: isLCS vectors are the same. sequence is the same except(maybe) its last vector.
		if (a.sequence[vecNum - 1] == b.sequence[vecNum - 1]) { // multistrings are identical, return one.
			return a;
		}
		else if (a.isLCS[vecNum - 1] == LCS_STR) { // last string vector is different and it is LCS, return empty multistring.
			return MultiString();
		}
		else { // last string vector is different and it is NON LCS, merge unique strings of last vector
			MultiString result = a;
			const vector<string>& aLastVec = a.sequence[vecNum - 1];
			const vector<string>& bLastVec = b.sequence[vecNum - 1];
			set<string> unique(aLastVec.begin(), aLastVec.end());
			unique.insert(bLastVec.begin(), bLastVec.end());
			//vector<string> merged(unique.begin(), unique.end());
			result.sequence[vecNum - 1] = vector<string>(unique.begin(), unique.end());
			return result;
		}
	}
}

vector<MultiString> MergeVectors(const vector<MultiString>& a, const vector<MultiString>& b) {
	vector<MultiString> result(a);
	vector<int> bIndexesToAdd;
	for (unsigned i = 0; i < b.size(); i++) {
		bool bAdded = false;
		for (unsigned j = 0; j < a.size(); j++) {
			MultiString temp = Merge(a[j], b[i]);
			if (not temp.empty()) { //a[i] and b[j] were merged
				result[j] = temp;
				bAdded = true;
				break;
			}
		}
		if (not bAdded) {
			bIndexesToAdd.push_back(i);
		}
	}
	for (unsigned k = 0; k < bIndexesToAdd.size(); k++) {
		result.push_back(b[bIndexesToAdd[k]]);
	}
	return result;
}

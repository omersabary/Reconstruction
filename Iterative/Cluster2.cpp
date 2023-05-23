#include <set>
#include <chrono>
#include <random>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <cctype>
#include <unordered_set>
#include "Cluster2.hpp"
#include "LongestPath.hpp"
using namespace std;

int CorrectSizeCloneCount(const vector<Clone>& clones, const int correctSize);
int CumCloneToClonesED(const int cloneIndex, const vector<Clone>& clones);

Cluster2::Cluster2(const unsigned strandLength, const int clonesNum, const double delProb, const double insProb,
		const double subProb, mt19937& generator) {

	// make original
	original = MakeStrand(strandLength, generator);

	unordered_set<Clone, CloneHasher> tmpSet; // remember  ordered set can help choose a clone with less inserts
	// make clones
	while ((int) tmpSet.size() < clonesNum) {
		tmpSet.insert(Clone(original, delProb, insProb, subProb, generator));
	}
	clones.insert(clones.end(), tmpSet.begin(), tmpSet.end());
	// shuffle to make order random in clones.
	shuffle(clones.begin(), clones.end(), generator);
	clonesBackup = clones;

}

Cluster2::Cluster2(const string& orig, const vector<string>& copies) {
	original = orig;
	for (unsigned i = 0; i < copies.size(); i++) {
		clones.push_back(Clone(orig, copies[i]));
	}
	clonesBackup = clones;
}

void Cluster2::ReverseClones() {
	for (unsigned i = 0; i < clones.size(); i++) {
		clones[i].Reverse();
	}
}

std::string Cluster2::Original() const {
	return original;
}

struct compare {
	bool operator()(const pair<int, int>& first, const pair<int, int>& second) {
		return first.first < second.first;
	}
};

struct compareBreakTies {
	bool operator()(const pair<pair<int, int>, int>& first, const pair<pair<int, int>, int>& second) {
		if (first.first.first == second.first.first) {
			return first.first.second < second.first.second;
		}
		else {
			return first.first.first < second.first.first;
		}
	}
};

struct compareStr {
	bool operator()(const pair<string, int>& first, const pair<string, int>& second) {
		return first.first < second.first;
	}
};

vector<Clone> SortClonesByED(const vector<Clone>& clones) {
	compare c;
	int clonesNum = clones.size();
	vector<pair<int, int> > cumED(clonesNum);
	for (int i = 0; i < clonesNum; i++) {
		cumED[i].first = CumCloneToClonesED(i, clones);
		cumED[i].second = i;
	}
	sort(cumED.begin(), cumED.end(), c);
	vector<Clone> sortedClones;
	for (int i = 0; i < clonesNum; i++) {
		int oldIndex = cumED[i].second;
		sortedClones.push_back(clones[oldIndex]);
	}
	return sortedClones;
}

vector<Clone> SortClonesByLen(const vector<Clone>& clones) {
	compare c;
	int clonesNum = clones.size();
	vector<pair<int, int> > clonesLen(clonesNum);
	for (int i = 0; i < clonesNum; i++) {
		clonesLen[i].first = clones[i].Len();
		clonesLen[i].second = i;
	}
	sort(clonesLen.begin(), clonesLen.end(), c);
	vector<Clone> sortedClones;
	for (int i = 0; i < clonesNum; i++) {
		int oldIndex = clonesLen[i].second;
		sortedClones.push_back(clones[oldIndex]);
	}
	return sortedClones;
}

vector<Clone> SortClonesByAbsLenDiff(const vector<Clone>& clones, const int correctSize) {
	compare c;
	int clonesNum = clones.size();
	vector<pair<int, int> > clonesLen(clonesNum);
	for (int i = 0; i < clonesNum; i++) {
		clonesLen[i].first = abs(clones[i].Len() - correctSize);
		clonesLen[i].second = i;
	}
	sort(clonesLen.begin(), clonesLen.end(), c);
	vector<Clone> sortedClones;
	for (int i = 0; i < clonesNum; i++) {
		int oldIndex = clonesLen[i].second;
		sortedClones.push_back(clones[oldIndex]);
	}
	return sortedClones;
}

vector<Clone> SortClonesLex(const vector<Clone>& clones) {
	compareStr c;
	int clonesNum = clones.size();
	vector<pair<string, int> > clonesLex(clonesNum);
	for (int i = 0; i < clonesNum; i++) {
		clonesLex[i].first = clones[i].String();
		clonesLex[i].second = i;
	}
	sort(clonesLex.begin(), clonesLex.end(), c);
	vector<Clone> sortedClones;
	for (int i = 0; i < clonesNum; i++) {
		int oldIndex = clonesLex[i].second;
		sortedClones.push_back(clones[oldIndex]);
	}
	return sortedClones;
}

vector<Clone> SortClonesLexInsideOut(const vector<Clone>& clones) {
	compareStr c;
	int clonesNum = clones.size();
	vector<pair<string, int> > clonesLex(clonesNum);
	for (int i = 0; i < clonesNum; i++) {
		clonesLex[i].first = clones[i].String();
		clonesLex[i].second = i;
	}
	sort(clonesLex.begin(), clonesLex.end(), c);
	vector<Clone> sortedClones;

	int middleIndex;
	if (clonesNum % 2 == 1) {
		middleIndex = (clonesNum + 1) / 2;
		sortedClones.push_back(clones[middleIndex - 1]);

	}
	else {
		middleIndex = clonesNum / 2;
	}

	for (int i = middleIndex; i < clonesNum; i++) {
		int oldIndex = clonesLex[i].second;
		sortedClones.push_back(clones[oldIndex]);
		oldIndex = clonesLex[clonesNum - i - 1].second;
		sortedClones.push_back(clones[oldIndex]);
	}
	return sortedClones;
}

vector<Clone> SortClonesLexInsideOutReverse(const vector<Clone>& clones) {
	compareStr c;
	int clonesNum = clones.size();
	vector<pair<string, int> > clonesLex(clonesNum);
	for (int i = 0; i < clonesNum; i++) {
		clonesLex[i].first = clones[i].String();
		clonesLex[i].second = i;
	}
	sort(clonesLex.begin(), clonesLex.end(), c);
	vector<Clone> sortedClones;

	int middleIndex;
	if (clonesNum % 2 == 1) {
		middleIndex = (clonesNum + 1) / 2;

	}
	else {
		middleIndex = clonesNum / 2;
	}

	for (int i = clonesNum - 1; i >= middleIndex; i--) {
		int oldIndex = clonesLex[i].second;
		sortedClones.push_back(clones[oldIndex]);
		oldIndex = clonesLex[clonesNum - i - 1].second;
		sortedClones.push_back(clones[oldIndex]);
	}
	if (clonesNum % 2 == 1) {
		sortedClones.push_back(clones[middleIndex - 1]);

	}
	return sortedClones;
}

int Cluster2::GetInsertNum(const int cloneIndex) const {
	return clones[cloneIndex].GetInsertNum();
}

int Cluster2::GetDeleteNum(const int cloneIndex) const {
	return clones[cloneIndex].GetDeleteNum();
}

void Cluster2::CmpLCS2(const unsigned int anchorIndex) {
	// compute lcs2Sets for clone 0 with all the others
	lcs2s.clear();
	for (unsigned int index = 0; index < clones.size(); index++) {
		if (anchorIndex != index) {
			lcs2s.push_back(ComputeLCS2(clones[anchorIndex].String(), clones[index].String()));
		}
	}
}

void Cluster2::ComputeEditLists(const unsigned int anchorIndex, const int priority, mt19937& generator) {
	editLists.clear();
	for (unsigned int index = 0; index < clones.size(); index++) {
		if (anchorIndex != index) {
			editLists.push_back(
					ComputeEditDistancePriority(clones[anchorIndex].String(), clones[index].String(), priority,
							generator));
//			editLists.push_back(
//					ComputeEditDistancePriorityReverse(clones[anchorIndex].String(), clones[index].String(), priority,
//							generator));
		}
	}
}

vector<int> Cluster2::InsArray(const int cloneIndex) {
	vector<int> cum(clones[cloneIndex].Len());
	CmpLCS2(cloneIndex);
	for (vector<LCS2>::iterator it = lcs2s.begin(); it != lcs2s.end(); it++) {
		vector<int> lcIndexes = it->String1LCSIndexes();
		for (unsigned int index = 0; index < lcIndexes.size(); index++) {
			cum[index] += lcIndexes[index];
		}
	}
	return cum;
}

vector<map<string, int> > Cluster2::CumEditOperationsNoInserts(const int cloneIndex, const int priority,
		mt19937& generator) {
	unsigned vectorSize = clones[cloneIndex].Len() + 1;
	vector<map<string, int> > cum(vectorSize);
	ComputeEditLists(cloneIndex, priority, generator);
	for (unsigned listIndex = 0; listIndex < editLists.size(); listIndex++) {
		for (unsigned letterIndex = 0; letterIndex < vectorSize; letterIndex++) {
			string opKey = editLists[listIndex][letterIndex].CDR;
			if (opKey == "R" or opKey == "C") {
				opKey += editLists[listIndex][letterIndex].opParam;
			}
			cum[letterIndex][opKey]++;
		}
	}
	return cum;
}

vector<map<string, int> > Cluster2::CumEditOperationsInserts(const int cloneIndex, const int priority,
		mt19937& generator) {
	unsigned vectorSize = clones[cloneIndex].Len() + 1;
	vector<map<string, int> > cum(vectorSize);
	ComputeEditLists(cloneIndex, priority, generator);
	for (unsigned listIndex = 0; listIndex < editLists.size(); listIndex++) {
		for (unsigned letterIndex = 0; letterIndex < vectorSize; letterIndex++) {
			string insString = editLists[listIndex][letterIndex].insert;
//			if (editLists[listIndex][letterIndex].CDR == "R") {
//				insString += editLists[listIndex][letterIndex].opParam;
//			}
			cum[letterIndex][insString]++;
		}
	}
	return cum;
}

vector<map<string, int> > Cluster2::RepStringLowercaseGapsDictionaries(const int cloneIndex, const int patternLen) {
	unsigned vectorSize = clones[cloneIndex].Len() + 1;
	vector<map<string, int> > dictionaries(vectorSize);
	CmpLCS2(cloneIndex);
	for (vector<LCS2>::iterator it = lcs2s.begin(); it != lcs2s.end(); it++) {
		vector<string> currentVector = it->RepStringWithLowercaseGapArray(patternLen);
		for (unsigned int index = 0; index < vectorSize; index++) {
			dictionaries[index][currentVector[index]]++;
		}
	}
	return dictionaries;
}

vector<string> GapsRDI(const vector<LetterOps>& operations) {
	int vectorSize = operations.size();
	vector<string> gaps(vectorSize);
	string cumGap;
	for (int index = 0; index < vectorSize - 1; index++) {
		if (operations[index].CDR == "D") {
			cumGap += operations[index].insert;
		}
		else if (operations[index].CDR == "R") {
			cumGap += operations[index].opParam;
		}
		else {
			gaps[index] = cumGap + operations[index].insert;
			cumGap.clear();
		}
	}
	// End Gap
	gaps[vectorSize - 1] = cumGap + operations[vectorSize - 1].insert;
	return gaps;
}

string RepStringRDI(const vector<LetterOps>& operations, const vector<string>& gaps, const int patternLen,
		const int index) {
	// DIFFERENT FROM PREVIOUS REP STRINGS - NOT_IN_LCS CAN HAVE DELETED STRING JUST BEFORE IT
	assert(patternLen > 1);
	string posStr;
	// handle first letter in string1
	if (index == 0) {
		posStr += 'S';
		for (int j = 0; j < patternLen - 1; j++) {
			string gap2 = operations[index + j].insert;
//			if (operations[index + j].CDR == "R") {
//				gap2 += operations[index + j].opParam;
//			}
			//string gap2 = gaps[index + j];
			transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
			posStr += gap2;
			if (operations[index + j].CDR == "D") {
				posStr += NOT_IN_LCS;
			}
			else { // CDR == "C" or "R"
				posStr += operations[index + j].opParam;
			}
		}
		return posStr;
	}

	//int lettersToEnd = string1.size() - index;
	int string1Size = operations.size() - 1;
	int lettersToEnd = string1Size - index;
	int maxPatternLen = ((patternLen - 1) < lettersToEnd) ? patternLen - 1 : lettersToEnd;
	for (int j = -1; j < maxPatternLen; j++) {
		if (j > -1) {
			string gap2 = operations[index + j].insert;
//			if (operations[index + j].CDR == "R") {
//				gap2 += operations[index + j].opParam;
//			}
			//string gap2 = gaps[index + j];
			transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
			posStr += gap2;
		}

		if (operations[index + j].CDR == "D") {
			posStr += NOT_IN_LCS;
		}
		else { // CDR == "C" or "R"
			posStr += operations[index + j].opParam;
		}
	}

	if (patternLen - 1 > lettersToEnd) {
		string gap2 = operations[string1Size].insert;
		//string gap2 = gaps[string1Size];
		transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
		posStr += gap2;
		posStr += 'S';
	}
	return posStr;
}

string RepStringRDISubstitutions(const vector<LetterOps>& operations, const int patternLen, const int index) {
	assert(patternLen > 1);
	string posStr;
	// handle first letter in string1
	if (index == 0) {
		posStr += 'S';
		for (int j = 0; j < patternLen - 1; j++) {
			string gap2 = operations[index + j].insert;
			transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
			posStr += gap2;
			if (operations[index + j].CDR == "D") {
				posStr += NOT_IN_LCS;
			}
			else { // CDR == "C" or "R"
				posStr += operations[index + j].opParam;
			}
		}
		return posStr;
	}

	int string1Size = operations.size() - 1;
	int lettersToEnd = string1Size - index;
	int maxPatternLen = ((patternLen - 1) < lettersToEnd) ? patternLen - 1 : lettersToEnd;
	for (int j = -1; j < maxPatternLen; j++) {
		if (j > -1) {
			string gap2 = operations[index + j].insert;
			transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
			posStr += gap2;
		}

		if (operations[index + j].CDR == "D") {
			posStr += NOT_IN_LCS;
		}
		else { // CDR == "C" or "R"
			posStr += operations[index + j].opParam;
		}
	}

	if (patternLen - 1 > lettersToEnd) {
		string gap2 = operations[string1Size].insert;
		transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
		posStr += gap2;
		posStr += 'S';
	}
	return posStr;
}

string RepStringRDISubstitutionsNoGaps(const vector<LetterOps>& operations, const int patternLen, const int index) {
	assert(patternLen > 1);
	string posStr;
	// handle first letter in string1
	if (index == 0) {
		posStr += 'S';
		for (int j = 0; j < patternLen - 1; j++) {
			if (operations[index + j].CDR == "D") {
				posStr += NOT_IN_LCS;
			}
			else { // CDR == "C" or "R"
				posStr += operations[index + j].opParam;
			}
		}
		return posStr;
	}

	int string1Size = operations.size() - 1;
	int lettersToEnd = string1Size - index;
	int maxPatternLen = ((patternLen - 1) < lettersToEnd) ? patternLen - 1 : lettersToEnd;
	for (int j = -1; j < maxPatternLen; j++) {

		if (operations[index + j].CDR == "D") {
			posStr += NOT_IN_LCS;
		}
		else { // CDR == "C" or "R"
			posStr += operations[index + j].opParam;
		}
	}

	if (patternLen - 1 > lettersToEnd) {
		posStr += 'S';
	}
	return posStr;
}

vector<string> RepStringRDIArray(const vector<LetterOps>& operations, const int patternLen) {
	unsigned vectorSize = operations.size();
	vector<string> repStrings(vectorSize);
	for (unsigned index = 0; index < vectorSize; index++) {
		//vector<LetterOps> noSubstitutions = RemoveSubstitutions(operations);
		repStrings[index] = RepStringRDI(operations, GapsRDI(operations), patternLen, index);
	}
	return repStrings;
}

vector<string> RepStringRDISubstitutionsArray(const vector<LetterOps>& operations, const int patternLen) {
	unsigned vectorSize = operations.size();
	vector<string> repStrings(vectorSize);
	for (unsigned index = 0; index < vectorSize; index++) {
		repStrings[index] = RepStringRDISubstitutions(operations, patternLen, index);
	}
	return repStrings;
}

vector<string> RepStringRDISubstitutionsArrayNoGaps(const vector<LetterOps>& operations, const int patternLen) {
	unsigned vectorSize = operations.size();
	vector<string> repStrings(vectorSize);
	for (unsigned index = 0; index < vectorSize; index++) {
		repStrings[index] = RepStringRDISubstitutionsNoGaps(operations, patternLen, index);
	}
	return repStrings;
}

vector<map<string, int> > Cluster2::RepStringRDIDictionaries(const int cloneIndex, const int patternLen,
		const int priority, mt19937& generator) {
	unsigned vectorSize = clones[cloneIndex].Len() + 1;
	vector<map<string, int> > dictionaries(vectorSize);
	ComputeEditLists(cloneIndex, priority, generator);
	for (unsigned listIndex = 0; listIndex < editLists.size(); listIndex++) {
		vector<string> currentVector = RepStringRDIArray(editLists[listIndex], patternLen);
		for (unsigned letterIndex = 0; letterIndex < vectorSize; letterIndex++) {
			dictionaries[letterIndex][currentVector[letterIndex]]++;
		}
	}
	return dictionaries;
}
vector<map<string, int> > Cluster2::RepStringRDISubDictionaries(const int cloneIndex, const int patternLen,
		const int priority, mt19937& generator) {
	unsigned vectorSize = clones[cloneIndex].Len() + 1;
	vector<map<string, int> > dictionaries(vectorSize);
	ComputeEditLists(cloneIndex, priority, generator);
	for (unsigned listIndex = 0; listIndex < editLists.size(); listIndex++) {
		vector<string> currentVector = RepStringRDISubstitutionsArray(editLists[listIndex], patternLen);
		for (unsigned letterIndex = 0; letterIndex < vectorSize; letterIndex++) {
			dictionaries[letterIndex][currentVector[letterIndex]]++;
		}
	}
	return dictionaries;
}

vector<map<string, int> > Cluster2::RepStringRDISubNoGapsDictionaries(const int cloneIndex, const int patternLen,
		const int priority, mt19937& generator) {
	unsigned vectorSize = clones[cloneIndex].Len() + 1;
	vector<map<string, int> > dictionaries(vectorSize);
	ComputeEditLists(cloneIndex, priority, generator);
	for (unsigned listIndex = 0; listIndex < editLists.size(); listIndex++) {
		vector<string> currentVector = RepStringRDISubstitutionsArrayNoGaps(editLists[listIndex], patternLen);
		for (unsigned letterIndex = 0; letterIndex < vectorSize; letterIndex++) {
			dictionaries[letterIndex][currentVector[letterIndex]]++;
		}
	}
	return dictionaries;
}
vector<map<string, int> > Cluster2::RepStringNumberGapsDictionaries(const int cloneIndex, const int patternLen) {
	unsigned vectorSize = clones[cloneIndex].Len() + 1;
	vector<map<string, int> > dictionaries(vectorSize);
	CmpLCS2(cloneIndex);
	for (vector<LCS2>::iterator it = lcs2s.begin(); it != lcs2s.end(); it++) {
		vector<string> currentVector = it->RepStringWithNumberGapArray(patternLen);
		for (unsigned int index = 0; index < vectorSize; index++) {
			dictionaries[index][currentVector[index]]++;
		}
	}
	return dictionaries;
}

vector<map<string, int> > Cluster2::RepStringNoGapsDictionaries(const int cloneIndex, const int patternLen) {
	unsigned vectorSize = clones[cloneIndex].Len() + 1;
	vector<map<string, int> > dictionaries(vectorSize);
	CmpLCS2(cloneIndex);
	for (vector<LCS2>::iterator it = lcs2s.begin(); it != lcs2s.end(); it++) {
		vector<string> currentVector = it->RepStringNoGapArray(patternLen);
		for (unsigned int index = 0; index < vectorSize; index++) {
			dictionaries[index][currentVector[index]]++;
		}
	}
	return dictionaries;
}

vector<int> Cluster2::GuessInsertByInsArray(const int cloneIndex, double insThreshold) {
	vector<int> inserts;
	vector<int> cum = InsArray(cloneIndex);
	for (int index = 0; index < (int) cum.size(); index++) {
		if ((double) cum[index] < (double) lcs2s.size() * insThreshold) {
			inserts.push_back(index);
		}
	}
	return inserts;
}

string max(const map<string, int>& m) {
	string maxKey = "";
	int maxValue = 0;
	for (map<string, int>::const_iterator it = m.begin(); it != m.end(); it++) {
		if (it->second > maxValue) {
			maxKey = it->first;
			maxValue = it->second;
		}
	}
	return maxKey;
}
string min(const map<string, int>& m) {
	string minKey = "";
	int minValue = 1000000;
	for (map<string, int>::const_iterator it = m.begin(); it != m.end(); it++) {
		if (it->second < minValue) {
			minKey = it->first;
			minValue = it->second;
		}
	}
	return minKey;
}

string SecondUppercase(const string& str) {
	int secondIndex = 1;
	while (not isupper(str[secondIndex])) {
		secondIndex++;
	}
	string result = str.substr(secondIndex, 1);
	return result;
}

vector<int> BuildInsertGuessesWithGaps(const int patternLen, const vector<string>& chosenStrings) {
	vector<int> guesses;
	for (unsigned int index = 0; index < chosenStrings.size() - 1; index++) {
		if (SecondUppercase(chosenStrings[index]) == NOT_IN_LCS) {
			guesses.push_back(index);
		}
	}
	return guesses;
}

vector<int> BuildInsertGuessesNoGaps(const int patternLen, const vector<string>& chosenStrings) {
	vector<int> guesses;
	for (unsigned int index = 0; index < chosenStrings.size() - 1; index++) {
		if (chosenStrings[index][1] == NOT_IN_LCS[0]) {
			guesses.push_back(index);
		}
	}
	return guesses;
}

vector<int> Cluster2::GuessInsertsByNumberGapsStringMax(const int cloneIndex, const int patternLen) {
	vector<map<string, int> > dictionariesArray = RepStringNumberGapsDictionaries(cloneIndex, patternLen);
	vector<string> chosenStrings;
	for (unsigned int index = 0; index < dictionariesArray.size(); index++) {
		chosenStrings.push_back(max(dictionariesArray[index]));
	}
	return BuildInsertGuessesWithGaps(patternLen, chosenStrings);
}

vector<int> Cluster2::GuessInsertsByNumberGapsStringPath(const int cloneIndex, const int patternLen) {
	vector<map<string, int> > dictionariesArray = RepStringNumberGapsDictionaries(cloneIndex, patternLen);
	vector<string> chosenStrings = DictLongestPath(dictionariesArray, patternLen);
	return BuildInsertGuessesWithGaps(patternLen, chosenStrings);
}

vector<int> Cluster2::GuessInsertsByLowercaseGapsStringMax(const int cloneIndex, const int patternLen) {
	vector<map<string, int> > dictionariesArray = RepStringLowercaseGapsDictionaries(cloneIndex, patternLen);
	vector<string> chosenStrings;
	for (unsigned int index = 0; index < dictionariesArray.size(); index++) {
		chosenStrings.push_back(max(dictionariesArray[index]));
	}
	return BuildInsertGuessesWithGaps(patternLen, chosenStrings);
}

vector<int> Cluster2::GuessInsertsByLowercaseGapsStringPath(const int cloneIndex, const int patternLen) {
	vector<map<string, int> > dictionariesArray = RepStringLowercaseGapsDictionaries(cloneIndex, patternLen);
	vector<string> chosenStrings = DictLongestPath(dictionariesArray, patternLen);
	return BuildInsertGuessesWithGaps(patternLen, chosenStrings);
}

vector<int> Cluster2::GuessInsertsByNoGapsStringMax(const int cloneIndex, const int patternLen) {
	vector<map<string, int> > dictionariesArray = RepStringNoGapsDictionaries(cloneIndex, patternLen);
	vector<string> chosenStrings;
	for (unsigned int index = 0; index < dictionariesArray.size(); index++) {
		chosenStrings.push_back(max(dictionariesArray[index]));
	}
	return BuildInsertGuessesNoGaps(patternLen, chosenStrings);
}

vector<int> Cluster2::GuessInsertsByNoGapsStringPath(const int cloneIndex, const int patternLen) {
	vector<map<string, int> > dictionariesArray = RepStringNoGapsDictionaries(cloneIndex, patternLen);
	vector<string> chosenStrings = DictLongestPath(dictionariesArray, patternLen);
	return BuildInsertGuessesNoGaps(patternLen, chosenStrings);
}

string FirstGap(const string& repString) {
	int strIndex = 2;
	while (islower(repString[strIndex])) {
		strIndex++;
	}
	string result = repString.substr(1, strIndex - 1);
	return result;
}

vector<pair<int, string> > BuildDeletionGuesses(const vector<string>& chosenStrings) {
	vector<pair<int, string> > guesses;
	for (unsigned int index = 0; index < chosenStrings.size(); index++) {
		if (not isupper(chosenStrings[index][1])) {
			string gap = FirstGap(chosenStrings[index]);
			transform(gap.begin(), gap.end(), gap.begin(), ::toupper);
			guesses.push_back(make_pair(index, gap));
		}
	}

	return guesses;
}

vector<pair<int, string> > Cluster2::GuessDeletionsByLowercaseGapsStringMax(const int cloneIndex,
		const int patternLen) {
	vector<map<string, int> > dictionariesArray = RepStringLowercaseGapsDictionaries(cloneIndex, patternLen);
	vector<string> chosenStrings;
	for (unsigned int index = 0; index < dictionariesArray.size(); index++) {
		// TODO: adjust max function for deletions as in python
		chosenStrings.push_back(max(dictionariesArray[index]));
	}
	return BuildDeletionGuesses(chosenStrings);
}

vector<pair<int, string> > Cluster2::GuessDeletionsByLowercaseGapsStringPath(const int cloneIndex,
		const int patternLen) {
	vector<map<string, int> > dictionariesArray = RepStringLowercaseGapsDictionaries(cloneIndex, patternLen);
	vector<string> chosenStrings = DictLongestPath(dictionariesArray, patternLen);
	return BuildDeletionGuesses(chosenStrings);
}

vector<pair<int, string> > Cluster2::GuessDeletionsByRDIStringPath(const int cloneIndex, const int patternLen,
		const int priority, mt19937& generator) {
	vector<map<string, int> > dictionariesArray = RepStringRDIDictionaries(cloneIndex, patternLen, priority, generator);
	vector<string> chosenStrings = DictLongestPath(dictionariesArray, patternLen);
	return BuildDeletionGuesses(chosenStrings);
}

vector<pair<int, string> > Cluster2::GuessSubstitutions(const int cloneIndex, const int priority, mt19937& generator) {
	vector<pair<int, string> > substitutions;
	vector<map<string, int> > dictionariesArray = CumEditOperationsNoInserts(cloneIndex, priority, generator);
	for (unsigned int index = 0; index < dictionariesArray.size(); index++) {
		string maxOp = max(dictionariesArray[index]);
		if (maxOp.empty()) {
			continue;
		}
		if (maxOp[0] == 'R') {
			substitutions.push_back(make_pair(index, string(1, maxOp[1])));
		}
	}
	return substitutions;
}

vector<int> Cluster2::GuessInsertsRDI(const int cloneIndex, const int priority, mt19937& generator) {
	vector<int> inserts;
	vector<map<string, int> > dictionariesArray = CumEditOperationsNoInserts(cloneIndex, priority, generator);
	for (unsigned int index = 0; index < dictionariesArray.size(); index++) {
		string maxOp = max(dictionariesArray[index]);
		if (maxOp.empty()) {
			continue;
		}
		if (maxOp[0] == 'D') {
			inserts.push_back(index);
		}
	}
	return inserts;
}

vector<pair<int, string> > Cluster2::GuessDeletionsRDI(const int cloneIndex, const int priority, mt19937& generator) {
	vector<pair<int, string> > deletions;
	vector<map<string, int> > dictionariesArray = CumEditOperationsInserts(cloneIndex, priority, generator);
	for (unsigned int index = 0; index < dictionariesArray.size(); index++) {
		string maxString = max(dictionariesArray[index]);
		if (not maxString.empty()) {
			deletions.push_back(make_pair(index, maxString));
		}
	}
	return deletions;
}

vector<pair<int, string> > BuildSubstitutionGuesses(const int patternLen, const vector<string>& chosenStrings,
		const string& copy) {
	vector<pair<int, string> > guesses;
	for (unsigned int index = 0; index < chosenStrings.size() - 1; index++) {
		string secondUppercase = SecondUppercase(chosenStrings[index]);
		if (secondUppercase != string(1, copy[index]) && secondUppercase != NOT_IN_LCS) {
			guesses.push_back(make_pair(index, secondUppercase));
		}
	}
	return guesses;
}

vector<pair<int, string> > BuildSubstitutionGuessesNoGaps(const int patternLen, const vector<string>& chosenStrings,
		const string& copy) {
	vector<pair<int, string> > guesses;
	for (unsigned int index = 0; index < chosenStrings.size() - 1; index++) {
		string secondUppercase(1, chosenStrings[index][1]);
		if (secondUppercase != string(1, copy[index]) && secondUppercase != NOT_IN_LCS) {
			guesses.push_back(make_pair(index, secondUppercase));
		}
	}
	return guesses;
}

vector<pair<int, string> > Cluster2::GuessSubstitutionsRDIPath(const int cloneIndex, const int patternLen,
		const int priority, mt19937& generator) {
	vector<map<string, int> > dictionariesArray = RepStringRDISubDictionaries(cloneIndex, patternLen, priority,
			generator);
	vector<string> chosenStrings = DictLongestPath(dictionariesArray, patternLen);
	return BuildSubstitutionGuesses(patternLen, chosenStrings, clones[cloneIndex].String());
}

vector<pair<int, string> > Cluster2::GuessSubstitutionsRDINoGapsPath(const int cloneIndex, const int patternLen,
		const int priority, mt19937& generator) {
	vector<map<string, int> > dictionariesArray = RepStringRDISubNoGapsDictionaries(cloneIndex, patternLen, priority,
			generator);
	vector<string> chosenStrings = DictLongestPath(dictionariesArray, patternLen);
	return BuildSubstitutionGuessesNoGaps(patternLen, chosenStrings, clones[cloneIndex].String());
}

void Cluster2::TestInsertGuesses(const int cloneIndex, const double insThreshold, const int patternLen,
		int* trueInsertGuesses, int* falseInsertGuesses, int* totalInserts, int* guessesNum) {
	int oldInsertNum = clones[cloneIndex].GetInsertNum();
	int oldDeleteNum = clones[cloneIndex].GetDeleteNum();
	//FixAllInserts(cloneIndex, patternLen);
	//vector<int> guesses = GuessInsertByInsArray(cloneIndex, insThreshold);
	//vector<int> guesses = GuessInsertsByLowercaseGapsStringMax(cloneIndex, patternLen);
	//vector<int> guesses = GuessInsertsByLowercaseGapsStringPath(cloneIndex, patternLen);
	//vector<int> guesses = GuessInsertsByNumberGapsStringMax(cloneIndex, patternLen);
	//vector<int> guesses = GuessInsertsByNumberGapsStringPath(cloneIndex, patternLen);
	//vector<int> guesses = GuessInsertsByNoGapsStringMax(cloneIndex, patternLen);
	vector<int> guesses = GuessInsertsByNoGapsStringPath(cloneIndex, patternLen);

	clones[cloneIndex].FixInsertGuesses(guesses);
	int newInsertNum = clones[cloneIndex].GetInsertNum();
	int newDeleteNum = clones[cloneIndex].GetDeleteNum();

	int insertNumDecrease = oldInsertNum - newInsertNum;
	assert(insertNumDecrease >= 0);
	int deletionNumIncrease = newDeleteNum - oldDeleteNum;

	*trueInsertGuesses = insertNumDecrease;
	*falseInsertGuesses = max(deletionNumIncrease, 0);
	*totalInserts = oldInsertNum;
	*guessesNum = guesses.size();
}

void Cluster2::TestDeletionGuesses(const int cloneIndex, const int patternLen, int* trueDeletionGuesses,
		int* falseDeletionGuesses, int* totalDeletions, int* guessesNum) {
	int oldInsertNum = clones[cloneIndex].GetInsertNum();
	int oldDeleteNum = clones[cloneIndex].GetDeleteNum();
	//FixDeletionsAllExcept(cloneIndex, patternLen);
	vector<pair<int, string> > guesses = GuessDeletionsByLowercaseGapsStringMax(cloneIndex, patternLen);
	//vector<pair<int, string> > guesses = GuessDeletionsByLowercaseGapsStringPath(cloneIndex, patternLen);

	clones[cloneIndex].FixDeletionGuesses(guesses);
	int newInsertNum = clones[cloneIndex].GetInsertNum();
	int newDeleteNum = clones[cloneIndex].GetDeleteNum();

	int insertNumIncrease = newInsertNum - oldInsertNum;
	assert(insertNumIncrease >= 0);
	int deletionNumDecrease = oldDeleteNum - newDeleteNum;
	assert(deletionNumDecrease >= 0);

	*trueDeletionGuesses = deletionNumDecrease;
	*falseDeletionGuesses = insertNumIncrease;
	*totalDeletions = oldDeleteNum;
	int countGuesses = 0;
	for (vector<pair<int, string> >::iterator it = guesses.begin(); it != guesses.end(); it++) {
		countGuesses += it->second.size();
	}
	*guessesNum = countGuesses;
}

void Cluster2::FixAllDeletions(const int patternLen) {
	for (int index = 0; index < (int) clones.size(); index++) {
		vector<pair<int, string> > guesses = GuessDeletionsByLowercaseGapsStringPath(index, patternLen);
		clones[index].FixDeletionGuesses(guesses);
	}
}

void Cluster2::FixAllInserts(const int insPriority, mt19937& generator) {
	for (int index = 0; index < (int) clones.size(); index++) {
		vector<int> guesses = GuessInsertsRDI(index, insPriority, generator);
		clones[index].FixInsertGuesses(guesses);
	}
}

void Cluster2::FixAllSubstitutions(const int subPriority, mt19937& generator) {
	for (int index = 0; index < (int) clones.size(); index++) {
		vector<pair<int, string> > guesses = GuessSubstitutions(index, subPriority, generator);
		clones[index].FixSubstitutionGuesses(guesses);
	}
}

// Fix deletions of all clones by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.

void Cluster2::FixAllDeletionsClonesLastInitial(const int delPatternLen, const int delPriority, mt19937& generator) {
	vector<Clone> computedClones, beforeFixClones = clones;
	for (int index = 0; index < (int) clones.size(); index++) {
		//clones[index] = clonesBackup[index];
		//vector<pair<int, string> > guesses = GuessDeletionsRDI(index, delPriority, generator);
		vector<pair<int, string> > guesses = GuessDeletionsByLowercaseGapsStringPath(index, delPatternLen);
		clones[index].FixDeletionGuesses(guesses);
		computedClones.push_back(clones[index]);
		clones[index] = beforeFixClones[index];
	}
	clones = computedClones;
}

//	Fix insertions of all clones by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.

void Cluster2::FixAllInsertsClonesLastInitial(const int insPriority, mt19937& generator) {
	vector<Clone> computedClones, beforeFixClones = clones;
	for (int index = 0; index < (int) clones.size(); index++) {
		//clones[index] = clonesBackup[index];
		vector<int> guesses = GuessInsertsRDI(index, insPriority, generator);
		clones[index].FixInsertGuesses(guesses);
		computedClones.push_back(clones[index]);
		clones[index] = beforeFixClones[index];
	}
	clones = computedClones;
}

// Fix substitutions of all clones by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.

void Cluster2::FixAllSubstitutionsClonesLastInitial(const int subPriority, mt19937& generator) {
	vector<Clone> computedClones, beforeFixClones = clones;
	for (int index = 0; index < (int) clones.size(); index++) {
		//clones[index] = clonesBackup[index];
		vector<pair<int, string> > guesses = GuessSubstitutions(index, subPriority, generator);
		clones[index].FixSubstitutionGuesses(guesses);
		computedClones.push_back(clones[index]);
		clones[index] = beforeFixClones[index];
	}
	clones = computedClones;
}

// Fix one clone. fix order: substitutions, deletions with STRING PATH, insertions.

void Cluster2::FixOneCopyRDIDelPath(const int cloneIndex, const int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator) {

	vector<pair<int, string> > guesses = GuessSubstitutions(cloneIndex, subPriority, generator);
//	if (cloneIndex == 0 and not guesses.empty()) {
//		cout << "Substitutions Fixed!" << endl;
//	}
	clones[cloneIndex].FixSubstitutionGuesses(guesses);

	vector<pair<int, string> > guesses1 = GuessDeletionsByLowercaseGapsStringPath(cloneIndex, delPatternLen);
//	vector<pair<int, string> > guesses1 = GuessDeletionsRDI(cloneIndex, delPriority, generator);
//	if (cloneIndex == 0 and not guesses1.empty()) {
//		cout << "Deletions Fixed!" << endl;
//	}
	clones[cloneIndex].FixDeletionGuesses(guesses1);

	vector<int> guesses2 = GuessInsertsRDI(cloneIndex, insPriority, generator);
//	if (cloneIndex == 0 and not guesses2.empty() and guesses1.empty()) {
//		cout << "Inserts Fixed without previous deletion fixed!" << endl;
//	}
//	if (cloneIndex == 0 and not guesses2.empty() and not guesses1.empty()) {
//		cout << "Inserts Fixed after deletion fixed!" << endl;
//	}
	clones[cloneIndex].FixInsertGuesses(guesses2);

}

// Fix one clone. fix order: substitutions, deletions RDI, insertions.

void Cluster2::FixOneCopyRDIAll(const int cloneIndex, const int subPatternLen, const int delPatternLen,
		const double insThreshold, const int subPriority, const int delPriority, const int insPriority,
		mt19937& generator) {

	vector<pair<int, string> > guesses = GuessSubstitutions(cloneIndex, subPriority, generator);
	clones[cloneIndex].FixSubstitutionGuesses(guesses);

	vector<pair<int, string> > guesses1 = GuessDeletionsRDI(cloneIndex, delPriority, generator);
	clones[cloneIndex].FixDeletionGuesses(guesses1);

	vector<int> guesses2 = GuessInsertsRDI(cloneIndex, insPriority, generator);
	clones[cloneIndex].FixInsertGuesses(guesses2);

}

vector<pair<int, string> > FilterDeletionGuesses(const vector<pair<int, string> >& guesses, const int maxIndex) {
	vector<pair<int, string> > filteredGuesses;
	for (unsigned int i = 0; i < guesses.size(); i++) {
		if (guesses[i].first >= maxIndex - (int) filteredGuesses.size()) {
			filteredGuesses.push_back(guesses[i]);
		}
	}
	return filteredGuesses;
}

vector<int> FilterInsertGuesses(const vector<int>& guesses, const int maxIndex) {
	vector<int> filteredGuesses;
	for (unsigned int i = 0; i < guesses.size(); i++) {
		if (guesses[i] >= maxIndex + (int) filteredGuesses.size()) {
			filteredGuesses.push_back(guesses[i]);
		}
	}
	return filteredGuesses;
}

// Fix all one by one: Fix one clone. fix order: substitutions, deletions with STRING PATH, insertions.

void Cluster2::FixAllOneByOneRDI(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator) {
	for (int index = 0; index < (int) clones.size(); index++) {
		FixOneCopyRDIDelPath(index, delPatternLen, subPriority, delPriority, insPriority, generator);
	}
}

//	Fix all clones one at a time. change clones only after all fixes computed.
//	Each clone: fix order: substitutions, deletions with STRING PATH, insertions.
//	All fixed clones will be stored in clones vector.

void Cluster2::FixAllOneByOneRDIChangeClonesLast(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator) {
	vector<Clone> computedClones, beforeFixClones = clones;
	for (int index = 0; index < (int) clones.size(); index++) {
		FixOneCopyRDIDelPath(index, delPatternLen, subPriority, delPriority, insPriority, generator);
		computedClones.push_back(clones[index]);
		clones[index] = beforeFixClones[index];

	}
	clones = computedClones;
}

//	Fix all clones one at a time. before fix restore clone to initial state in backup. change clones only after all fixes computed.
//	Each clone: fix order: substitutions, deletions with STRING PATH, insertions.
//	All fixed clones will be stored in clones vector.

void Cluster2::FixAllOneByOneRDIChangeClonesLastInitial(const int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator) {
	vector<Clone> computedClones, beforeFixClones = clones;
	for (int index = 0; index < (int) clones.size(); index++) {
		clones[index] = clonesBackup[index];
		FixOneCopyRDIDelPath(index, delPatternLen, subPriority, delPriority, insPriority, generator);
		computedClones.push_back(clones[index]);
		clones[index] = beforeFixClones[index];
	}
	clones = computedClones;
}

int CountCorrectSize(const vector<Clone>& clones, unsigned correctSize) {
	int correctSizeNum = 0;
	for (unsigned i = 0; i < clones.size(); i++) {
		if (clones[i].String().size() == correctSize) {
			correctSizeNum++;
		}
	}
	return correctSizeNum;
}

// Repeat maxReps times:
// Fix all one by one: Fix one clone. fix order: substitutions, deletions with STRING PATH, insertions.

void Cluster2::FixAllOneByOneRDIRepeat(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator, const int maxReps) {

	for (int i = 0; i < maxReps; i++) {
		FixAllOneByOneRDI(delPatternLen, subPriority, delPriority, insPriority, generator);
	}
}

//	Repeat maxReps times.
//	Fix all clones one at a time. change clones only after all fixes computed.
//	Each clone: fix order: substitutions, deletions with STRING PATH, insertions.
//	All fixed clones will be stored in clones vector.

void Cluster2::FixAllOneByOneRDIChangeClonesLastRepeat(const int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int minEqual, const int maxReps) {

	for (int i = 0; i < maxReps; i++) {
		FixAllOneByOneRDIChangeClonesLast(delPatternLen, subPriority, delPriority, insPriority, generator);
	}
}

//	Repeat maxReps times:
//	Fix all clones one at a time. before fix restore clone to initial state in backup. change clones only after all fixes computed.
//	Each clone: fix order: substitutions, deletions, insertions.
//	All fixed clones will be stored in clones vector.

void Cluster2::FixAllOneByOneRDIChangeClonesLastInitialRepeat(const int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int maxReps) {

	for (int i = 0; i < maxReps; i++) {
		FixAllOneByOneRDIChangeClonesLastInitial(delPatternLen, subPriority, delPriority, insPriority, generator);
	}
}

void Cluster2::TestSubstitutionGuesses(const int clone_index, map<string, double>& roundCountBefore,
		map<string, double>& roundCountAfter, double& distNoReplaceBefore, double& distNoReplaceAfter,
		const int priority, mt19937& generator) {

	FixAllSubstitutions(priority, generator);

	vector<LetterOps> opListAfter = ComputeEditDistancePriority(original, clones[clone_index].String(), priority,
			generator);
	roundCountAfter = CountOperations(opListAfter);
	distNoReplaceAfter = ComputeEditDistWithoutReplace(clones[clone_index].String(), original);
}

int CumCloneToClonesED(const int cloneIndex, const vector<Clone>& clones) {
	int clonesNum = clones.size();
	int cumCloneED = 0;
	for (int i = 0; i < clonesNum; i++) {
		if (i != cloneIndex) {
			cumCloneED += ComputeEditDistanceNum(clones[i].String(), clones[cloneIndex].String());
		}
	}
	return cumCloneED;
}
vector<int> AllCumCloneToClonesED(const vector<Clone>& clones) {
	int clonesNum = clones.size();
	vector<int> clonesED(clonesNum);
	for (int i = 0; i < clonesNum; i++) {
		for (int j = i + 1; j < clonesNum; j++) {
			int edij = ComputeEditDistanceNum(clones[i].String(), clones[j].String());
			clonesED[i] += edij;
			clonesED[j] += edij;
		}
	}
	return clonesED;
}

map<string, int> UniqueClonesMap(const vector<Clone>& clones) {
	int clonesNum = clones.size();
	map<string, int> clonesMap;
	for (int i = 0; i < clonesNum; i++) {
		clonesMap[clones[i].String()]++;
	}
	return clonesMap;
}

string MapCorrectSizeMaxString(const vector<Clone>& clones, const int correctSize) {
	string maxString;
	int maxNum = 0;
	map<string, int> clonesMap = UniqueClonesMap(clones);
	for (map<string, int>::iterator it = clonesMap.begin(); it != clonesMap.end(); it++) {
		if ((int) it->first.size() == correctSize and it->second > maxNum) {
			maxNum = it->second;
			maxString = it->first;
		}
	}
	return maxString;
}

int CorrectSizeCloneCount(const vector<Clone>& clones, const int correctSize) {
	int clonesNum = clones.size();
	int correctSizeNum = 0;
	for (int i = 0; i < clonesNum; i++) {
		if (clones[i].Len() == correctSize) {
			correctSizeNum++;
		}
	}
	return correctSizeNum;
}

void Stats4(const vector<Clone>& clones, const string& original, vector<int>& originalEDS, int& cumCloneOriginalED,
		int& cumCorrectSizeCloneOriginalED, int& correctSizeNum, int& bingoCount, int& firstHalfBingoCount,
		int& secondHalfBingoCount, int& minED, int& maxED) {
	originalEDS.clear();
	cumCloneOriginalED = 0;
	cumCorrectSizeCloneOriginalED = 0;
	correctSizeNum = 0;
	bingoCount = 0;
	firstHalfBingoCount = 0;
	secondHalfBingoCount = 0;
	minED = original.size();
	maxED = 0;
	int originalSize = original.size();
	string original1stHalf = original.substr(0, originalSize / 2);
	string original2ndHalf = original.substr(originalSize / 2);
	string cloneStr, clone1stHalf, clone2ndHalf;
	for (unsigned int i = 0; i < clones.size(); i++) {
		int editDist = ComputeEditDistanceNum(original, clones[i].String());
		originalEDS.push_back(editDist);
		if (editDist == 0) {
			bingoCount++;
		}
		if (editDist < minED) {
			minED = editDist;
		}
		if (editDist > maxED) {
			maxED = editDist;
		}
		cumCloneOriginalED += editDist;
		if (clones[i].Len() == (int) original.size()) {
			cloneStr = clones[i].String();
			clone1stHalf = cloneStr.substr(0, originalSize / 2);
			clone2ndHalf = cloneStr.substr(originalSize / 2);
			if (clone1stHalf == original1stHalf) {
				firstHalfBingoCount++;
			}
			if (clone2ndHalf == original2ndHalf) {
				secondHalfBingoCount++;
			}
			cumCorrectSizeCloneOriginalED += editDist;
			correctSizeNum++;
		}
	}
}

char Max(const map<char, int>& count) {
	char maxLetter = 0;
	int maxValue = 0;
	for (map<char, int>::const_iterator it = count.begin(); it != count.end(); it++) {
		if (it->second > maxValue) {
			maxValue = it->second;
			maxLetter = it->first;
		}
	}
	return maxLetter;
}

string LetterwiseMajority(const vector<string>& strings) {
	int stringsLen = strings[0].size();
	map<char, int> count;
	string majority;
	for (int i = 0; i < stringsLen; i++) {
		count.clear();
		for (int j = 0; j < (int) strings.size(); j++) {
			count[strings[j][i]]++;
		}
		char maxLetter = Max(count);
		majority += maxLetter;
	}
	return majority;
}

string Majority(const vector<string>& strings) {
	int stringsLen = strings[0].size();
	map<string, int> count;
	for (int j = 0; j < (int) strings.size(); j++) {
		string firstHalf = strings[j].substr(0, stringsLen / 2);
		count[firstHalf]++;
	}
	string maxFirstHalf = max(count);
	count.clear();
	for (int j = 0; j < (int) strings.size(); j++) {
		string secondHalf = strings[j].substr(stringsLen / 2);
		count[secondHalf]++;
	}
	string maxSecondHalf = max(count);
	return maxFirstHalf + maxSecondHalf;
}

string BestClone(const vector<Clone>& clones, int correctSize) {
//	return clones[0].String();
	return clones[clones.size() - 1].String();

	vector<int> cumCloneToClonesED = AllCumCloneToClonesED(clones);

	int minValue = cumCloneToClonesED[0];
	string minString = clones[0].String();
	for (unsigned i = 0; i < clones.size(); i++) {
		if (cumCloneToClonesED[i] < minValue) {
			minValue = cumCloneToClonesED[i];
			minString = clones[i].String();
		}
		else if (cumCloneToClonesED[i] == minValue) {
			int minStrLenDif = abs((int) minString.size() - correctSize);
			int candidateLenDif = abs(clones[i].Len() - correctSize);
			if (candidateLenDif < minStrLenDif) {
				minValue = cumCloneToClonesED[i];
				minString = clones[i].String();
			}
		}
	}
	return minString;
}

void Cluster2::SortClonesByInitial(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator) {
	FixAllOneByOneRDIChangeClonesLast(delPatternLen, subPriority, delPriority, insPriority, generator);
	//compare c;
	compare c;
	int clonesNum = clones.size();
	vector<pair<int, int> > cumED(clonesNum);
	for (int i = 0; i < clonesNum; i++) {
		cumED[i].first = CumCloneToClonesED(i, clones);
		cumED[i].second = i;
	}

	sort(cumED.begin(), cumED.end(), c);
	vector<Clone> sortedClones;
	for (int i = 0; i < clonesNum; i++) {
		int oldIndex = cumED[i].second;
		sortedClones.push_back(clonesBackup[oldIndex]);
	}
	clones = sortedClones;
	clonesBackup = clones;
}

// Repeat maxReps times:
// Fix all one by one: Fix one clone. fix order: substitutions, deletions with STRING PATH, insertions.
// Return fixed clones in normalResult.
// Reset clones to original state in backup

void Cluster2::FixAllNormal(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator, const int maxReps, vector<Clone>& normalResult,
		vector<Clone>& reverseResult) {
	FixAllOneByOneRDIRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, maxReps);
	normalResult = clones;
	clones = clonesBackup;
}

set<string> CorrectSizeSet(const vector<Clone>& clones, const unsigned correctSize) {
	set<string> result;
	for (unsigned i = 0; i < clones.size(); i++) {
		if (clones[i].String().size() == correctSize) {
			result.insert(clones[i].String());
		}
	}
	return result;
}

void NearClonesMaps(const vector<Clone>& clones, const unsigned correctSize, map<string, int>& correctSizeClones,
		map<string, int>& oneLetterDiffClones) {
	correctSizeClones.clear();
	oneLetterDiffClones.clear();
	for (unsigned i = 0; i < clones.size(); i++) {
		if (clones[i].String().size() == correctSize) {
			correctSizeClones[clones[i].String()]++;
		}
		else if (clones[i].String().size() == correctSize + 1 or clones[i].String().size() == correctSize - 1) {
			oneLetterDiffClones[clones[i].String()]++;
		}
	}
}

string MaxSet(const set<string>& allGuesses, const vector<Clone>& originalClones) {
	if (allGuesses.size() == 1) {
		return *allGuesses.begin();
	}
	vector<string> vec(allGuesses.begin(), allGuesses.end());
	map<string, int> cumEDS;
	for (unsigned i = 0; i < vec.size(); i++) {
		for (unsigned j = 0; j < originalClones.size(); j++) {
			int editDistance = ComputeEditDistanceNum(vec[i], originalClones[j].String());
			cumEDS[vec[i]] += editDistance;
		}
	}
	return min(cumEDS);
}

string FinalGuess(const vector<Clone>& clones, const vector<Clone>& clonesBackup, const int correctSize) {
	map<string, int> correctSizeClones, oneLetterDiffClones;
	set<string> correctSizeAll = CorrectSizeSet(clones, correctSize);
	string finalGuess;
	// count occurrences in normal and reverse together and guess the most frequent.
	NearClonesMaps(clones, correctSize, correctSizeClones, oneLetterDiffClones);
	if (not correctSizeClones.empty()) {
		finalGuess = MaxSet(correctSizeAll, clonesBackup);
		//finalGuess = max(correctSizeClones);
	}
	else if (not oneLetterDiffClones.empty()) {
		finalGuess = max(oneLetterDiffClones);
	}
	else {
		finalGuess = clones[0].String();
	}
	return finalGuess;
}

int FindIndex(const vector<Clone>& clones, const string& guess) {
	for (unsigned i = 0; i < clones.size(); i++) {
		if (clones[i].String() == guess) {
			return i;
		}
	}
	assert(0);
	return 0;
}

int ConvGuess(const vector<Clone>& clones, const vector<Clone>& clonesBackup, const int correctSize) {
//	map<string, int> correctSizeClones, oneLetterDiffClones;
//	set<string> correctSizeAll = CorrectSizeSet(clones, correctSize);
//	string finalGuess;
//	// count occurrences in normal and reverse together and guess the most frequent.
//	NearClonesMaps(clones, correctSize, correctSizeClones, oneLetterDiffClones);
//	if (not correctSizeClones.empty()) {
//		finalGuess = MaxSet(correctSizeAll, clonesBackup);
//	}
//	else if (not oneLetterDiffClones.empty()) {
//		finalGuess = max(oneLetterDiffClones);
//	}
//	else {
//		finalGuess = clones[0].String();
//	}
//	return FindIndex(clones, finalGuess);
	int maxLenDiff = abs(correctSize - clones[0].Len());
	int maxLenIndex = 0;
	for (unsigned i = 1; i < clones.size(); i++) {
		int currLenDiff = abs(correctSize - clones[i].Len());
		if (currLenDiff > maxLenDiff) {
			maxLenDiff = currLenDiff;
			maxLenIndex = i;
		}
	}
	return maxLenIndex;
}

void Cluster2::TestFixAll(const int delPatternLen, vector<int>& cloneOriginalEDS, double& avgCloneEditDist,
		int& cumCorrectSizeCloneEditDist, int& correctSizeCloneNum, int& bingoNum, int& firstHalfBingoNum,
		int& secondHalfBingoNum, int& minED, int& maxED, int& finalGuessEditDist, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int maxReps) {

	vector<Clone> normalResult, reverseResult, joinedResult;
	FixAllNormal(delPatternLen, subPriority, delPriority, insPriority, generator, maxReps, normalResult, reverseResult);
	// join normal clones and reverse clones
	joinedResult = normalResult;
	joinedResult.insert(joinedResult.end(), reverseResult.begin(), reverseResult.end());

	int cumCloneOriginalED, correctSizeNum, bingoCount, firstHalfBingoCount, secondHalfBingoCount,
			cumCorrectSizeCloneOriginalED;

	Stats4(joinedResult, original, cloneOriginalEDS, cumCloneOriginalED, cumCorrectSizeCloneOriginalED, correctSizeNum,
			bingoCount, firstHalfBingoCount, secondHalfBingoCount, minED, maxED);
	avgCloneEditDist = (double) cumCloneOriginalED / joinedResult.size();
	cumCorrectSizeCloneEditDist = cumCorrectSizeCloneOriginalED;
	correctSizeCloneNum = correctSizeNum;
	bingoNum = bingoCount;
	firstHalfBingoNum = firstHalfBingoCount;
	secondHalfBingoNum = secondHalfBingoCount;

	string finalGuess = FinalGuess(joinedResult, clonesBackup, original.size());

	finalGuessEditDist = ComputeEditDistanceNum(original, finalGuess);
}

//	Repeat maxReps times.
//	Fix all clones one at a time. change clones only after all fixes computed.
//	Each clone: fix order: substitutions, deletions with STRING PATH, insertions.
//	All fixed clones will be stored in normalResult vector.
//	Restore clones to initial state in backup.

void Cluster2::FixAllClonesLastNormal(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator, const int minEqual, const int maxReps, vector<Clone>& normalResult,
		vector<Clone>& reverseResult) {
	FixAllOneByOneRDIChangeClonesLastRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, minEqual,
			maxReps);
	normalResult = clones;
	clones = clonesBackup;
}

void Cluster2::TestFixAllClonesLast(const int delPatternLen, vector<int>& cloneOriginalEDS, double& avgCloneEditDist,
		int& cumCorrectSizeCloneEditDist, int& correctSizeCloneNum, int& bingoNum, int& firstHalfBingoNum,
		int& secondHalfBingoNum, int& minED, int& maxED, int& finalGuessEditDist, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int minEqual, const int maxReps) {

	vector<Clone> normalResult, reverseResult, joinedResult;
	FixAllClonesLastNormal(delPatternLen, subPriority, delPriority, insPriority, generator, minEqual, maxReps,
			normalResult, reverseResult);
	// join normal clones and reverse clones
	joinedResult = normalResult;
	joinedResult.insert(joinedResult.end(), reverseResult.begin(), reverseResult.end());

	int cumCloneOriginalED, correctSizeNum, bingoCount, firstHalfBingoCount, secondHalfBingoCount,
			cumCorrectSizeCloneOriginalED;

	Stats4(joinedResult, original, cloneOriginalEDS, cumCloneOriginalED, cumCorrectSizeCloneOriginalED, correctSizeNum,
			bingoCount, firstHalfBingoCount, secondHalfBingoCount, minED, maxED);
	avgCloneEditDist = (double) cumCloneOriginalED / joinedResult.size();
	cumCorrectSizeCloneEditDist = cumCorrectSizeCloneOriginalED;
	correctSizeCloneNum = correctSizeNum;
	bingoNum = bingoCount;
	firstHalfBingoNum = firstHalfBingoCount;
	secondHalfBingoNum = secondHalfBingoCount;

	string finalGuess = FinalGuess(joinedResult, clonesBackup, original.size());

	finalGuessEditDist = ComputeEditDistanceNum(original, finalGuess);
}

void NearClonesMaps(const vector<Clone>& normalResult, const vector<Clone>& reverseResult, const unsigned correctSize,
		map<string, int>& correctSizeClones, map<string, int>& oneLetterDiffClones) {
	correctSizeClones.clear();
	oneLetterDiffClones.clear();
	for (unsigned i = 0; i < normalResult.size(); i++) {
		if (normalResult[i].String().size() == correctSize) {
			correctSizeClones[normalResult[i].String()]++;
		}
		else if (normalResult[i].String().size() == correctSize + 1
				or normalResult[i].String().size() == correctSize - 1) {
			oneLetterDiffClones[normalResult[i].String()]++;
		}
	}
	for (unsigned i = 0; i < reverseResult.size(); i++) {
		if (reverseResult[i].String().size() == correctSize) {
			correctSizeClones[reverseResult[i].String()]++;
		}
		else if (reverseResult[i].String().size() == correctSize + 1
				or reverseResult[i].String().size() == correctSize - 1) {
			oneLetterDiffClones[reverseResult[i].String()]++;
		}
	}
}

//	Vertical stage:
//	Fix substitutions in all clones, then deletions and then insertions.

//	Horizontal stage:
//	Fix all clones one at a time. before fix restore clone to initial state in backup.
//	Each clone: fix order: substitutions, deletions, insertions.

//	Store fixed clones in normalResult.
//	Reset clones to original state in backup.

void Cluster2::FixAllClonesVerticalHorizontalNormal(const int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int maxReps,
		vector<Clone>& normalResult, vector<Clone>& reverseResult) {
	FixAllSubstitutionsClonesLastInitial(subPriority, generator);
	FixAllDeletionsClonesLastInitial(delPatternLen, delPriority, generator);
	FixAllInsertsClonesLastInitial(insPriority, generator);

	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, 1);
	normalResult = clones;
	clones = clonesBackup;
}

//	Normal phase:
//	Vertical stage:
//	Fix substitutions in all clones, then deletions and then insertions.

//	Horizontal stage:
//	Fix all clones one at a time. before fix restore clone to initial state in backup.
//	Each clone: fix order: substitutions, deletions, insertions.

//	Store fixed clones in normalResult.
//	Reset clones to original state in backup.

//	Reverse phase:
//	Reverse all clones. Do same as in normal phase.
//	Store fixed clones in reverseResult
//	Reset clones to original backup state.

void Cluster2::FixAllClonesVerticalHorizontalNormalAndReverse(const int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int maxReps,
		vector<Clone>& normalResult, vector<Clone>& reverseResult) {
	vector<Clone> temp = clonesBackup;

	FixAllSubstitutionsClonesLastInitial(subPriority, generator);
	FixAllDeletionsClonesLastInitial(delPatternLen, delPriority, generator);
	FixAllInsertsClonesLastInitial(insPriority, generator);
	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, 1);
	normalResult = clones;
	clones = clonesBackup;
	ReverseClones();
	clonesBackup = clones;

	FixAllSubstitutionsClonesLastInitial(subPriority, generator);
	FixAllDeletionsClonesLastInitial(delPatternLen, delPriority, generator);
	FixAllInsertsClonesLastInitial(insPriority, generator);
	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, 1);

	ReverseClones();
	reverseResult = clones;
	// restore originals
	clonesBackup = temp;
	clones = clonesBackup;
}

//	Horizontal stage:
//	Fix all clones one at a time. before fix restore clone to initial state in backup.
//	Each clone: fix order: substitutions, deletions, insertions.

//	Vertical stage:
//	Fix substitutions in all clones, then deletions and then insertions.

//	Store fixed clones in normalResult.
//	Reset clones to original state in backup.

void Cluster2::FixAllClonesHorizontalVerticalNormal(const int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int maxReps,
		vector<Clone>& normalResult, vector<Clone>& reverseResult) {
	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, 1);

	FixAllSubstitutionsClonesLastInitial(subPriority, generator);
	FixAllDeletionsClonesLastInitial(delPatternLen, delPriority, generator);
	FixAllInsertsClonesLastInitial(insPriority, generator);

	normalResult = clones;
	clones = clonesBackup;
}

//	Normal phase:

//	Horizontal stage:
//	Fix all clones one at a time. before fix restore clone to initial state in backup.
//	Each clone: fix order: substitutions, deletions, insertions.

//	Vertical stage:
//	Fix substitutions in all clones, then deletions and then insertions.

//	Store fixed clones in normalResult.
//	Reset clones to original state in backup.

//	Reverse phase:
//	Reverse all clones. Do same as in normal phase.
//	Store fixed clones in reverseResult
//	Reset clones to original backup state.

void Cluster2::FixAllClonesHorizontalVerticalNormalAndReverse(const int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int maxReps,
		vector<Clone>& normalResult, vector<Clone>& reverseResult) {
	vector<Clone> temp = clonesBackup;

	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, 1);
	FixAllSubstitutionsClonesLastInitial(subPriority, generator);
	FixAllDeletionsClonesLastInitial(delPatternLen, delPriority, generator);
	FixAllInsertsClonesLastInitial(insPriority, generator);

	normalResult = clones;
	clones = clonesBackup;
	ReverseClones();
	clonesBackup = clones;

	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, 1);
	FixAllSubstitutionsClonesLastInitial(subPriority, generator);
	FixAllDeletionsClonesLastInitial(delPatternLen, delPriority, generator);
	FixAllInsertsClonesLastInitial(insPriority, generator);

	ReverseClones();
	reverseResult = clones;
	// restore originals
	clonesBackup = temp;
	clones = clonesBackup;
}

//	Repeat maxReps times:
//	Fix all clones one at a time. before fix restore clone to initial state in backup.
//	Each clone: fix order: substitutions, deletions, insertions.
//	All fixed clones will be stored in normalResult.
//	clones will be reset to original backup state.

void Cluster2::FixAllClonesLastInitialNormal(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator, const int maxReps, vector<Clone>& normalResult,
		vector<Clone>& reverseResult) {
	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator,
			maxReps);
	normalResult = clones;
	clones = clonesBackup;
}

//	Normal stage:
//	Repeat maxReps times:
//	Fix all clones one at a time. before fix restore clone to initial state in backup.
//	Each clone: fix order: substitutions, deletions, insertions.
//	store fixed clones in normalResult.
//	clones will be reset to original backup state.

//	Reverse stage:
//	Reverse all clones. Do same as in normal stage.
//	Store fixed clones in reverseResult
//	Reset clones to original backup state.

void Cluster2::FixAllClonesLastInitialNormalAndReverse(const int delPatternLen, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int maxReps,
		vector<Clone>& normalResult, vector<Clone>& reverseResult) {
	vector<Clone> temp = clonesBackup;
	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator,
			maxReps);
	normalResult = clones;
	clones = clonesBackup;
	ReverseClones();
	clonesBackup = clones;
	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator,
			maxReps);
	ReverseClones();
	reverseResult = clones;
	// restore originals
	clonesBackup = temp;
	clones = clonesBackup;
}

void Cluster2::TestFixAllClonesLastInitial(const int delPatternLen, vector<int>& cloneOriginalEDS,
		double& avgCloneEditDist, int& cumCorrectSizeCloneEditDist, int& correctSizeCloneNum, int& bingoNum,
		int& firstHalfBingoNum, int& secondHalfBingoNum, int& minED, int& maxED, int& finalGuessEditDist,
		const int subPriority, const int delPriority, const int insPriority, mt19937& generator, const int maxReps) {

	vector<Clone> normalResult, reverseResult, joinedResult;
	FixAllClonesLastInitialNormal(delPatternLen, subPriority, delPriority, insPriority, generator, maxReps,
			normalResult, reverseResult);
	// join normal clones and reverse clones
	joinedResult = normalResult;
	joinedResult.insert(joinedResult.end(), reverseResult.begin(), reverseResult.end());

	int cumCloneOriginalED, correctSizeNum, bingoCount, firstHalfBingoCount, secondHalfBingoCount,
			cumCorrectSizeCloneOriginalED;

	Stats4(joinedResult, original, cloneOriginalEDS, cumCloneOriginalED, cumCorrectSizeCloneOriginalED, correctSizeNum,
			bingoCount, firstHalfBingoCount, secondHalfBingoCount, minED, maxED);
	avgCloneEditDist = (double) cumCloneOriginalED / joinedResult.size();
	cumCorrectSizeCloneEditDist = cumCorrectSizeCloneOriginalED;
	correctSizeCloneNum = correctSizeNum;
	bingoNum = bingoCount;
	firstHalfBingoNum = firstHalfBingoCount;
	secondHalfBingoNum = secondHalfBingoCount;
	string finalGuess = FinalGuess(joinedResult, clonesBackup, original.size());

	finalGuessEditDist = ComputeEditDistanceNum(original, finalGuess);
}

//	Fix substitutions of all clones by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.
//	Fix deletion of all clones by string path by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.
//	Fix insertions of all clones by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.

void Cluster2::FixAllByErrorTypeRepeat(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator, const int maxReps) {
	for (int i = 0; i < maxReps; i++) {
		FixAllSubstitutionsClonesLastInitial(subPriority, generator);
		FixAllDeletionsClonesLastInitial(delPatternLen, delPriority, generator);
		FixAllInsertsClonesLastInitial(insPriority, generator);
	}
}

//	Fix substitutions of all clones by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.
//	Fix deletion of all clones by string path by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.
//	Fix insertions of all clones by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones returned in nomralResult. Reset clones to initial state in backup

void Cluster2::FixAllByErrorTypeNormal(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator, const int maxReps, vector<Clone>& normalResult,
		vector<Clone>& reverseResult) {
	FixAllByErrorTypeRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, maxReps);
	normalResult = clones;
	clones = clonesBackup;
}

//	Normal phase:
//	Fix substitutions of all clones by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.
//	Fix deletion of all clones by string path by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones stored in clones.
//	Fix insertions of all clones by the current state of clones. CHANGE TO BACKUP IS COMMENTED
//	Fixed clones returned in nomralResult. Reset clones to initial state in backup

//	Reverse phase:
//	Reverse clones. Do same as above on reversed. store fixed clones in reverseResult.
//	Reset clones to initial state in backup

void Cluster2::FixAllByErrorTypeNormalAndReverse(const int delPatternLen, const int subPriority, const int delPriority,
		const int insPriority, mt19937& generator, vector<Clone>& normalResult, vector<Clone>& reverseResult) {
	vector<Clone> temp = clonesBackup;
	FixAllSubstitutionsClonesLastInitial(subPriority, generator);
	FixAllDeletionsClonesLastInitial(delPatternLen, delPriority, generator);
	FixAllInsertsClonesLastInitial(insPriority, generator);
	normalResult = clones;

	// reverse clones and backup
	clones = clonesBackup;
	ReverseClones();
	clonesBackup = clones;

	FixAllSubstitutions(subPriority, generator);
	FixAllDeletions(delPatternLen);
	FixAllInserts(insPriority, generator);
	ReverseClones();
	reverseResult = clones;

	// restore originals
	clonesBackup = temp;
	clones = clonesBackup;
}

void Cluster2::TestFixAllByErrorType(const int delPatternLen, const double insThreshold, vector<int>& cloneOriginalEDS,
		double& avgCloneEditDist, int& cumCorrectSizeCloneEditDist, int& correctSizeCloneNum, int& bingoNum,
		int& firstHalfBingoNum, int& secondHalfBingoNum, int& minED, int& maxED, int& finalGuessEditDist,
		const int subPriority, const int delPriority, const int insPriority, mt19937& generator, const int minEqual,
		const int maxReps, const int convReps) {

	vector<Clone> normalResult, reverseResult, joinedResult;
//	FixAllByErrorTypeConverge(subPatternLen, delPatternLen, subPriority, delPriority, insPriority, generator, maxReps,
//			convReps, normalResult, reverseResult);
	FixAllByErrorTypeNormal(delPatternLen, subPriority, delPriority, insPriority, generator, maxReps, normalResult,
			reverseResult);
	// join normal clones and reverse clones
	joinedResult = normalResult;
	joinedResult.insert(joinedResult.end(), reverseResult.begin(), reverseResult.end());

	int cumCloneOriginalED, correctSizeNum, bingoCount, firstHalfBingoCount, secondHalfBingoCount,
			cumCorrectSizeCloneOriginalED;

	Stats4(joinedResult, original, cloneOriginalEDS, cumCloneOriginalED, cumCorrectSizeCloneOriginalED, correctSizeNum,
			bingoCount, firstHalfBingoCount, secondHalfBingoCount, minED, maxED);
	avgCloneEditDist = (double) cumCloneOriginalED / joinedResult.size();
	cumCorrectSizeCloneEditDist = cumCorrectSizeCloneOriginalED;
	correctSizeCloneNum = correctSizeNum;
	bingoNum = bingoCount;
	firstHalfBingoNum = firstHalfBingoCount;
	secondHalfBingoNum = secondHalfBingoCount;
	string finalGuess = FinalGuess(joinedResult, clonesBackup, original.size());

	finalGuessEditDist = ComputeEditDistanceNum(original, finalGuess);
}

//	compute final guess by one method. add to clones. fix with another method.

void Cluster2::FixAllMixed(const int delPatternLen, const int subPriority, const int delPriority, const int insPriority,
		mt19937& generator, const int maxReps) {

	vector<Clone> normalClones, reverseClones;
	FixAllOneByOneRDIChangeClonesLastInitialRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, 2);
//	FixAllByErrorTypeRepeat(subPatternLen, delPatternLen, subPriority, delPriority, insPriority, generator, maxReps);

	string finalGuess = FinalGuess(clones, clonesBackup, original.size());
	Clone newClone(original, finalGuess);
	clones = clonesBackup;
	clones.insert(clones.begin(), newClone);
	FixAllOneByOneRDIRepeat(delPatternLen, subPriority, delPriority, insPriority, generator, 2);

}

void Cluster2::TestFixAllMixed(const int delPatternLen, vector<int>& cloneOriginalEDS, double& avgCloneEditDist,
		int& cumCorrectSizeCloneEditDist, int& correctSizeCloneNum, int& bingoNum, int& firstHalfBingoNum,
		int& secondHalfBingoNum, int& minED, int& maxED, int& finalGuessEditDist, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int maxReps) {

	FixAllMixed(delPatternLen, subPriority, delPriority, insPriority, generator, maxReps);

	map<string, int> correctSizeClones, oneLetterDiffClones;
	int cumCloneOriginalED, correctSizeNum, bingoCount, firstHalfBingoCount, secondHalfBingoCount,
			cumCorrectSizeCloneOriginalED;
	Stats4(clones, original, cloneOriginalEDS, cumCloneOriginalED, cumCorrectSizeCloneOriginalED, correctSizeNum,
			bingoCount, firstHalfBingoCount, secondHalfBingoCount, minED, maxED);
	avgCloneEditDist = (double) cumCloneOriginalED / clones.size();
	cumCorrectSizeCloneEditDist = cumCorrectSizeCloneOriginalED;
	correctSizeCloneNum = correctSizeNum;
	bingoNum = bingoCount;
	firstHalfBingoNum = firstHalfBingoCount;
	secondHalfBingoNum = secondHalfBingoCount;

	string finalGuess = FinalGuess(clones, clonesBackup, original.size());
	//string finalGuess = clones[0].String();
	finalGuessEditDist = ComputeEditDistanceNum(original, finalGuess);
}

//	Best Algorithm so far:
//	fix clones first initial RDI NORMAL AND REVERSE. 2 rounds.
//	fix by error type NORMAL AND REVERSE 1 round.
//	join results.

string Cluster2::TestBest(const int delPatternLen, int& finalGuessEditDist, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator, const int maxReps) {
	vector<Clone> normalResult, reverseResult, joinedResult;
	FixAllClonesLastInitialNormalAndReverse(delPatternLen, subPriority, delPriority, insPriority, generator, maxReps,
			normalResult, reverseResult);
	joinedResult = normalResult;
	joinedResult.insert(joinedResult.end(), reverseResult.begin(), reverseResult.end());
	FixAllByErrorTypeNormalAndReverse(delPatternLen, subPriority, delPriority, insPriority, generator, normalResult,
			reverseResult);
	joinedResult.insert(joinedResult.end(), normalResult.begin(), normalResult.end());
	joinedResult.insert(joinedResult.end(), reverseResult.begin(), reverseResult.end());

//	map<string, int> correctSizeClones, oneLetterDiffClones;
//	int cumCloneOriginalED, correctSizeNum, bingoCount, firstHalfBingoCount, secondHalfBingoCount,
//			cumCorrectSizeCloneOriginalED;
//	Stats4(joinedResult, original, cloneOriginalEDS, cumCloneOriginalED, cumCorrectSizeCloneOriginalED, correctSizeNum,
//			bingoCount, firstHalfBingoCount, secondHalfBingoCount, minED, maxED);
//	avgCloneEditDist = (double) cumCloneOriginalED / clones.size();
//	cumCorrectSizeCloneEditDist = cumCorrectSizeCloneOriginalED;
//	correctSizeCloneNum = correctSizeNum;
//	bingoNum = bingoCount;
//	firstHalfBingoNum = firstHalfBingoCount;
//	secondHalfBingoNum = secondHalfBingoCount;

	string finalGuess = FinalGuess(joinedResult, clonesBackup, original.size());
	//string finalGuess = clones[0].String();
//	finalGuessEditDist = ComputeEditDistanceNum(original, finalGuess);
	return finalGuess;
}

void Cluster2::TestOriginalRetention(const int delPatternLen, int& finalGuessEditDist, const int subPriority,
		const int delPriority, const int insPriority, mt19937& generator) {
	Clone org(original, original);
	clones[0] = org;
	clonesBackup[0] = org;
//	FixAllOneByOneRDIChangeClonesLastInitial(subPatternLen, delPatternLen, insThreshold, subPriority, delPriority,
//			insPriority, generator);
	FixOneCopyRDIDelPath(0, delPatternLen, subPriority, delPriority, insPriority, generator);
//	FixAllOneByOneRDIChangeClonesLast(subPatternLen, delPatternLen, insThreshold, subPriority, delPriority,
//				insPriority, generator);
//	FixAllSubstitutionsClonesLastInitial(subPatternLen, subPriority, generator);
//	FixAllDeletionsClonesLastInitial(delPatternLen);
//	FixAllInsertsClonesLastInitial(insPriority, generator);

	string finalGuess = clones[0].String();
	finalGuessEditDist = ComputeEditDistanceNum(original, finalGuess);
}

void Cluster2::CloneStats(vector<int>& lengths, vector<int>& originalED, vector<int>& clonesED) const {
	int clonesNum = clones.size();
	lengths = vector<int>(clonesNum);
	originalED = vector<int>(clonesNum);
	clonesED = vector<int>(clonesNum);
	for (int index = 0; index < clonesNum; index++) {
		lengths[index] = clones[index].Len();
		originalED[index] = ComputeEditDistanceNum(original, clones[index].String());
	}
	for (int i = 0; i < clonesNum; i++) {
		for (int j = i + 1; j < clonesNum; j++) {
			int edij = ComputeEditDistanceNum(clones[i].String(), clones[j].String());
			clonesED[i] += edij;
			clonesED[j] += edij;
		}
	}
}

void Cluster2::PrintStats() const {
	int clonesNum = clones.size();
	vector<int> lengths;
	vector<int> originalED;
	vector<int> clonesED;
	CloneStats(lengths, originalED, clonesED);
	for (int index = 0; index < clonesNum; index++) {
		cout << lengths[index] << "\t" << originalED[index] << "\t" << (double) clonesED[index] / (clonesNum - 1)
				<< endl;
	}
	cout << "*************************************************" << endl;
}

void Cluster2::Stats(const int delPatternLen, const int subPriority, const int delPriority, const int insPriority,
		mt19937& generator) {
	PrintStats();
	FixAllOneByOneRDIChangeClonesLast(delPatternLen, subPriority, delPriority, insPriority, generator);
	PrintStats();
	FixAllOneByOneRDIChangeClonesLast(delPatternLen, subPriority, delPriority, insPriority, generator);

	PrintStats();
	FixAllOneByOneRDIChangeClonesLast(delPatternLen, subPriority, delPriority, insPriority, generator);

	PrintStats();
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

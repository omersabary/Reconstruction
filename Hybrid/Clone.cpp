#include "Clone.hpp"
#include "LCS2.hpp"
#include <iostream>
#include <algorithm>
#include <map>
#include <chrono>
using namespace std;

string MakeStrand(const unsigned length, mt19937& generator) {
	string strand;
	uniform_int_distribution<int> distribution(0, 3);
	string letters("GTAC");
	for (unsigned i = 0; i < length; i++) {
		strand += letters[distribution(generator)];
	}
	return strand;
}

string CopyStrand(const string& source, mt19937& generator, const double delProb, const double insProb,
		const double subProb) {
	uniform_real_distribution<double> distributionReal(0.0, 1.0);
	uniform_int_distribution<int> distributionInt(0, 3);
	uniform_int_distribution<int> distributionInt3(0, 2);
	string letters("GTAC");
	map<char, string> substitutions = { { 'A', "CGT" }, { 'C', "AGT" }, { 'G', "ACT" }, { 'T', "ACG" } };
	string copy;
	for (size_t i = 0; i < source.size(); i++) {
		if (distributionReal(generator) <= insProb) { // chose to insert
			copy += letters[distributionInt(generator)]; // choose random letter from GTAC
		}
		if (distributionReal(generator) > delProb) { // don't delete
			if (distributionReal(generator) <= subProb) { // substitute
				copy += substitutions[source[i]][distributionInt3(generator)];
			}
			else {	// just copy
				copy += source[i];
			}
		}
		else { // delete
			if (distributionReal(generator) <= subProb) { // delete + substitute = substitute
				copy += substitutions[source[i]][distributionInt3(generator)];
			}
		}
	}
	// insert at the end
	if (distributionReal(generator) <= insProb) {
		copy += letters[distributionInt(generator)];
	}
	return copy;
}

Clone::Clone(const string& original, const string& clone) :
		original(original), clone(clone), insertNum(0), deleteNum(0) {
}

Clone::Clone(const string& original, const double delProb, const double insProb, const double subProb,
		std::mt19937& generator) :
		original(original), clone(), insertNum(0), deleteNum(0) {
	clone = CopyStrand(original, generator, delProb, insProb, subProb);
	UpdateInsertDelete();
}

bool Clone::operator<(const Clone& a) const {
	return clone < a.clone;
}
bool Clone::operator==(const Clone& a) const {
	return clone == a.clone;
}

size_t CloneHasher::operator()(const Clone& clone) const {
	return hash<string>()(clone.String());
}

string Clone::String() const {
	return clone;
}
int Clone::Len() const {
	return clone.length();
}

void Clone::UpdateInsertDelete() {
	LCS2 lcs = ComputeLCS2(original, clone);
	deleteNum = original.length() - lcs.Len();
	insertNum = clone.length() - lcs.Len();
}

int Clone::GetInsertNum() const {
	return insertNum;
}

int Clone::GetDeleteNum() const {
	return deleteNum;
}

void Clone::Reverse(){
	reverse(clone.begin(),clone.end());
}

//fix insert guesses by deleting them from copy
void Clone::FixInsertGuesses(const vector<int>& guesses) {
	for (int index = (int) guesses.size() - 1; index >= 0; index--) {
		clone.erase(guesses[index], 1);
	}
	UpdateInsertDelete();
}

void Clone::FixDeletionGuesses(const vector<pair<int, string> >& guesses) {
	for (int index = (int) guesses.size() - 1; index >= 0; index--) {
		int position = guesses[index].first;
		string deleted = guesses[index].second;
		clone.insert(position, deleted);
	}
	UpdateInsertDelete();
}

void Clone::FixSubstitutionGuesses(const vector<pair<int, string> >& guesses) {
	for (int index = (int) guesses.size() - 1; index >= 0; index--) {
		int position = guesses[index].first;
		string rep = guesses[index].second;
		clone.replace(position, 1, rep);
	}
}

ostream& operator<<(std::ostream& os, const Clone& a) {
	os << a.clone;
	return os;
}

void TestCopyStrand() {
	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);
	string strand = MakeStrand(20, generator);
	cout << CopyStrand(strand, generator, 0.1, 0.1, 0.1) << endl;
}

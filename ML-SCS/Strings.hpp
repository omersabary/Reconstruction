#ifndef STRINGS_HPP_
#define STRINGS_HPP_

#include <string>
#include <random>
#include <set>
#include <map>

std::string MakeRandomString4(const int len, std::mt19937& generator);
std::string CopyStrandWithDeletion(const std::string& source, std::mt19937& generator, double delProb);
void RandomStrings(std::vector<std::string>& strings, const int stringNum, const int stringLen,
		std::mt19937& generator);
void RandomStringsRandomLen(std::vector<std::string>& strings, const int stringNum, const int minStrLen,
		const int maxStrLen, std::mt19937& generator);
std::string RandomCopies2(std::vector<std::string>& copies, const int originalLen, const int copiesNum,
		const double delProb, std::mt19937& generator);
std::string RandomCopies4(std::vector<std::string>& copies, const int originalLen, const int copiesNum,
		const double delProb, std::mt19937& generator);
int CountEmbeddings(const std::string& a, const std::string& b);
std::vector<std::string> AllMostLikely(const std::string& X, const std::string& Y,
		const std::vector<std::string>& SCSs);
int EditDistance(const std::string& s, const std::string& t);
int EditDistanceDI(const std::string& s, const std::string& t);
std::string MostLikelySCS(const std::vector<std::string>& inputStrings, const long double p, const int originalLen,
		std::mt19937& generator);
bool HasPairWithCorrectSizeSCS(const std::vector<std::string>& inputStrings, const int originalLen);
int CorrectSizeKtupleNumSortByLetterMax(const std::vector<std::string>& inputStrings, const int k, const int originalLen);
int CorrectSizeKtupleNumSortBySumLen(const std::vector<std::string>& inputStrings, const int k, const int originalLen);
std::string MostLikelyFilteredSCS(const std::vector<std::string>& inputStrings, const int k, const long double p,
		const int originalLen);
std::string MinEditDistanceFilteredSCS(const std::vector<std::string>& inputStrings, const int k, const long double p,
		const int originalLen);
std::string MostLikelyCorrectSizeSCS(const std::vector<std::string>& inputStrings, const int k, const long double p,
		const int originalLen);
std::string MostLikelyCorrectSizeSCS(const std::vector<std::string>& inputStrings, const int k, const long double p,
		const int originalLen, int& foundKtupleIndex);
std::string MostLikelyMaxLenSCS(const std::vector<std::string>& inputStrings, const int k, const long double p,
		const int originalLen);
std::string MostLikelyMaxLenSCS(const std::vector<std::string>& inputStrings, const int k, const long double p,
		const int originalLen, bool& foundCorrectSize, int& foundKtupleIndex);
void MostLikelyEitan(const std::vector<std::string>& inputStrings, const int originalLen, int& scsNum,
		int& mostLikelyNum);
std::string MostLikelyEitanFull(const std::vector<std::string>& inputStrings, const int originalLen, int& scsNum,
		int& mostLikelyNum);
std::vector<std::string> AllMostLikelyEitan(const std::vector<std::string>& inputStrings, const int originalLen);
void MostLikelyNewAlgorithmVsOld(const std::vector<std::string>& inputStrings);
std::vector<std::map<char, int>> MaxLetterFrequency(const std::vector<std::string>& inputStrings,
		const int originalLen);
bool IsSubSequence(const std::string& str1, const std::string& str2);
std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& a);
std::ostream& operator<<(std::ostream& os, const std::set<std::string>& a);
std::ostream& operator<<(std::ostream& os, const std::vector<int>& a);
int LCS2Len(const std::string& X, const std::string& Y);
std::vector<std::string> LCS2(const std::string& X, const std::string& Y);
std::vector<std::string> LCS2DP(const std::string& X, const std::string& Y);
std::vector<std::vector<std::pair<int, int>>>LCS2Indexes(const std::string& X, const std::string& Y);
std::vector<std::string> LCS2Pointers(const std::string& X, const std::string& Y);
std::string Word(const std::vector<int>& indexes, const std::vector<char>& alphabet);
int min(const int a, const int b, const int c);
#endif /* STRINGS_HPP_ */

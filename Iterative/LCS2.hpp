#ifndef LCS2_HPP_
#define LCS2_HPP_
#include "CommonSubstring2.hpp"
#include <vector>
#include <set>

const std::string NOT_IN_LCS = "X";

class LCS2 {
	std::string string1;
	std::string string2;
	std::vector<CommonSubstring2> commonSubstrings;
	int commonSubstringsIndex;
	bool commonActive;
public:
	LCS2(const std::string& string1, const std::string& string2);
	void AddLetter(const char& letter, const int firstIndex, const int secondIndex);
	int Len() const;
	std::vector<std::vector<int> > Edges() const;
	std::vector<int> String1LCSIndexes() const;
	std::vector<int> String1Mirrors() const;
	std::vector<std::string> DelStringGaps(const std::vector<int>& string1Mirrors) const;
	std::string RepStringWithLowercaseGap(const int patternLen, const int index, const std::vector<int>& string1Mirrors,
			const std::vector<std::string>& string1MirrorsGaps) const;
	std::string RepStringWithNumberGap(const int patternLen, const int index, const std::vector<int>& string1Mirrors,
				const std::vector<std::string>& string1MirrorsGaps) const;
	std::string RepStringNoGap(const int patternLen, const int index, const std::vector<int>& string1Mirrors) const;
	std::vector<std::string> RepStringWithLowercaseGapArray(const int patternLen) const;
	std::vector<std::string> RepStringWithNumberGapArray(const int patternLen) const;
	std::vector<std::string> RepStringNoGapArray(const int patternLen) const;
	friend std::ostream& operator<<(std::ostream& os, const LCS2& a);
};

std::ostream& operator<<(std::ostream& os, const LCS2& a);
LCS2 ComputeLCS2(const std::string& X, const std::string& Y);
int ComputeEditDistWithoutReplace(const std::string& X, const std::string& Y);
#endif /* LCS2_HPP_ */

#ifndef MULTISTRING_HPP_
#define MULTISTRING_HPP_

#include <vector>
#include <string>

class MultiString {
	// sequence of vectors of strings.
	// sequence: each entry is a vector of strings:
	// 		for a string which is part of LCS: the vector will contain this string only
	//		UNIQUE strings which are not part of LCS: vector will contain them.
	// isLCS: each entry indicates if the corresponding vector in sequence is part of LCS or not.

	// example: the Common supersequences AB-XYZ-CD, AB-YXZ-CD will be represented by the multi string:
	// sequence[0]: ["AB"]			isLCS[0]: true
	// sequence[1]: ["XYZ","YXZ"]	isLCS[1]: false
	// sequence[2]: ["CD"]			isLCS[2]: true

	std::vector<std::vector<std::string>> sequence;
	std::vector<bool> isLCS;

public:
	void AddLetterLCS(const char& letter);
	void AddLetterNONLCS(const char& letter);
	void AddStringNONLCS(const std::string& str);
	bool empty() const;
	void Print() const;
	friend MultiString Merge(const MultiString& a, const MultiString& b);
};

std::vector<MultiString> MergeVectors(const std::vector<MultiString>& a, const std::vector<MultiString>& b);

#endif /* MULTISTRING_HPP_ */

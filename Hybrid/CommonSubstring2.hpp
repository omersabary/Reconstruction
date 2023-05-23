#ifndef COMMONSUBSTRING2_HPP_
#define COMMONSUBSTRING2_HPP_

#include <string>
#include <vector>
class CommonSubstring2{
	std::string substring;
	int startIn1;
	int endIn1;
	int startIn2;
	int endIn2;
public:
	CommonSubstring2(const char& letter, const int firstIndex, const int secondIndex);
	void AddLetter(const char& letter);
	void StartEnd(int* start1, int* end1, int* start2) const;
	std::pair<int, int> Range1() const;
	int Len() const;
	friend std::ostream& operator<<(std::ostream& os, const CommonSubstring2& a);
};

std::ostream& operator<<(std::ostream& os, const CommonSubstring2& a);

#endif /* COMMONSUBSTRING2_HPP_ */

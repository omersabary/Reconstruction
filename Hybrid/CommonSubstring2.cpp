#include "CommonSubstring2.hpp"
#include <iostream>
using namespace std;
CommonSubstring2::CommonSubstring2(const char& letter, const int firstIndex,
		const int secondIndex) :
		substring(1, letter), startIn1(firstIndex), endIn1(firstIndex + 1), startIn2(
				secondIndex), endIn2(secondIndex + 1) {
}

void CommonSubstring2::AddLetter(const char& letter) {
	substring += letter;
	endIn1++;
	endIn2++;
}

void CommonSubstring2::StartEnd(int* start1, int* end1, int* start2) const {
	*start1 = startIn1;
	*end1 = endIn1;
	*start2 = startIn2;
	return;
}

pair<int, int> CommonSubstring2::Range1() const{
	return make_pair(startIn1, endIn1);
}

int CommonSubstring2::Len() const{
	return substring.length();
}

std::ostream& operator<<(std::ostream& os, const CommonSubstring2& a) {
	os << a.substring;
	return os;
}

#ifndef CLONE_HPP_
#define CLONE_HPP_
#include <string>
#include <random>
#include <list>

class Clone {
	std::string original;
	std::string clone;
	int insertNum;
	int deleteNum;
	// list of operations
	// an operation: pair<char, char>
	// insert:		'I', letter inserted
	// deletion:	'D', letter deleted
	// copy: 		'C', letter copied
	// std::vector<std::pair<char, char> > operations;
	// trace a letter of the copy back to the operation
	// std::vector<int> letterToOperation;
public:
	Clone(const std::string& original, const std::string& clone);
	Clone(const std::string& original, const double delProb, const double insProb, const double subProb,
			std::mt19937& generator);
	bool operator<(const Clone& a) const;
	bool operator==(const Clone& a) const;
	std::string String() const;
	int Len() const;
	void UpdateInsertDelete();
	void Reverse();
	int GetInsertNum() const;
	int GetDeleteNum() const;
	void FixInsertGuesses(const std::vector<int>& guesses);
	void FixDeletionGuesses(const std::vector<std::pair<int, std::string> >& guesses);
	void FixSubstitutionGuesses(const std::vector<std::pair<int, std::string> >& guesses);
	friend std::ostream& operator<<(std::ostream& os, const Clone& a);
};

struct CloneHasher {
	size_t operator()(const Clone& clone) const;
};

std::string MakeStrand(const unsigned length, std::mt19937& generator);
std::string CopyStrand(const std::string& source, std::mt19937& generator, const double delProb, const double insProb, const double subProb);
void TestCopyStrand();

#endif /* CLONE_HPP_ */

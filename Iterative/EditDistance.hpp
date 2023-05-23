/*
 * EditDistance.hpp
 *
 *  Created on: 31 αιεμι 2019
 *      Author: Alex
 */

#ifndef EDITDISTANCE_HPP_
#define EDITDISTANCE_HPP_
#include <string>
#include <vector>
#include <map>
#include <random>
struct LetterOps {
	std::string insert; // string to insert before letter
	std::string CDR; // Copy or Delete or Substitute letter
	std::string opParam; // Which letter to copy or which to delete or which to substitute with
};
std::vector<LetterOps> ComputeEditDistancePriority(const std::string& X, const std::string& Y, const int priority,
		std::mt19937& generator);
std::vector<LetterOps> ComputeEditDistancePriorityReverse(const std::string& X, const std::string& Y,
		const int priority, std::mt19937& generator);
int ComputeEditDistanceNum(const std::string& X, const std::string& Y);
std::map<std::string, double> CountOperations(const std::vector<LetterOps>& opList);
std::vector<LetterOps> RemoveSubstitutions(const std::vector<LetterOps>& ops);

#endif /* EDITDISTANCE_HPP_ */

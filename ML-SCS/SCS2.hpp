
#ifndef SCS2_HPP_
#define SCS2_HPP_

#include <vector>
#include <string>
#include <random>
#include "MultiString.hpp"

std::vector<std::string> SCS2(const std::string& X, const std::string& Y);
std::string SCS2WithWildCards(const std::string& X, const std::string& Y);
std::string SCS2WithWildCardsRandom(const std::string& X, const std::string& Y, std::mt19937& generator);
int EstimateSCSNum(const std::string& scsWithWildCards);
//std::vector<std::string> MostLikelySCS2(const std::string& X, const std::string& Y);
int SCS2Len(const std::string& X, const std::string& Y);
std::vector<std::string> SCS2Fast(const std::string& X, const std::string& Y, const int originalLen);
std::vector<std::string> SCS2Fastest(const std::string& X, const std::string& Y, const int originalLen);
int SCS2LenFast(const std::string& X, const std::string& Y, const int originalLen);

#endif /* SCS2_HPP_ */

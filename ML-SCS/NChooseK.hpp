#ifndef NCHOOSEK_HPP_
#define NCHOOSEK_HPP_
#include <vector>
#include <string>
std::vector<std::vector<int> > Combinations(const int n, const int k);
std::vector<std::pair<std::string, std::string> > StringPairings(const std::vector<std::string>& strings);
std::vector<std::vector<std::string> > StringKtuples(const std::vector<std::string>& strings, const int k);
std::vector<std::vector<std::string> > SortedStringKtuplesBySumLen(const std::vector<std::string>& strings, const int k);
std::vector<std::vector<std::string> > SortedStringKtuplesBySumLetterMax(const std::vector<std::string>& strings, const int k);
#endif /* NCHOOSEK_HPP_ */

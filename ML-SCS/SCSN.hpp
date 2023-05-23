#ifndef SCSN_HPP_
#define SCSN_HPP_

#include <string>
#include "IndexVector.hpp"

int SCSNLen(const std::vector<std::string>& strings);
std::vector<std::string> SCSN(const std::vector<std::string>& strings);
int SCSNLenFast(const std::vector<std::string>& strings, const int originalLen);
int SCSNLenFastLessMem(const std::vector<std::string>& strings, const int originalLen);
std::vector<std::string> SCSNFast(const std::vector<std::string>& strings, const int originalLen);
#endif /* SCSN_HPP_ */

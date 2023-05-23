//
// Created by User on 26/03/2020.
//

#ifndef DNA_RECONSTRUCTION_INDEXVECTORD_H
#define DNA_RECONSTRUCTION_INDEXVECTORD_H
using namespace std;
#include <vector>


class IndexVectorD {
    int dimNum;
    int maxMemorySize;
    int absoluteIndex;
    std::vector<int> indexes;
    const std::vector<int> dimSizes;
    std::vector<int> dimSizePowers;
    const std::vector<int> D;
    std::vector<int> lowerBounds;
    std::vector<int> upperBounds;
public:
    IndexVectorD(const std::vector<int>& maxIndexes, const std::vector<int>& D);
    void operator++(int);
    bool IsEnd() const;
    const std::vector<int>& Vector() const;
    int MaxMemorySize() const;
    void UpdateAbsoluteIndex();
    int AbsInd() const;
    int AbsIndBeta(const std::vector<int>& K) const;
    std::vector<int> BetaIndexesOnly(const std::vector<int>& K) const;
    IndexVectorD BetaClass(const std::vector<int>& K) const;
    IndexVectorD& Last();
};

#endif //DNA_RECONSTRUCTION_INDEXVECTORD_H

#ifndef INDEXVECTOR_HPP_
#define INDEXVECTOR_HPP_

#include <vector>

class IndexVector {
	int dimNum;
	std::vector<int> dimSizes;
	std::vector<int> dimSizePowers;
	std::vector<int> indexVector;
	int absoluteArrayIndex;
	int arraySize;

public:
	IndexVector(const std::vector<int>& dims);
	void operator++(int);
	IndexVector Beta(const std::vector<int>& K);
	IndexVector& Last();
	std::vector<int> Vector() const {
		return indexVector;
	}
	bool IsPastEnd() const {
		return absoluteArrayIndex < 0 or arraySize == 0;
	}
	int AbsInd() const {
		return absoluteArrayIndex;
	}
	int ArraySize() const {
		return arraySize;
	}
};

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

#endif /* INDEXVECTOR_HPP_ */

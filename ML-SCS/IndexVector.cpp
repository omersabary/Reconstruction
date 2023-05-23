#include "IndexVector.hpp"
#include <iostream>
#include <cassert>
#include <algorithm>
using namespace std;

IndexVector::IndexVector(const std::vector<int>& dims) :
		dimNum(dims.size()), dimSizes(dims), dimSizePowers(dimNum), indexVector(dimNum), absoluteArrayIndex(0), arraySize(
				0) {
	if (dimNum == 0) {
		absoluteArrayIndex = -1;
		return;
	}
	dimSizePowers[dimNum - 1] = 1;
	for (int i = dimNum - 2; i >= 0; i--) {
		dimSizePowers[i] = dimSizePowers[i + 1] * dimSizes[i + 1];
	}

	arraySize = dimSizePowers[0] * dimSizes[0];
	//cout << "Array dimension:\t" << arraySize << endl;
}

void IndexVector::operator++(int) {
	if (absoluteArrayIndex == arraySize - 1) { // we are at array end
		absoluteArrayIndex = -1; // signal past end vector
	}
	else { // we are not at array end
		absoluteArrayIndex++;
		int currentDimIndex = dimNum - 1;
		//set to zero all indexes which are at max
		while (indexVector[currentDimIndex] == dimSizes[currentDimIndex] - 1) {
			indexVector[currentDimIndex] = 0;
			currentDimIndex--;
		}
		// increment first non max index
		indexVector[currentDimIndex]++;
	}
	return;
}

IndexVector IndexVector::Beta(const vector<int>& K) {
	IndexVector result = *this;
	for (vector<int>::const_iterator itvec = K.begin(); itvec != K.end(); itvec++) {
		int currentDimIndex = *itvec;
		assert(result.indexVector[currentDimIndex] > 0);
		result.indexVector[currentDimIndex]--;
		result.absoluteArrayIndex -= dimSizePowers[currentDimIndex];
	}
	return result;
}

IndexVector& IndexVector::Last() {
	absoluteArrayIndex = arraySize - 1;
	for (int i = 0; i < dimNum; i++) {
		indexVector[i] = dimSizes[i] - 1;
	}
	return *this;
}

// IndexVector D

IndexVectorD::IndexVectorD(const vector<int>& dims, const vector<int>& D) :
		dimNum(dims.size()), maxMemorySize(0), absoluteIndex(0), indexes(dimNum), dimSizes(dims), dimSizePowers(dimNum), D(
				D), lowerBounds(dimNum), upperBounds(dimNum) {
	// Start vector is [0,0...]. Always in range.
	// calculate upper bounds for all indexes zero
	vector<int> limits(D);
	// limits for dim[0]
	limits[0] = dimSizes[0] - 1;
	for (int dim = 1; dim < dimNum; dim++) {
		limits[dim] += dimSizes[dim] - 1;
	}

	upperBounds[0] = *min_element(limits.begin(), limits.end());
	// go over other dims. from one iteration to the next only limits for current dim and previous have to be changed
	for (int dim = 1; dim < dimNum; dim++) {
		limits[dim] = dimSizes[dim] - 1;
		limits[dim - 1] = D[dim - 1];
		upperBounds[dim] = *min_element(limits.begin(), limits.end());
	}
	// max memory size
	if (dimNum != 0) {
		dimSizePowers[dimNum - 1] = 1;
		for (int i = dimNum - 2; i >= 0; i--) {
			dimSizePowers[i] = dimSizePowers[i + 1] * dimSizes[i + 1];
		}

		maxMemorySize = dimSizePowers[0] * dimSizes[0];
	}
}

void IndexVectorD::operator++(int) {
	assert(indexes.size() != 0);
	if (absoluteIndex == -1) {		// past end. no need to increment
		return;
	}
	int currentDim = dimNum - 1;
	while (currentDim >= 0 and indexes[currentDim] == upperBounds[currentDim]) {
		currentDim--;
	}

	if (currentDim == -1) {		// we are at the end. next vector is past end
		absoluteIndex = -1;
		return;
	}
	// increment current dim index
	indexes[currentDim]++;
	// update bounds of dimensions past it. and set them to lower end of bound
	for (int dim = currentDim + 1; dim < dimNum; dim++) {
		// lower bounds
		vector<int> lowerLimits(dim + 1);
		for (int prevDim = 0; prevDim < dim; prevDim++) {
			lowerLimits[prevDim] = indexes[prevDim] - D[dim];
		}
		lowerLimits[dim] = 0;
		lowerBounds[dim] = *max_element(lowerLimits.begin(), lowerLimits.end());
		indexes[dim] = lowerBounds[dim];

		// upper bounds
		vector<int> upperLimits(D);
		for (int prevDim = 0; prevDim < dim; prevDim++) {
			upperLimits[prevDim] += indexes[prevDim];
		}
		upperLimits[dim] = dimSizes[dim] - 1;
		for (int nextDim = dim + 1; nextDim < dimNum; nextDim++) {
			upperLimits[nextDim] += dimSizes[nextDim] - 1;
		}
		upperBounds[dim] = *min_element(upperLimits.begin(), upperLimits.end());
	}
	UpdateAbsoluteIndex();
}

int IndexVectorD::AbsIndBeta(const vector<int>& K) const {
	int kAbsoluteIndex = absoluteIndex;
	for (vector<int>::const_iterator itvec = K.begin(); itvec != K.end(); itvec++) {
		int currentDimIndex = *itvec;
		assert(indexes[currentDimIndex] > 0);
		kAbsoluteIndex -= dimSizePowers[currentDimIndex];
	}
	return kAbsoluteIndex;
}

vector<int> IndexVectorD::BetaIndexesOnly(const vector<int>& K) const {
	vector<int> result(indexes);
	for (vector<int>::const_iterator itvec = K.begin(); itvec != K.end(); itvec++) {
		int currentDimIndex = *itvec;
		assert(indexes[currentDimIndex] > 0);
		result[currentDimIndex]--;
	}
	return result;
}

IndexVectorD IndexVectorD::BetaClass(const vector<int>& K) const {
	IndexVectorD result = *this;
	for (vector<int>::const_iterator itvec = K.begin(); itvec != K.end(); itvec++) {
		int currentDimIndex = *itvec;
		assert(result.indexes[currentDimIndex] > 0);
		result.indexes[currentDimIndex]--;
		result.absoluteIndex -= dimSizePowers[currentDimIndex];
	}
	return result;
}

IndexVectorD& IndexVectorD::Last(){
	absoluteIndex = maxMemorySize - 1;
	for (int i = 0; i < dimNum; i++) {
		indexes[i] = dimSizes[i] - 1;
	}
	return *this;
}

bool IndexVectorD::IsEnd() const {
	assert(indexes.size() > 0);
	return absoluteIndex == -1;
}
const vector<int>& IndexVectorD::Vector() const {
	return indexes;
}
int IndexVectorD::MaxMemorySize() const {
	return maxMemorySize;
}
void IndexVectorD::UpdateAbsoluteIndex() {
	absoluteIndex = 0;
	for (int i = 0; i < dimNum; i++) {
		absoluteIndex += indexes[i] * dimSizePowers[i];
	}
}
int IndexVectorD::AbsInd() const {
	return absoluteIndex;
}

#include "NChooseK.hpp"
#include <algorithm>
#include <map>
using namespace std;

/* arr[]  ---> Input Array
 data[] ---> Temporary array to store current combination
 start & end ---> Staring and Ending indexes in arr[]
 index  ---> Current index in data[]
 k ---> Size of a combination to be created */
void combinationUtil(int arr[], int data[], int start, int end, int index, int k, vector<vector<int> >& ktuples) {
	// Current combination is complete
	if (index == k) {
		vector<int> currentTuple(k);
		for (int j = 0; j < k; j++)
			currentTuple[j] = data[j];
		ktuples.push_back(currentTuple);
		return;
	}

	// replace index with all possible elements. The condition "end-i+1 >= r-index" makes sure that including one
	// element at index will make a combination with remaining elements at remaining positions
	for (int i = start; i <= end && end - i + 1 >= k - index; i++) {
		data[index] = arr[i];
		combinationUtil(arr, data, i + 1, end, index + 1, k, ktuples);
	}
}

// Create all possible combinations of size k from [0,...,n-1]
vector<vector<int> > Combinations(const int n, const int k) {
	vector<vector<int> > ktuples;
	int arr[n];
	for (int i = 0; i < n; i++) {
		arr[i] = i;
	}

	// A temporary array to store all combination one by one
	int data[k];

	// Create all combination using temporary array 'data[]'
	combinationUtil(arr, data, 0, n - 1, 0, k, ktuples);
	return ktuples;
}

// vector of all possible string pairings

vector<pair<string, string> > StringPairings(const vector<string>& strings) {
	int n = strings.size();
	vector<vector<int> > ktuplesIndexes = Combinations(n, 2);
	int pairNum = ktuplesIndexes.size();
	vector<pair<string, string> > pairs(pairNum);
	for (int i = 0; i < pairNum; i++) {
		int firstStringIndex = ktuplesIndexes[i][0];
		int secondStringIndex = ktuplesIndexes[i][1];
		pairs[i] = make_pair(strings[firstStringIndex], strings[secondStringIndex]);
	}
	return pairs;
}

vector<vector<string> > StringKtuples(const vector<string>& strings, const int k) {
	int n = strings.size();
	vector<vector<int> > ktuplesIndexes = Combinations(n, k);
	int ktuplesNum = ktuplesIndexes.size();
	vector<vector<string>> ktuples(ktuplesNum, vector<string>(k));
	for (int i = 0; i < ktuplesNum; i++) {
		for (int j = 0; j < k; j++) {
			ktuples[i][j] = strings[ktuplesIndexes[i][j]];
		}
	}
	return ktuples;
}

bool SortKTupleByLen(const string& a, const string& b) {
	return a.size() > b.size();
}

bool SortKTuplesBySumLen(const vector<string>& a, const vector<string>& b) {
	int k = a.size();
	unsigned aLenSum = 0, bLenSum = 0;
	for (int i = 0; i < k; i++) {
		aLenSum += a[i].size();
		bLenSum += b[i].size();
	}
	if (aLenSum < bLenSum) {
		return false;
	}
	else if (aLenSum > bLenSum) {
		return true;
	}
	else {	//sum of lengths of strings is the same.
		for (int j = 0; j < k; j++) {
			if (a[j].size() != b[j].size()) {
				return a[j].size() > b[j].size();
			}
		}
		// all lengths the same
		return false;
	}
}

bool SortKTuplesBySumLetterMax(const vector<string>& a, const vector<string>& b) {
	int k = a.size();
	map<char, vector<int> > letterFreqa, letterFreqb;
	letterFreqa['A'] = vector<int>(k);
	letterFreqa['C'] = vector<int>(k);
	letterFreqa['G'] = vector<int>(k);
	letterFreqa['T'] = vector<int>(k);
	letterFreqb['A'] = vector<int>(k);
	letterFreqb['C'] = vector<int>(k);
	letterFreqb['G'] = vector<int>(k);
	letterFreqb['T'] = vector<int>(k);

	for (int i = 0; i < k; i++) {
		for (unsigned j = 0; j < a[i].size(); j++) {
			letterFreqa[a[i][j]][i]++;
		}
		for (unsigned l = 0; l < b[i].size(); l++) {
			letterFreqb[b[i][l]][i]++;
		}
	}
	int maxAa = *max_element(letterFreqa['A'].begin(), letterFreqa['A'].end());
	int maxCa = *max_element(letterFreqa['C'].begin(), letterFreqa['C'].end());
	int maxGa = *max_element(letterFreqa['G'].begin(), letterFreqa['G'].end());
	int maxTa = *max_element(letterFreqa['T'].begin(), letterFreqa['T'].end());
	int sumMaxa = maxAa + maxCa + maxGa + maxTa;

	int maxAb = *max_element(letterFreqb['A'].begin(), letterFreqb['A'].end());
	int maxCb = *max_element(letterFreqb['C'].begin(), letterFreqb['C'].end());
	int maxGb = *max_element(letterFreqb['G'].begin(), letterFreqb['G'].end());
	int maxTb = *max_element(letterFreqb['T'].begin(), letterFreqb['T'].end());
	int sumMaxb = maxAb + maxCb + maxGb + maxTb;
	if (sumMaxa != sumMaxb){
		return sumMaxa > sumMaxb;
	}
	else {
		return SortKTuplesBySumLen(a,b);
	}
}

vector<vector<string> > SortedStringKtuplesBySumLetterMax(const vector<string>& strings, const int k) {
	int n = strings.size();
	vector<vector<int> > ktuplesIndexes = Combinations(n, k);
	int ktuplesNum = ktuplesIndexes.size();
	vector<vector<string>> ktuples(ktuplesNum, vector<string>(k));
	for (int i = 0; i < ktuplesNum; i++) {
		for (int j = 0; j < k; j++) {
			ktuples[i][j] = strings[ktuplesIndexes[i][j]];
		}
		sort(ktuples[i].begin(), ktuples[i].end(), SortKTupleByLen);
	}
	sort(ktuples.begin(), ktuples.end(), SortKTuplesBySumLetterMax);

	return ktuples;
}

vector<vector<string> > SortedStringKtuplesBySumLen(const vector<string>& strings, const int k) {
	int n = strings.size();
	vector<vector<int> > ktuplesIndexes = Combinations(n, k);
	int ktuplesNum = ktuplesIndexes.size();
	vector<vector<string>> ktuples(ktuplesNum, vector<string>(k));
	for (int i = 0; i < ktuplesNum; i++) {
		for (int j = 0; j < k; j++) {
			ktuples[i][j] = strings[ktuplesIndexes[i][j]];
		}
		sort(ktuples[i].begin(), ktuples[i].end(), SortKTupleByLen);
	}
	sort(ktuples.begin(), ktuples.end(), SortKTuplesBySumLen);

	return ktuples;
}

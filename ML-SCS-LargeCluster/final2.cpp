
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <map>
#include <set>

#include <string>
#include <iterator>
#include <time.h>
#include <math.h>
#include <random>
#include <chrono>
#include <cstdint>
#include <climits>
#include <fstream>

using namespace std;
#define DES_LEN 200
#define HISTO_LEN 1000

string _correct_binary_indel(int n, int m, int a, string y);
string _correct_q_ary_indel(int n, int m, int a, int b, int q, string y);
string convert_y_to_alpha(string y);
string genRandomSeq(int len);
vector<int> compute_indices(char y1, string x);
int count_embeddings(string a, string b);
string maximum_likelihood_old(set<string> &candidates, const vector<string> &cluster);
string maximum_likelihood(set<string> &candidates, const vector<string> &cluster);
void createNoisyCluster(vector<string> &cluster, string original, int cluster_size, float del_prob);
string reconstruction_new(vector<string> &cluster, int des_len);
string reconstruction_new_alg4(vector<string> &cluster, int des_len);
void reomoveByLength(set<string> &candidates, int len);
void removeNonSuperseq_withLen(set<string> &candidates, vector<string> &cluster, int len);
void removeNonSuperseq(set<string> &candidates, vector<string> &cluster, int len);
int compute_syndrome_binary(int m, int a, string y);
string convert_y_to_alpha(string y);
int compute_syndrome_q_ary(int m, int a, int b, int q, string y, int &sec);
bool is_codeword(string &y, int q);
string decode_VT(string y, int q, int n, int m , int a, int b);

int design_len=DES_LEN;
int total_dis_tests = 0;
float prob = 0.14;
int clus_size = 6;
int num_of_test = 100000;
int count=0;
int q=4;
int scs_size=0;
int scs_after_filter_size=0;
int scs_most_likely_size=0;
int most_likelihood_count=0;
int pairs_histo[HISTO_LEN]={0};
int three_histo[HISTO_LEN]={0};
int four_histo[HISTO_LEN]={0};
int SR_histo[10]={0};
int SR_histo_BMA[10]={0};

set<string> ML;
map<int, int> scs_sizes;
map<int, int> scs_sizes_count;
map<int, int> mst_like;
map<int, int> mst_like_count;
//int table3[4][4]={{0}};
//int table4[4][4]={{0}};
//int table_l_3={0,1,3,4}=;
//int table_r_3=;
double error_rate_by_length[DES_LEN+50]={0};
double count_by_length[DES_LEN+50]={0};
map < string, int > table_1_rev;
map < string, int > table_2_rev;
vector<int> systematic_positions_step_1;
int lucky = 0;
int two = 0;
int three = 0;
int four = 0;
int nothing = 0;
ofstream myfile;



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

unsigned int edit_distance(const std::string& s1, const std::string& s2)
{
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<std::vector<unsigned int>> d(len1 + 1, std::vector<unsigned int>(len2 + 1));

    d[0][0] = 0;
    for(unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
    for(unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

    for(unsigned int i = 1; i <= len1; ++i)
        for(unsigned int j = 1; j <= len2; ++j)
                      // note that std::min({arg1, arg2, arg3}) works only in C++11,
                      // for C++98 use std::min(std::min(arg1, arg2), arg3)
                      d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
    return d[len1][len2];
}

/* only ins - del */
unsigned int levenstein_dis(const std::string& s1, const std::string& s2)
{
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<std::vector<unsigned int>> d(len1 + 1, std::vector<unsigned int>(len2 + 1));

    d[0][0] = 0;
    for(unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
    for(unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

    for(unsigned int i = 1; i <= len1; ++i)
        for(unsigned int j = 1; j <= len2; ++j)
                      // note that std::min({arg1, arg2, arg3}) works only in C++11,
                      // for C++98 use std::min(std::min(arg1, arg2), arg3)
                      d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1000) });
    return d[len1][len2];
}


void reverseStr(string& str)
{
    int n = str.length();
  
    // Swap character starting from two
    // corners
    for (int i = 0; i < n / 2; i++)
        swap(str[i], str[n - i - 1]);
}


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

map<char, vector<int> > LettersMap(const vector<string>& strings, const vector<int>& index) {
	map<char, vector<int> > lettersMap;
	for (int i = 0; i < (int) strings.size(); i++) {
		if (index[i] == 0) {
			lettersMap['$'].push_back(i);
		}
		else {
			lettersMap[strings[i][index[i] - 1]].push_back(i);
		}
	}
	return lettersMap;
}

vector<string> BacktrackSCSFast(const vector<string>& strings, IndexVectorD index, const vector<int>& L) {
	int n = strings.size();
	map<char, vector<int> > lettersMap = LettersMap(strings, index.Vector());
	if ((int) lettersMap['$'].size() >= n - 1) { // at most 1 non zero index.
		string currentStr;
		// find non zero dimension
		for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
			if (it->first != '$') {
				int nonZeroDimension = it->second[0];
				vector<int> indexes = index.Vector();
				currentStr = strings[nonZeroDimension].substr(0, indexes[nonZeroDimension]);
				break;
			}
		}
		return vector<string>( { currentStr });
	}
	else { // at least 2 non zero indexes.
		lettersMap.erase('$');
		//int minLen = INT_MAX, currentValue;
		int minValue = L[index.AbsInd()] - 1, currentValue;
		vector<char> lettersToAdd;
		vector<IndexVectorD> minIndexes;
		for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
			vector<int> K = it->second;
			currentValue = L[index.AbsIndBeta(K)];
			if (currentValue == minValue) {
				lettersToAdd.push_back(it->first);
				minIndexes.push_back(index.BetaClass(K));
			}
		}
		vector<string> result;
		for (unsigned int i = 0; i < minIndexes.size(); i++) {
			vector<string> currentResult = BacktrackSCSFast(strings, minIndexes[i], L);
			for (vector<string>::iterator it = currentResult.begin(); it != currentResult.end(); it++) {
				it->push_back(lettersToAdd[i]);
			}
			result.insert(result.end(), currentResult.begin(), currentResult.end());
		}
		return result;
	}
}

// function to calculate the len of SCS
int SCSNLenFast(const vector<string>& strings, vector<int>& L, const IndexVectorD& firstIndex) {
	int n = strings.size();
	L[0] = 0;
	for (IndexVectorD currentIndex = firstIndex; not currentIndex.IsEnd(); currentIndex++) {
		map<char, vector<int> > lettersMap = LettersMap(strings, currentIndex.Vector());
		if ((int) lettersMap['$'].size() >= n - 1) { // at most 1 non zero index.
			for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
				if (it->first != '$') {
					int nonZeroDimension = it->second[0];
					vector<int> indexes = currentIndex.Vector();
					L[currentIndex.AbsInd()] = indexes[nonZeroDimension];
					break;
				}
			}
		}
		else { // at least 2 non zero indexes.
			lettersMap.erase('$');
			int minLen = INT_MAX, currentValue;
			for (map<char, vector<int> >::iterator it = lettersMap.begin(); it != lettersMap.end(); it++) {
				vector<int> K = it->second;
				currentValue = L[currentIndex.AbsIndBeta(K)];
				if (currentValue < minLen) {
					minLen = currentValue;
				}
			}

			// if minLen==INT_MAX at this stage there is no legit path to this vector, so leave it at INT_MAX
			if (minLen < INT_MAX) {
				L[currentIndex.AbsInd()] = minLen + 1;
			}
		}
	}
	assert(L[firstIndex.MaxMemorySize() - 1] != INT_MAX); // There is always a path
	return L[firstIndex.MaxMemorySize() - 1];
}

// Omer: the function for all SCSs.
vector<string> SCSNFast(const vector<string>& strings, const int originalLen) {
	int n = strings.size();
	vector<int> dimSizes(n);
	vector<int> D(n);
	for (int i = 0; i < n; i++) {
		dimSizes[i] = strings[i].size() + 1;
		D[i] = originalLen - strings[i].size();
	}
	IndexVectorD indexVector(dimSizes, D);
	vector<int> L(indexVector.MaxMemorySize(), INT_MAX);
	SCSNLenFast(strings, L, indexVector);
	return BacktrackSCSFast(strings, indexVector.Last(), L);
}
// Omer: the function for SCS length only
int SCSNLenFast(const vector<string>& strings, const int originalLen) {
	int n = strings.size();
	vector<int> dimSizes(n);
	vector<int> D(n);
	for (int i = 0; i < n; i++) {
		dimSizes[i] = strings[i].size() + 1;
		D[i] = originalLen - strings[i].size();
	}
	IndexVectorD indexVector(dimSizes, D);
	vector<int> L(indexVector.MaxMemorySize(), INT_MAX);
	return SCSNLenFast(strings, L, indexVector);
}


/* create a random string of a given length on alphabet sized q */
string genRandomSeq(int len, int q){
    //random_device random_device; // create object for seeding
    //mt19937 engine{random_device()}; // create engine and seed it
    unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
    mt19937 engine(sd);
    uniform_int_distribution<> dist(0, q-1); // create distribution for integers with [1; 9] range
    string res = "";
    for(int i=0; i<len; i++){
        int iSecret = dist(engine); // generate pseudo random int
        //char c=intToDNA(iSecret);
        char c = iSecret + '0';
        //cout << iSecret;
        //cout << c;
        res.push_back(c);
    }
    return res;
}

// given a and b count num of embeddings of b in a
// Iterative DP function to find the number of times
// the second string occurs in the first string,
// whether continuous or discontinuous
int count_embeddings(string a, string b)
{
    int m = a.length();
    int n = b.length();
    //int lookup[MAX][MAX];

    // Create a table to store results of sub-problems
   // int lookup[m + 1][n + 1] = { { 0 } };
    int** lookup = new int*[m+1];
    for(int i = 0; i < m+1; ++i){
        lookup[i] = new int[n+1];
        for(int j=0; j < n+1; j++){
            lookup[i][j]=0;
        }
    }
    // If first string is empty
    for (int i = 0; i <= n; ++i)
        lookup[0][i] = 0;
  
    // If second string is empty
    for (int i = 0; i <= m; ++i)
        lookup[i][0] = 1;
  
    // Fill lookup[][] in bottom up manner
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            // If last characters are same, we have two
            // options -
            // 1. consider last characters of both strings
            //    in solution
            // 2. ignore last character of first string
            if (a[i - 1] == b[j - 1])
                lookup[i][j] = lookup[i - 1][j - 1] +
                               lookup[i - 1][j];
                  
            else
                // If last character are different, ignore
                // last character of first string
                lookup[i][j] = lookup[i - 1][j];
        }
    }
    int res = lookup[m][n];

    for(int i = 0; i < m+1; ++i)
        delete lookup[i];
    delete[] lookup;
    return res;
}

// might be better for me
/* gets a list of candidates - SCS with the lengths, and return the maximum likelihood to be the original string*/
string maximum_likelihood_old(set<string> &candidates, const vector<string> &cluster){
    double max_prob = 0.0;
    string maxStr = "";
    if(candidates.empty()){
        return "ERROR";
    }
    for(string candid: candidates){
        if(candid=="")
            continue;
        double prob = 1;
        auto it = cluster.begin();
        bool ok=false;
        while( it != cluster.end()){
            ok=true;
            prob*=count_embeddings(candid, *it);
            it++;
        }
        if (prob >= max_prob && ok==true){
            max_prob=prob;
            maxStr=candid;
        }
    }
    return maxStr;

}


/* gets a list of candidates - SCS with the lengths, and return the maximum likelihood to be the original string*/
string maximum_likelihood(set<string> &candidates, const vector<string> &cluster){
    double max_prob = 0.0;
    string maxStr = *(candidates.begin());
    if(candidates.empty()){
        return "ERROR";
    }
    int candid_len;
    for(string candid: candidates){
        if(candid=="")
            continue;
        candid_len=candid.length();
        double prob = 1;
        auto it = cluster.begin();
        bool ok=false;
        while( it != cluster.end()){
            ok=true;
            prob*=count_embeddings(candid, *it);
            it++;
        }
        if (prob >= max_prob && ok==true){
            max_prob=prob;
            maxStr=candid;
        }
    }
    int counter=0;
    // count number of most likely hood
    for(string candid: candidates){
        if(candid=="")
            continue;
        double prob = 1;
        auto it = cluster.begin();
        bool ok=false;
        while( it != cluster.end()){
            ok=true;
            prob*=count_embeddings(candid, *it);
            it++;
        }
        if (prob == max_prob){
            ML.insert(candid);
            counter++;
        }
    }
    if(mst_like.find(candid_len)==mst_like.end()){
        mst_like[candid_len]=candidates.size();
        mst_like_count[candid_len]=1;
    }
    else{
        mst_like[candid_len]+=candidates.size();
        mst_like_count[candid_len]+=1;
    }
    
    return maxStr;
}

char generateRandomDNA(float prob){
    if(prob <=0.25){
        return '0';
    }
    if(prob <= 0.5){
        return '1';
    }
    if(prob <= 0.75){
        return '2';
    }
    if(prob <= 1.0){
        return '3';
    }
    return 'X';
}

// gets a string and creates noisy cluster //
void createNoisyCluster(vector<string> &cluster, string original, int cluster_size, float del_prob){
    /* initialize random seed: */
    //srand (time(NULL));
    //random_device random_device; // create object for seeding
    //mt19937 engine{random_device()}; // create engine and seed it
    //uniform_int_distribution<> dist(1,100); // create distribution for integers with [1; 9] range
    unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
    mt19937 engine(sd);
    std::uniform_real_distribution<> dist(0.0, 1.0);
    double iSecret;
    float prob;
    for (int i =0; i < cluster_size ; i++){
        string new_string  ="";
        int count=0;
        for(int j=0 ; j<original.length(); j++ ){
            iSecret = dist(engine); // generate pseudo random number
            prob=iSecret;
            if(prob<=1-del_prob){
                new_string.push_back(original[j]);
            }else{
                count++;
                continue;
            }
        }
        cluster.push_back(new_string);
        new_string = "";
    }
}

// this function remove candidate out of set of candidates  by their length
void reomoveByLength(set<string> &candidates, int len){
    auto it = candidates.begin();
    while (it != candidates.end()){
        string candid = *it;
        
        if (candid.length() != len){
            it = candidates.erase(it);
        }
        else{
            it++;
        }
    }
    
    it = candidates.begin();
    while (it != candidates.end()){
        string candid = *it;
        
        if (candid.length() != len){
            cout << "prob" <<endl;
        }
        it++;
    }

}

// Returns true if str1[] is a subsequence of str2[]. m is
// length of str1 and n is length of str2
bool isSubSequence(string str1, string str2, int m, int n)
{
    int j = 0; // For index of str1 (or subsequence
    
    // Traverse str2 and str1, and compare current character
    // of str2 with first unmatched char of str1, if matched
    // then move ahead in str1
    for (int i=0; i<n&&j<m; i++)
        if (str1[j] == str2[i])
            j++;
    
    // If all characters of str1 were found in str2
    return (j==m);
}

// this function remove candidate out of set of candidates  by their length
void removeNonSuperseq_withLen(set<string> &candidates, vector<string> &cluster, int len){
    auto it = candidates.begin();
    bool del = false;
    while (it != candidates.end()){
        del = false;
        string candid = *it;
        for(string read: cluster){
            if(read.length() > candid.length())
                continue;
            if (!isSubSequence(read, candid, read.length(), candid.length())){
                it = candidates.erase(it);
                del = true;
                break;
            }
        }
        if(!del){
            it++;
        }
    }

}



// this function remove candidate out of set of candidates  by their length
void removeNonSuperseq(set<string> &candidates, vector<string> &cluster, int len){
    auto it = candidates.begin();
    bool del = false;
    while (it != candidates.end()){
        del = false;
        string candid = *it;
        for(string read: cluster){
            if (!isSubSequence(read, candid, read.length(), candid.length())){
                it = candidates.erase(it);
                del = true;
                break;
            }
        }
        if(!del){
            it++;
        }
    }

}


string recon_alex_omer(const vector<string> &cluster, const int des_len, int &origSize){
    //cout << "bef len" << endl;

    int scs_len = SCSNLenFast(cluster, des_len);
    //cout << "after len" << endl;
    origSize=scs_len;
    if(scs_len==des_len){
        vector<string> SCS_set = SCSNFast(cluster, des_len);
        //cout << "after SCS in len" << endl;
        scs_size=SCS_set.size();
        if(scs_sizes.find(des_len)==scs_sizes.end()){
            scs_sizes[des_len]=0;
            scs_sizes_count[des_len]=0;
            scs_sizes_count[des_len]++;
            scs_sizes[des_len]+=SCS_set.size();
        }
        else{
            scs_sizes[des_len]+=SCS_set.size();
            scs_sizes_count[des_len]++;
        }
        set<string> scs_set;
        for (string &i: SCS_set) {
            scs_set.insert(i);
        }
        return maximum_likelihood(scs_set, cluster);
    }
    else{
        vector<string> SCS_set = SCSNFast(cluster, des_len);
        //cout << "after SCS not len" << endl;
        scs_size=SCS_set.size();
        set<string> scs_set;
        for (string &i: SCS_set) {
            scs_set.insert(i);
        }
        return maximum_likelihood(scs_set, cluster);
        
        //cout << "not the length " << scs_len << endl;
        //return "";
    }
}

int find_index_of_short(const vector<string> &cluster){
    int i=0;
    int min_len = cluster[0].length();
    int min_ind=0;
    for(string str: cluster){
        if(str.length()<=min_len){
            min_ind=i;
            min_len=str.length();
        }
        i++;
    }
    return min_ind;
}

int find_index_of_long(const vector<string> &cluster){
    int i=0;
    int max_len = cluster[0].length();
    int max_ind=0;
    for(string str: cluster){
        if(str.length()>=max_len){
            max_len=i;
            
        }
        i++;
    }
    return max_len;
}

/* check if str is super sequence of all sequence in a given cluster*/
bool isSupSeq(string str, const vector<string> &cluster){
    for(string cp:cluster){
        if(!isSubSequence(cp, str, cp.length(), str.length()))
            return false;
    }
    return true;
}

string recon_2(const vector<string> &cluster, const int des_len, int &origSize){
    for(string cp1: cluster){
        for(string cp2: cluster){
            vector<string> cluster_pairs;
            cluster_pairs.push_back(cp1);
            cluster_pairs.push_back(cp2);
            int scs_len = SCSNLenFast(cluster_pairs, des_len);
            if(scs_len== des_len){
                vector<string> SCS_set = SCSNFast(cluster_pairs, des_len);
                scs_size=SCS_set.size();
                if(scs_sizes.find(des_len)==scs_sizes.end()){
                    scs_sizes[des_len]=0;
                    scs_sizes_count[des_len]=0;
                    scs_sizes_count[des_len]++;
                    scs_sizes[des_len]+=SCS_set.size();
                }
                else{
                    scs_sizes[des_len]+=SCS_set.size();
                    scs_sizes_count[des_len]++;
                }
                set<string> scs_set;
                for (string &i: SCS_set) {
                    if(isSupSeq(i, cluster))
                        scs_set.insert(i);
                }
                if(scs_set.size()>0)
                    return maximum_likelihood(scs_set, cluster);
            }
        }
    }
    return "";
        
}

string recon_3(const vector<string> &cluster, const int des_len, int &origSize){
    for(string cp1: cluster){
        for(string cp2: cluster){
            for(string cp3: cluster){
                vector<string> cluster_three;
                cluster_three.push_back(cp1);
                cluster_three.push_back(cp2);
                cluster_three.push_back(cp3);
                int scs_len = SCSNLenFast(cluster_three, des_len);
                if(scs_len== des_len){
                    vector<string> SCS_set = SCSNFast(cluster_three, des_len);
                    scs_size=SCS_set.size();
                    if(scs_sizes.find(des_len)==scs_sizes.end()){
                        scs_sizes[des_len]=0;
                        scs_sizes_count[des_len]=0;
                        scs_sizes_count[des_len]++;
                        scs_sizes[des_len]+=SCS_set.size();
                    }
                    else{
                        scs_sizes[des_len]+=SCS_set.size();
                        scs_sizes_count[des_len]++;
                    }
                    set<string> scs_set;
                    for (string &i: SCS_set) {
                        if(isSupSeq(i, cluster))
                            scs_set.insert(i);
                    }
                    if(scs_set.size()>0)
                        return maximum_likelihood(scs_set, cluster);
                }
            }
        }
    }
    return "";
}

string recon_4(const vector<string> &cluster, const int des_len, int &origSize){
    for(string cp1: cluster){
        for(string cp2: cluster){
            for(string cp3: cluster){
                for(string cp4: cluster){
                    vector<string> cluster_four;
                    cluster_four.push_back(cp1);
                    cluster_four.push_back(cp2);
                    cluster_four.push_back(cp3);
                    cluster_four.push_back(cp4);

                    int scs_len = SCSNLenFast(cluster_four, des_len);
                    if(scs_len== des_len){
                        vector<string> SCS_set = SCSNFast(cluster_four, des_len);
                        scs_size=SCS_set.size();
                        if(scs_sizes.find(des_len)==scs_sizes.end()){
                            scs_sizes[des_len]=0;
                            scs_sizes_count[des_len]=0;
                            scs_sizes_count[des_len]++;
                            scs_sizes[des_len]+=SCS_set.size();
                        }
                        else{
                            scs_sizes[des_len]+=SCS_set.size();
                            scs_sizes_count[des_len]++;
                        }
                        set<string> scs_set;
                        for (string &i: SCS_set) {
                            if(isSupSeq(i, cluster))
                                scs_set.insert(i);
                        }
                        if(scs_set.size()>0)
                            return maximum_likelihood(scs_set, cluster);
                    }
                }
            }
        }
    }
    return "";
}

int hamming_dis(string s1, string s2){
    int dis=0;
    int i=0;
    for(i=0; i<s1.length() && i<s2.length(); i++){
        if(s1[i]!=s2[i]){
            dis++;
        }
    }
    while(i<s1.length()){
        dis++;
        i++;
    }
    while(i<s2.length()){
        dis++;
        i++;
    }
    return dis;
    
}

vector<string> findMaxHamming(const vector<string> &cluster){
    int max_dis=0;
    string max_cp1, max_cp2;
    vector<string> max_hamm_pair;
    for(string cp1: cluster){
        for(string cp2: cluster){
            int h_dis=hamming_dis(cp1, cp2);
            if(h_dis>=max_dis){
                max_dis=h_dis;
                max_cp1=cp1;
                max_cp2=cp2;
            }
        }
    }
    max_hamm_pair.push_back(max_cp1);
    max_hamm_pair.push_back(max_cp2);
    return max_hamm_pair;
}

vector<string> findMaxLev(const vector<string> &cluster){
    int max_dis=0;
    string max_cp1, max_cp2;
    vector<string> max_hamm_pair;
    for(string cp1: cluster){
        for(string cp2: cluster){
            int h_dis=levenstein_dis(cp1, cp2);
            if(h_dis>=max_dis){
                max_dis=h_dis;
                max_cp1=cp1;
                max_cp2=cp2;
            }
        }
    }
    max_hamm_pair.push_back(max_cp1);
    max_hamm_pair.push_back(max_cp2);
    return max_hamm_pair;
}

bool compareLen(const std::string& a, const std::string& b)
{
    if(a.size()==b.size()){
        return a>b;
    }
    return (a.size() > b.size());
}

bool compareLen_sets(const vector <string> &a, const vector <string> &b)
{
    int sum_a=0, sum_b=0;
    for(string s: a){
        sum_a+=s.length();
    }
    for(string s: b){
        sum_b+=s.length();
    }
    if(sum_a < sum_b){
        return false;
    }
    if(sum_a>sum_b){
        return true;
    }
    int v_size=a.size();
    //cout << v_size << endl;
    if(sum_a==sum_b){
        if(a[0].length()>b[0].length()){
            return true;
        }
        if(a[1].length()>b[1].length()){
            return true;
        }
        return false;

    }
    return false;

}

string minED(const vector<string> &cluster, int des_len){
    vector<string> clus=cluster;
    sort(clus.begin(), clus.end(), compareLen);
    vector<string> maxFour;
    maxFour.push_back(clus[0]);
    maxFour.push_back(clus[1]);
    maxFour.push_back(clus[2]);
    maxFour.push_back(clus[3]);
    vector<string> SCS_set=SCSNFast(maxFour, des_len);
    int min_dis=100;
    int curr_dis=0;
    string min_SCS=SCS_set[0];
    for(string scs: SCS_set){
        for(string cp: cluster){
            curr_dis+=levenstein_dis(scs, cp);
            
        }
        if(curr_dis<min_dis){
            min_dis=curr_dis;
            min_SCS=scs;
        }
        curr_dis=0;
    }
    return min_SCS;
}

string min_ed_candidate(const vector<string> &cluster, set<string> &candidates, int des_len){
    int curr_dis=0;
    string min_SCS=*(candidates.begin());
    int min_dis=0;
    for(string cp: cluster){
        min_dis+=levenstein_dis(*(candidates.begin()), cp);
    }
    for(string scs: candidates){
        for(string cp: cluster){
            curr_dis+=levenstein_dis(scs, cp);
        }
        if(curr_dis<min_dis){
            min_dis=curr_dis;
            min_SCS=scs;
        }
        curr_dis=0;
    }
    return min_SCS;
}

vector<string> findTwoLong(const vector<string> &cluster){
    int len1=0;
    string longest1, longest2;
    vector<string> twolngest;

    for(string cp1: cluster){
        if(cp1.length()>=len1){
            longest1=cp1;
            len1=cp1.length();
        }
    }
    int len2=0;
    for(string cp1: cluster){
        if(cp1.length()>=len2 && cp1!=longest1){
            longest2=cp1;
            len2=cp1.length();
        }
    }
    twolngest.push_back(longest1);
    twolngest.push_back(longest2);
    return twolngest;
}



string ML_SuperSeqSCS(vector<string> f_cluster, const vector<string> &cluster, int &emd, int des_len){
    emd = 0;
    vector<string> SCS_set = SCSNFast(f_cluster, des_len);
    set<string> scs_set;
    int nonSuperSeq=0;
    int minNonSuperSeq=0;
    for (string &i: SCS_set) {
        nonSuperSeq = 0;
        for(string cp: cluster){
            if(!isSubSequence(cp, i, cp.length(), i.length())){
                nonSuperSeq++;
            }
        }
        if(nonSuperSeq==0){
            scs_set.insert(i);
        }
    }
    if(scs_set.size()>0){
        //lucky++;
        return maximum_likelihood(scs_set, cluster);
    }
    return "";
}

vector < vector < string> > create_sets(const vector<string> &cluster, int set_size){
    int set_counter=0;
    vector < vector < string > > vec;
    
    if (set_size == 2){
        vector <string> v;
        for(string cp1: cluster){
            for(string cp2: cluster){
                v.push_back(cp1);
                v.push_back(cp2);
                sort( v.begin(), v.end() , compareLen);
                v.erase( unique( v.begin(), v.end() ), v.end() );
                if(v.size()!=set_size){
                    v.clear();
                }
                else{
                    vec.push_back(v);
                    v.clear();
                }
            }
            
        }
        sort(vec.begin(), vec.end() );
        vec.erase(unique( vec.begin(), vec.end() ), vec.end() );
    }
    
    if (set_size == 3){
        vector <string> v;
        for(string cp1: cluster){
            for(string cp2: cluster){
                for(string cp3: cluster){
                    v.push_back(cp1);
                    v.push_back(cp2);
                    v.push_back(cp3);
                    sort( v.begin(), v.end() , compareLen);
                    v.erase( unique( v.begin(), v.end() ), v.end() );
                    if(v.size()!=set_size){
                        v.clear();
                    }
                    else{
                        vec.push_back(v);
                        v.clear();
                    }
                }
            }
        }
        sort(vec.begin(), vec.end() );

        vec.erase(unique( vec.begin(), vec.end() ), vec.end() );
    }
    
    if (set_size == 4){
        vector <string> v;
        for(string cp1: cluster){
            for(string cp2: cluster){
                for(string cp3: cluster){
                    for(string cp4: cluster){
                        v.push_back(cp1);
                        v.push_back(cp2);
                        v.push_back(cp3);
                        v.push_back(cp4);
                        sort( v.begin(), v.end() , compareLen);
                        v.erase( unique( v.begin(), v.end() ), v.end() );
                        if(v.size()!=set_size){
                            v.clear();
                        }
                        else{
                            vec.push_back(v);
                            v.clear();
                        }

                    }
                }
            }
        }
        sort(vec.begin(), vec.end() );

        vec.erase(unique( vec.begin(), vec.end() ), vec.end() );
    }

    return vec;
}

string _recon(const vector<string> &cluster, const int des_len, int &origSize){
    vector<string> max_hamm;
    max_hamm=findTwoLong(cluster);
    int scs_len = SCSNLenFast(max_hamm, des_len);
    origSize=scs_len;
    int emd=cluster.size();
    string guess="";
    if(scs_len==des_len){
        vector<string> SCS_set = SCSNFast(max_hamm, des_len);
        string res=ML_SuperSeqSCS(max_hamm, cluster, emd, des_len);
        if(res!=""){
            guess=res;
        }
        cout << guess << " " << emd;
        
    }
    string str = recon_2(cluster, des_len, origSize);
    if(str!=""){
=        two++;
        return str;
    }
    str = recon_3(cluster, des_len, origSize);
    if(str!=""){
        three++;
        return str;
    }
    str = recon_4(cluster, des_len, origSize);
    if(str!=""){
=        four++;
        return str;
    }
    nothing++;

    return minED(cluster, des_len);
}


string _recon_final(const vector<string> &cluster, const int des_len, int &origSize){
    int emd=cluster.size();
    string guess="";
    string res="";
    // PAIRS
    int max_len = 0 ;
    vector <vector <string> >  max_set;
    vector <vector <string> > vec = create_sets(cluster, 2);
    int count=0;
    for(vector<string> v: vec){
        count++;
        int scs_len = SCSNLenFast(v, des_len);
        if(scs_len==des_len){
            //cout << "*/";
            vector<string> SCS_set = SCSNFast(v, des_len);
            int original_emd=emd;
            res=ML_SuperSeqSCS(v, cluster, emd, des_len);
            if(res!=""){
                two++;
                if(count<HISTO_LEN)
                    pairs_histo[count]++;
                else
                    pairs_histo[HISTO_LEN-1]++;
                return res;
            }
        }
        else{
            if(scs_len > max_len){
                max_len=scs_len;
                max_set.clear();
                max_set.push_back(v);
            }
            if(scs_len == max_len){
                max_set.push_back(v);
            }
        }
    }


    // THREE
    vec = create_sets(cluster, 3);
    //cout << "three" << endl;
    //cout << vec.size() << endl;
    count=0;
    for(vector<string> v: vec){
        count++;
        int scs_len = SCSNLenFast(v, des_len);
        if(scs_len==des_len){
            //cout << "*/";
            vector<string> SCS_set = SCSNFast(v, des_len);
            int original_emd=emd;
            res=ML_SuperSeqSCS(v, cluster, emd, des_len);
            if(res!=""){
                three++;
                if(count<HISTO_LEN)
                    three_histo[count]++;
                else
                    three_histo[HISTO_LEN-1]++;
                return res;
            }
        }
        else{
            if(scs_len > max_len){
                max_len=scs_len;
                max_set.clear();
                max_set.push_back(v);
            }
            if(scs_len == max_len){
                max_set.push_back(v);
            }
        }
    }
    
    // FOUR
    //cout << "four" << endl;
    vec = create_sets(cluster, 4);
    //cout << vec.size() << endl;
    count=0;

    for(vector<string> v: vec){
        count++;
        int scs_len = SCSNLenFast(v, des_len);
        if(scs_len==des_len){
            vector<string> SCS_set = SCSNFast(v, des_len);
            int original_emd=emd;
            res=ML_SuperSeqSCS(v, cluster, emd, des_len);
            if(res!=""){
                four++;
                if(count<HISTO_LEN)
                    four_histo[count]++;
                else
                    four_histo[HISTO_LEN-1]++;
                return res;
            }
        }
        else{
            if(scs_len > max_len){
                max_len=scs_len;
                max_set.clear();
                max_set.push_back(v);
            }
            if(scs_len == max_len){
                max_set.push_back(v);
            }
        }
    }
    nothing++;
    if(res==""){
        set<string> total_candids;
        for(vector<string> ms:max_set){
            vector<string> SCS_set = SCSNFast(ms, des_len);
            for(string i: SCS_set){
                total_candids.insert(i);
            }
        }
        set <string> candidates;
        for(string i: total_candids){
            int nonSuperSeq = 0;
            for(string cp: cluster){
                if(!isSubSequence(cp, i, cp.length(), i.length())){
                    nonSuperSeq++;
                }
            }
            if(nonSuperSeq==0){
                candidates.insert(i);
            }
        }
        if(candidates.size()>0){
            return maximum_likelihood(candidates, cluster);
        }
        else{
            return min_ed_candidate(cluster, total_candids, des_len);
        }
    }
    return guess;

}



/// BMA algo

void padCluster(vector<string> &cluster, int des_len){
    for(string str: cluster){
        int mis = des_len - str.length();
        while(mis>0){
            str.push_back('-');
            mis--;
        }
    }
}


void createPointers(vector<int> &pointers, int cluster_size){
    for(int i=0; i<cluster_size; i++){
        pointers.push_back(0);
    }
    
}

char BMAMajority(vector<string> &cluster, vector<int> &pointers){
    //cout << "set  size is " << sequences.size() << endl;
    int majority[5] = {0};
    int max = majority[0];
    int max_i =0;
    auto  it  =cluster.begin();
    int current_seq=0;
    while(it!=cluster.end()){
        string  seq = *it;
        if(seq.empty()){
            current_seq++;
            it++;
            continue;
        }
        if(seq[pointers[current_seq]]==0){
            current_seq++;
            it++;
            continue;
        }
        if((seq[pointers[current_seq]])-'0' <= -1 || (seq[pointers[current_seq]])-'0' >=4){
            cerr << "darlin " << " " << pointers[current_seq] << " " << seq[pointers[current_seq]] << endl;
            //throw MyException("INVALID DNA Char ");
            //cout << "ERROERERE  " << endl;
        }
        else{
            majority[(seq[pointers[current_seq]])-'0']++;
        }
        current_seq++;
        it++;
    }


    for(int i=0; i<4; i++){
        if(majority[i] > max) {
            max_i= i;
            max = majority[i];
        }
    }
    //cout  << intToDNA(max_i) << endl;
    return (max_i)+'0';
}

void updateLetter(vector<string> &cluster, vector<int> &pointers, string &res, char c){
    res.push_back(c);
    int current=0;
    for(string seq:cluster){
        if(seq[pointers[current]]==c){
            pointers[current]++;
        }
        current++;
    }
    return;
}


string BMA(vector<string> &cluster, int des_len, int &origSize){

    string res = "";
    vector<int> pointers;
    int cluster_size = cluster.size();

    padCluster(cluster, des_len);

    createPointers(pointers, cluster_size);
    for(int i=0; i<des_len-1; i++){
        int c=BMAMajority(cluster, pointers);
        //cout << c;
        updateLetter(cluster, pointers, res, c);
        //cout <<endl;
    }
    //cout << res << endl;
    return res;
}

int NCR(int n, int r)
{
    if (r == 0) return 1;

    /*
     Extra computation saving for large R,
     using property:
     N choose R = N choose (N-R)
    */
    if (r > n / 2) return NCR(n, n - r);

    long res = 1;

    for (int k = 1; k <= r; ++k)
    {
        res *= n - k + 1;
        res /= k;
    }

    return res;
}

// main function
int main()
{
    double succes_rate=0.0;
    try{
        int counterOfGoodLen=0;
        string original;
    vector<string> cluster;
    double succes_rate_BMA=0.0;
    int total_dis_tests_BMA = 0;

    for(int j=0; j<1; j++){

        cout << "Cluster size :  " << clus_size << " Deletion Prob : " << prob << " Number of tests " << num_of_test << " str len : "<< design_len <<" lower bound : " << 2.2*prob*prob << " q: " << q  << endl << endl;
        int origSize=0;
        for(int i =0; i <num_of_test; i++){
            /*do{
                original=genRandomSeq(design_len, q);
            }
            while(!is_codeword(original, q));*/
            //clus_size = 6;
            original=genRandomSeq(design_len, q);
            createNoisyCluster(cluster, original, clus_size, prob);
            sort(cluster.begin(), cluster.end(), compareLen);
            string recon = BMA(cluster, design_len, origSize);
            int  dis = int(levenstein_dis(recon, original));
            for (int i=dis; i<10; i++){
                SR_histo_BMA[i]++;
                //cout << "upd " << i << endl;
            }
            total_dis_tests_BMA += dis;
            succes_rate_BMA+= dis > 0 ? 0 : 1;
            scs_size=0;
            ML.clear();
            recon = _recon_final(cluster, design_len, origSize);

            dis = int(levenstein_dis(recon, original));
            for (int i=dis; i<10; i++){
                SR_histo[i]++;
                //cout << "upd " << i << endl;
            }
            total_dis_tests += dis;
            succes_rate+= dis > 0 ? 0 : 1;
            scs_size=0;
            ML.clear();
            if (i%100 == 0){
                cout << i << endl;
                cout << "******* SCS ALG ********"<< endl;

                double res = (double)(total_dis_tests) / (i+1);
                res = res / original.length();
                cout << "Error rate: " << res << endl;
                cout << "Succes rate: " << succes_rate/(i+1) << endl;
                cout << "Lucky: " << lucky << endl;
                cout << "Two: " << two << endl;
                cout << "Three: " << three << endl;
                cout << "Four: " << four << endl;
                cout << "Nothing: " << nothing << endl;
                cout << "SR histo: " << endl;
                for(int i=0;i<10; i++){
                    //if(i%10==0)
                    //    cout << endl;
                    cout << i << " : " << SR_histo[i] << "\t";
                }
                cout << endl << endl;
                cout << "two histo: " << endl;
                int ncr =NCR(clus_size, 2)+1;
                for(int i=0;i<ncr && i<HISTO_LEN; i++){
                    if(i%10==0)
                        cout << endl;
                    cout << i << " - " << pairs_histo[i] << "\t";
                }
                cout << endl << endl;
                cout << "three histo: " << endl;
                ncr =NCR(clus_size, 3)+1;

                for(int i=0;i<ncr&& i<HISTO_LEN; i++){
                    if(i%10==0)
                        cout << endl;
                    cout << i << " - " << three_histo[i] << "\t";
                }
                cout << endl << endl;
                cout << "four histo: " << endl;
                ncr =NCR(clus_size, 4)+1;

                for(int i=0;i<ncr&& i<HISTO_LEN; i++){
                    if(i%10==0)
                        cout << endl;
                    cout << i << " - " << four_histo[i] << "\t";
                }
                cout << endl << endl;
                double res_BMA = (double)(total_dis_tests_BMA) / (i+1);
                res_BMA = res_BMA / original.length();
                cout << "******* BMA ALG ********"<< endl;

                cout << "Error rate BMA: " << res_BMA << endl;
                cout << "Succes rate BMA: " << succes_rate_BMA/(i+1) << endl;
                cout << "SR histo BMA: " << endl;
                for(int i=0;i<10; i++){
                    //if(i%10==0)
                    //    cout << endl;
                    cout << i << " : " << SR_histo_BMA[i] << "\t";
                }
                cout << endl << endl;
            }
            if ( (double)dis/ original.length() == 1){
                
            }
            if(recon.length()<design_len+50 ){
                error_rate_by_length[origSize]+=((double)dis / original.length());
                count_by_length[origSize]+=1;
            }
            cluster.clear();
            if(recon.length()==design_len){
                counterOfGoodLen++;
            }
        }
        cout << endl;
        double res = (double)(total_dis_tests) / (num_of_test);
        res = res / original.length();
        cout << "******* SCS ALG ********"<< endl;
        cout << "Error rate: " << res << endl;
        cout << "Succes rate: " << succes_rate/num_of_test << endl;
        cout << "Lucky: " << lucky << endl;
        cout << "Two: " << two << endl;
        cout << "Three: " << three << endl;
        cout << "Four: " << four << endl;
        cout << "Nothing: " << nothing << endl;
        cout << "SR histo: " << endl;
        for(int i=0;i<10; i++){
            //if(i%10==0)
            //    cout << endl;
            cout << i << " : " << SR_histo[i] << "\t";
        }
        cout << endl << endl;
        cout << "two histo: " << endl;
        int ncr =NCR(clus_size, 2)+1;
        for(int i=0;i<ncr && i<HISTO_LEN; i++){
            if(i%10==0)
                cout << endl;
            cout << i << " - " << pairs_histo[i] << "\t";
        }
        cout << endl << endl;
        cout << "three histo: " << endl;
        ncr =NCR(clus_size, 3)+1;

        for(int i=0;i<ncr&& i<HISTO_LEN; i++){
            if(i%10==0)
                cout << endl;
            cout << i << " - " << three_histo[i] << "\t";
        }
        cout << endl << endl;
        cout << "four histo: " << endl;
        ncr =NCR(clus_size, 4)+1;

        for(int i=0;i<ncr&& i<HISTO_LEN; i++){
            if(i%10==0)
                cout << endl;
            cout << i << " - " << four_histo[i] << "\t";
        }
        cout << endl << endl;

        cout << "******* BMA ALG ********"<< endl;

        double res_BMA = (double)(total_dis_tests_BMA) / (num_of_test);
        res_BMA = res_BMA / original.length();
        cout << "Error rate BMA: " << res_BMA << endl;
        cout << "Succes rate BMA: " << succes_rate_BMA/num_of_test << endl;
        cout << "SR histo BMA: " << endl;
        for(int i=0;i<10; i++){
            //if(i%10==0)
            //    cout << endl;
            cout << i << " : " << SR_histo_BMA[i] << "\t";
        }
        cout << endl << endl;

    }

    }
    catch(exception& e){
        cout<<e.what()<<endl;

        cerr<<e.what()<<endl;

    }

    exit(0);
    return 0;
}

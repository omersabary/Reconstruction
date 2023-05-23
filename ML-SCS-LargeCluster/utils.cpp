//
// Created by User on 26/03/2020.
//

#include "utils.h"
using namespace std;
char intToDNA(int num){
    switch(num){
        case 0 : return 'A';
        case 1 : return 'C';
        case 2 : return 'G';
        case 3 : return 'T';
    }
    return 'N';
}

int DNAtoInt(char c){
    if (c=='A')
        return 0;
    if (c=='C')
        return 1;
    if (c=='G')
        return 2;
    if (c=='T')
        return 3;
    return -1;
}
void reverseStr(string& str)
{
    int n = str.length();

    // Swap character starting from two
    // corners
    for (int i = 0; i < n / 2; i++)
        swap(str[i], str[n - i - 1]);
}



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

template<typename T>
vector<T> create_copy(const vector<T> &vec) {
    vector<T> v(vec);
    return v;
}

int k_string_to_int(string str) {
    vector<char> vec(str.begin(), str.end());
    vector<int> digits;
    for (char c : vec){
        int digit = c - 48;
        digits.push_back(digit);
    }

    int binary_rep [2 * digits.size()];
    for (int i = 0; i < digits.size(); i++){
        int tens = int(digits[i] >= 2);
        int ones = int (digits[i] % 2 == 1);
        binary_rep[2 * i] = tens;
        binary_rep[2 * i + 1]= ones;
    }
    int result = 0;
    int index = 0;
    for (int i = 2 * digits.size() - 1; i >= 0; i--){
        result += binary_rep[i] * pow (2,index);
        index++;
    }
    return result;
}

vector<int> get_kmer(string str, int k) {
    vector<int> k_mer(pow(4,k), 0);
    for (int i = 0; i < str.length() - k + 1; i++){
        string k_str = str.substr(i, k);
        k_mer.at(k_string_to_int(k_str)) += 1;
    }
    return k_mer;
}

vector<vector<int>> cluster_to_kmer(vector<string> &cluster) {
    vector<vector<int>> kmer_cluster;
    for (string strand : cluster){
        kmer_cluster.push_back(get_kmer(strand,2));
    }
    return kmer_cluster;
}

int vector_distance(vector<int> &a, vector<int> &b) {
    int diff = 0;
    for (int i = 0; i < a.size(); i++){
        diff += abs(a.at(i) - b.at(i));
    }
    return diff;
}
int vector_distance(vector<int> &a, vector<int> &b, vector<int> & c) {
    int diff = 0;
    for (int i = 0; i < a.size(); i++){
        diff += abs(a.at(i) - b.at(i)) + abs(a.at(i) - c.at(i)) + abs(b.at(i) - c.at(i));
    }
    return diff;
}
vector<string> two_farthest(vector<tuple<string, vector<int>>>& kmer_cluster){
    int max_dist = -1;
    int total_dist = 0;
    int count = 0;
    int curr_dist = -1;
    vector<int> distances;
    vector<string> result;
    for (int i = 0; i < kmer_cluster.size() - 1; i++){
        for (int j = i + 1; j < kmer_cluster.size(); j++){
            curr_dist = vector_distance(get<1>(kmer_cluster.at(i)), get<1>(kmer_cluster.at(j)));
            if (curr_dist <= 31)
                count ++;
            total_dist += curr_dist;
            distances.push_back(curr_dist);
            if (curr_dist >= max_dist){
                max_dist = curr_dist;
                result.clear();
                result.push_back(get<0>(kmer_cluster.at(i)));
                result.push_back(get<0>(kmer_cluster.at(j)));
            }
        }
    }
    return result;
}

int per_dist_pairs(vector<tuple<string, vector<int>>> &kmer_cluster, int per) {
    vector<int> distances;
    int curr_dist;
    for (int i = 0; i < kmer_cluster.size() - 1; i++) {
        for (int j = i + 1; j < kmer_cluster.size(); j++) {
            curr_dist = vector_distance(get<1>(kmer_cluster.at(i)), get<1>(kmer_cluster.at(j)));;
            distances.push_back(curr_dist);
        }
    }
    sort(distances.begin(), distances.end(), compareValue);
    return distances[distances.size() / per];
}
bool compareValue(int a, int b)
{
    return a > b;
}

vector<string> remove_below_per_pairs(vector<tuple<string, vector<int>>> &kmer_cluster, int per) {
    int curr_dist = -1;
    int median = per_dist_pairs(kmer_cluster, per);
    vector<string> result;
    for (int i = 0; i < kmer_cluster.size() - 1; i++){
        for (int j = i + 1; j < kmer_cluster.size(); j++){
            curr_dist = vector_distance(get<1>(kmer_cluster.at(i)), get<1>(kmer_cluster.at(j)));
            if (curr_dist >= median){
                result.push_back(get<0>(kmer_cluster.at(i)));
                result.push_back(get<0>(kmer_cluster.at(j)));
            }
        }
    }
    return result;
}


int percent_dist_threes(vector<tuple<string, vector<int>>> &kmer_cluster, int per) {
    vector<int> distances;
    int curr_dist;
    for (int i = 0; i < kmer_cluster.size() - 2; i++) {
        for (int j = i + 1; j < kmer_cluster.size() - 1; j++) {
            for (int k = j + 1; k < kmer_cluster.size(); k++ ) {
                curr_dist = vector_distance(get<1>(kmer_cluster.at(i)), get<1>(kmer_cluster.at(j)));;
                distances.push_back(curr_dist);
            }
        }
    }
    sort(distances.begin(), distances.end(), compareValue);
    return distances[distances.size() / per];
}

vector<string> remove_below_per_threes(vector<tuple<string, vector<int>>> &kmer_cluster, int per) {
    int curr_dist = -1;
    int median = percent_dist_threes(kmer_cluster, per);
    vector<string> result;
    for (int i = 0; i < kmer_cluster.size() - 2; i++){
        for (int j = i + 1; j < kmer_cluster.size() - 1 ; j++){
            for (int k = j + 1; k < kmer_cluster.size(); k++) {
                curr_dist = vector_distance(get<1>(kmer_cluster.at(i)), get<1>(kmer_cluster.at(j)));
                if (curr_dist >= median) {
                    result.push_back(get<0>(kmer_cluster.at(i)));
                    result.push_back(get<0>(kmer_cluster.at(j)));
                    result.push_back(get<0>(kmer_cluster.at(k)));
                }
            }
        }
    }
    return result;
}

//
// Created by User on 26/03/2020.
//

#ifndef DNA_RECONSTRUCTION_UTILS_H
#define DNA_RECONSTRUCTION_UTILS_H

using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include "tuple"
#include <string>
#include <iterator>
#include <time.h>
#include <math.h>
#include <random>
#include <chrono>
#include <cstdint>
#include <climits>
#include <fstream>

int DNAtoInt(char c);
char intToDNA(int num);
void reverseStr(string& str);
unsigned int levenstein_dis(const std::string& s1, const std::string& s2);
unsigned int edit_distance(const std::string& s1, const std::string& s2);
template <typename T>
vector<T>create_copy(const vector<T>& vec);

/***
 *
 * @param str: input string
 * @return translate the 4-based string to binary and get the value for hash function
 */
int k_string_to_int(string str);
/***
 * @param str : given string
 * @param k : value for k-mer
 * @return corespoding k-mer rep for the given k
 */
vector<int> get_kmer(string str, int k=2);

/***
 * @param cluster : given cluster of DNA strands
 * @return vector of vectors. each inner vector is a matching kmer of a strand in the clusters
 */
vector<vector<int>> cluster_to_kmer(vector<string> & cluster);

/***
 * @params : a,b two given vectors
 * @return cell wise L1 distance
 */
int vector_distance (vector<int>& a, vector<int>& b);
/**
 *
 * @param a , b ,c kmer vectors
 * @return kmer distance of group
 */
int vector_distance(vector<int> &a, vector<int> &b, vector<int> & c);
/**
 * @param kmer_cluster
 * @return returns the two most farthest strings in the cluster (k-mer wise = L1 for k-mer)
 */

vector<string> two_farthest(vector<tuple<string, vector<int>>>& kmer_cluster);
/***
 * @param kmer_cluster
 * @return the median distance between two strands in the cluster (k-mer wise = L1 for kmer)
 */
int per_dist_pairs(vector<tuple<string, vector<int>>>& kmer_cluster, int per);
bool compareValue(int a, int b);

vector<string> remove_below_per_pairs(vector<tuple<string, vector<int>>>& kmer_cluster, int per);
vector<string> remove_below_per_threes(vector<tuple<string, vector<int>>>& kmer_cluster, int per);

/*/
 * Same as median_dist_pairs but for group size = 3
 */
int percent_dist_threes(vector<tuple<string, vector<int>>> &kmer_cluster, int per);
#endif //DNA_RECONSTRUCTION_UTILS_H

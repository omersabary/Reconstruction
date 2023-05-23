//
// Created by User on 27/03/2020.
//

#ifndef DNA_RECONSTRUCTION_RECONSTRUCT_H
#define DNA_RECONSTRUCTION_RECONSTRUCT_H
using namespace std;
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
#include "utils.h"

/***
 * @param contain : a string to check if the other input is contained in
 * @param contained
 * @return true if contained and false otherwise
 */
bool check_contained(string contain, string contained);
/**
 * @param cluster : the given cluster
 * @param res : empty vector
 * @return removes strands from the cluster that are fully contained in other strands (thus their SCS is the containing strand) in res

 */
void remove_fully_contained(const vector<string>& cluster, vector<string>& res);
/***
 * @param cluster : given cluster of strands the originate from same source
 * @param des_len : the source of the origin strand
 * @param origSize : redundant input for backwards comptibilty
 * @return the original source if found, if not found returns the longest strand that we found and are
 *  sure about it being contained in the source
 */
string reconstruct_big_cluster(const vector<string> &cluster, const int des_len, int &origSize);

#endif //DNA_RECONSTRUCTION_RECONSTRUCT_H

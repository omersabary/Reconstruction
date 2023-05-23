//
// Created by User on 27/03/2020.
//

#include "reconstruct.h"

bool check_contained(string contain, string contained) {
    int size = contained.length();
    int count = 0;
    for (int i = 0; i < contain.length(); i ++ ){
        if (contain[i] == contained[count]){
            count ++;
            if (count == size){
                break;
            }
        }
    }
    return count == size;
}

void remove_fully_contained(const vector <string> &cluster, vector<string>&res) {
//    vector<string>res = vector<string>();
    bool contained_flag = false;
    for(size_t i = cluster.size() - 1; i > 0; i --){
        for (int j = i - 1; j >= 0 ; j --){
            if (check_contained(cluster[j], cluster[i]))
                contained_flag = true;
        }
        if (!contained_flag){
            res.push_back(cluster[i]);
        }
        contained_flag = false;
    }
    res.push_back(cluster[0]);
    reverse(res.begin(), res.end());
}


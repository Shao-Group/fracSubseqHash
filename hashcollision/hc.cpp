/*
  For a list of pairs of strings generated by simulate.cpp, calculate
  the frequency of hash collision using fracSubseqHash with parameter k.

  By: Ke@PSU
  Last edited: 02/08/2023
*/

#include "../util.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define MAXN 50
#define THRESHOLDFACTOR ((1lu<<30)+(1lu<<29))*0

bool getSeed(string &s, int k, RandTableCell* table,
	     DPCell* dp, double threshold, kmer* seed_s){
    int n = s.length();
    const char* ss = s.c_str();
    fillDPTable(ss, n, k, table, dp);
    
    int q = access2d(k+1, n, k);
    double score = fabs(dp[q].min);
    if(dp[q].max > score) score = dp[q].max;
    if(score >= threshold){
	backtrackDPTable(ss, n, k, dp, seed_s);
	return true;
    }
    return false;
}

int main(int argc, const char* argv[]){
    if(argc != 4){
	printf("usage: hc.out stringFile k randTableFile\n");
	return 1;
    }

    int k = atoi(argv[2]);
    double threshold = THRESHOLDFACTOR*k*(20-k)/20.0;

    RandTableCell table[k*ALPHABETSIZE];
    //char table_filename[200];
    //sprintf(table_filename, "%s-k%d", argv[3], k);
    const char* table_filename = argv[3];

    struct stat test_table;
    if(stat(table_filename, &test_table) == 0){//file exists
	loadRandTable(table_filename, k, table);
    }else{
	initRandTable(k, table);
	saveRandTable(table_filename, k, table);
    }

    ifstream fin(argv[1], ifstream::in);

    string s, t;
    kmer seed_s, seed_t;
    int n = 0, invalid = 0, match = 0, no_seed = 0;

    DPCell dp[MAXN * (k+1)];
    
    while(getline(fin, s)){
	getline(fin, t);
	++ n;
	if(k > t.length()) ++ invalid;
	else if(getSeed(s, k, table, dp, threshold, &seed_s) &&
		getSeed(t, k, table, dp, threshold, &seed_t)){
	    if(seed_s == seed_t) ++ match;
	}else ++ no_seed;
    }

    printf("%.2f, %d, %d, %d\n", match * 100.0 / (n - invalid),
	   match, invalid, no_seed);
    
    return 0;
}

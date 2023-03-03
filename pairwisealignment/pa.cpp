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

#define THRESHOLDFACTOR ((1lu<<30)+(1lu<<29))*0.785

struct Seed{
    kmer s;
    int* pos;

    Seed(int k){
	pos = (int*) malloc(sizeof *pos * k);
    }

    Seed(Seed && o){
	s = o.s;
	pos = o.pos;
	o.pos = nullptr;
    }

    Seed& operator=(Seed && o){
	if(this != &o){
	    free(pos);
	    pos = o.pos;
	    s = o.s;
	    o.pos = nullptr;
	}
	return *this;
    }

    ~Seed(){
	if(pos) free(pos);
    }

    friend bool operator<(const Seed &a, const Seed &b){
	if(a.s == b.s) return a.pos[0] < b.pos[0];
	else return a.s < b.s;
    }
};

static inline double getScoreFromDPTable(const int n, const int k,
					 const DPCell* dp){
    int q = access2d(k+1, n, k);
    double score = fabs(dp[q].min);
    if(score < dp[q].max) return dp[q].max;
    else return score;
}

void getSeedsThresholdWithPos(const string &read, const int n, const int k,
			      const RandTableCell* tp, const double threshold,
			      vector<Seed> &seeds){
    int len = read.length();
    int i;
    char cur[n+1];
    double score = threshold;
    DPCell dp[(n+2)*(k+1)];

    for(i=0; i<len-n; i+=1){
	read.copy(cur, n+1, i);
	fillDPTable(cur, n+1, k, tp, dp);
	//printf("called at %d\n", i);

	//get seed from pos i
	score = getScoreFromDPTable(n, k, dp);
	if(score >= threshold){
	    seeds.emplace_back(k);
	    backtrackDPTableWithPos(cur, n, k, dp, &seeds.back().s, i, seeds.back().pos);
	}

	score = getScoreFromDPTable(n+1, k, dp);
	if(score >= threshold){
	    seeds.emplace_back(k);
	    if(!backtrackDPTableWithPos(cur, n+1, k, dp, &seeds.back().s, i, seeds.back().pos)){//first char not used
		++i; //skip recalculation of next position
	    }else{
		seeds.pop_back();
	    }
	}else{//score less than threshold even with extra column, actual score can only be lower
	    ++i;
	}
    }

    if(i == len - n){
	read.copy(cur, n, i);
	fillDPTable(cur, n, k, tp, dp);
	//printf("called at %d\n", i);
	score = getScoreFromDPTable(n, k, dp);
	if(score >= threshold){
	    seeds.emplace_back(k);
	    backtrackDPTableWithPos(cur, n, k, dp, &seeds.back().s, i, seeds.back().pos);
	}
    }
}

void getMatches(const string &s, const string &t,
		const int n, const int k,
		const RandTableCell* tp, const double threshold,
		const int* align, double* result){
    vector<Seed> s_seeds, t_seeds;
    getSeedsThresholdWithPos(s, n, k, tp, threshold, s_seeds);
    sort(s_seeds.begin(), s_seeds.end());
    getSeedsThresholdWithPos(t, n, k, tp, threshold, t_seeds);
    sort(t_seeds.begin(), t_seeds.end());


    int s_len = s.length(), t_len = t.length();
    int s_total = s_seeds.size(), t_total = t_seeds.size();
    bool s_cover[2][s_len];//0 - false, 1 - true
    bool t_cover[2][t_len];
    memset(s_cover, 0, sizeof **s_cover * 2 * s_len);
    memset(t_cover, 0, sizeof **t_cover * 2 * t_len);

    int is = 0, it = 0, inner, i;
    int correct_ct;
    int match_total = 0, match_true = 0;
    int classification;
    kmer cur;
    while(is < s_total && it < t_total){
	cur = s_seeds[is].s;
	if(cur == t_seeds[it].s){
	    inner = it;
	    while(inner < t_total && cur == t_seeds[inner].s){
		++ match_total;
		correct_ct = 0;
		for(i=0; i<k; ++i){
		    if(align[s_seeds[is].pos[i]] == t_seeds[inner].pos[i]){
			++ correct_ct;
		    }
		}
		if(correct_ct >= (k>>1)){//>=50%
		    ++ match_true;
		    classification = 1;
		}else classification = 0;

		for(i=0; i<k; ++i){
		    s_cover[classification][s_seeds[is].pos[i]] = true;
		}
		for(i=0; i<k; ++i){
		    t_cover[classification][t_seeds[inner].pos[i]] = true;
		}
		++ inner;
	    }
	    ++ is;
	}else if(cur < t_seeds[it].s){
	    ++ is;
	}else{
	    ++ it;
	}
    }
    
    result[0] = match_total;
    result[1] = match_true;

    int coverage_ct = 0;
    for(i=0; i<s_len; ++i){
	if(s_cover[0][i]) ++ coverage_ct;
    }
    for(i=0; i<t_len; ++i){
	if(t_cover[0][i]) ++ coverage_ct;
    }

    result[3] = ((double)coverage_ct)/(s_len+t_len);

    coverage_ct = 0;
    for(i=0; i<s_len; ++i){
	if(s_cover[1][i]) ++ coverage_ct;
    }
    for(i=0; i<t_len; ++i){
	if(t_cover[1][i]) ++ coverage_ct;
    }

    result[2] = ((double)coverage_ct)/(s_len+t_len);
    
}

int main(int argc, const char* argv[]){
    if(argc != 5){
	printf("usage: pa.out stringFile n k randTableFile\n");
	return 1;
    }

    int n = atoi(argv[2]);
    int k = atoi(argv[3]);
    double threshold = THRESHOLDFACTOR*k;//*(20-k)/20.0;

    RandTableCell table[k*ALPHABETSIZE];
    //char table_filename[200];
    //sprintf(table_filename, "%s-k%d", argv[3], k);
    const char* table_filename = argv[4];

    struct stat test_table;
    if(stat(table_filename, &test_table) == 0){//file exists
	loadRandTable(table_filename, k, table);
    }else{
	initRandTable(k, table);
	saveRandTable(table_filename, k, table);
    }

    ifstream fin(argv[1], ifstream::in);

    string s, t;
    int len, i, pair_ct = 0;
    double cur_result[4]; //num_match, num_true_match, true_coverage, false_coverage
    double result[4];

    memset(result, 0, sizeof *result *4);
    
    int* align;
    
    
    while(getline(fin, s)){
	getline(fin, t);

	len = s.length();
	align = (int*) malloc(sizeof *align * len);
	for(i=0; i<len; ++i){
	    //fscanf(fin, "%d ", &align[i]);
	    fin >> align[i];
	}
	
	++ pair_ct;

	getMatches(s, t, n, k, table, threshold, align, cur_result);

	for(i=0; i<4; ++i){
	    result[i] += cur_result[i];
	}
	
	free(align);
	getline(fin, s);
    }

    printf("%d/%d, %.2f, %.2f, %.4f, %.4f, %.4f\n", n, k,
	   result[0], result[1], result[1]/result[0],
	   result[2]/pair_ct, result[3]/pair_ct);
    
    return 0;
}

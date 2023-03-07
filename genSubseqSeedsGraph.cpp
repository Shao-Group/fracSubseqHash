/*
  Given a set of reads, extract their FracSubseqHash seeds as vertices.
  Two seeds a and b are connected by a directed edge if they are adjacent
  seeds (in this order) on some read.

  Each read is therefore represented by a path in the resulting graph.

  After all reads are processed, remove nodes (seeds) that only appear
  in one read.

  Output the graph in dot format.
  
  By: Ke@PSU
  Last edited: 03/06/2023
*/

#include "util.h"
#include "SeedsGraph.hpp"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <functional>

using namespace std;

#define NUMTHREADS 15
#define EXPECTEDVALUE ((1lu<<30)+(1lu<<29))
#define THRESHOLDFACTOR 0.785

//from brewer color scheme spectral11
const char* COLORS[] = {"#910142", "#5e4fa2", "#f46d43",
	"#66c2a5", "#fee08b", "#e6f598", "#3288bd",
	"#d53e4f", "#abdda4", "#fdae61", "#ffffbf"};

typedef SeedsGraph<kmer> Graph;
typedef Graph::Node Node;

string kmerToString(const kmer& x, int k, char* buf){
    decode(x, k, buf);
    return string(buf);
}

struct ReadPath{
    size_t read_idx;
    Node *head, *tail;

    ReadPath(size_t read_idx):read_idx(read_idx),
			      head(nullptr), tail(nullptr){};
};

struct Read{
    string seq;
    size_t idx;

    Read(string&& s, size_t i): seq(move(s)), idx(i) {};
    Read(Read&& o): seq(move(o.seq)), idx(exchange(o.idx, 0)) {};
};

class SeedFactory{
    const int n;
    const int k;
    const RandTableCell* table;
    const double threshold;
    Graph &graph;
    
    queue<Read> jobs;
    vector<thread> minions;
    bool done;
    mutex door;
    condition_variable trumpet;

    Node* storeSeedWithPosInGraph(const kmer seed, const size_t read_idx,
				  const size_t cur_pos, size_t* prev_pos,
				  Node* prev, Graph& g);

    double getScoreFromDPTable(const int n, const int k, const DPCell* dp);
    
    void getAndSaveSubseqSeeds(const Read &r);
    void atWork(int x);

public:
    SeedFactory(const int n, const int k, const RandTableCell* table,
		const double threshold, Graph& g):
	n(n), k(k), table(table), threshold(threshold),
	graph(g), done(false){

	minions.reserve(NUMTHREADS);
	for(int i=0; i<NUMTHREADS; ++i){
	    minions.emplace_back(bind(&SeedFactory::atWork, this, i));
	}
    }

    ~SeedFactory(){
	unique_lock<mutex> lock(door);
	done = true;
	lock.unlock();
	trumpet.notify_all();

	for(auto& x : minions){
	    x.join();
	}
    }

    void addJob(string&& r, size_t idx){
	unique_lock<mutex> lock(door);
	jobs.emplace(move(r), idx);
	trumpet.notify_one();
    }
};

int main(int argc, const char * argv[])
{
    if(argc != 5){
	printf("usage: genSubseqGraph.out readFile n k randTableFile\n");
	return 1;
    }

    int n = atoi(argv[2]);
    int k = atoi(argv[3]);
    double threshold = THRESHOLDFACTOR * EXPECTEDVALUE * k;

    //load table
    RandTableCell table[k*ALPHABETSIZE];
    const char* table_filename = argv[4];

    struct stat test_table;
    if(stat(table_filename, &test_table) == 0){//file exists
	loadRandTable(table_filename, k, table);
    }else{
	initRandTable(k, table);
	saveRandTable(table_filename, k, table);
    }

    //input reads and process
    ifstream fin(argv[1], ifstream::in);

    Graph g;
    //vector<ReadPath> paths;
    SeedFactory factory(n, k, table, threshold, g);
    
    string read;
    size_t read_idx = 0;
    
    while(fin.get() == '>'){
	//fin >> read_idx;
	//skip the header
	fin.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(fin, read);
	++ read_idx;
	factory.addJob(move(read), read_idx);
    }

    //only keep reads that appear on multiple distinct reads
    g.removeUniqSeeds();

    //output
    char output_filename[200];
    int output_len = strstr(argv[1], ".efa") - argv[1];
    int tablename_st = strlen(table_filename) - 1;
    for(; tablename_st>=0; --tablename_st){
	if(table_filename[tablename_st] == '/'){
	    break;
	}
    }
    ++tablename_st;
    output_len = sprintf(output_filename, "%.*s-%s-t%f.graph",
			 output_len, argv[1],
			 table_filename+tablename_st,
			 THRESHOLDFACTOR);
    
    g.saveGraph(output_filename);
    
    return 0;
}


/*** implementation of SeedFactory functions ***/

inline Node* SeedFactory::storeSeedWithPosInGraph(
    kmer seed, const size_t read_idx, const size_t cur_pos,
    size_t* prev_pos, Node* prev, Graph& g){
    //avoid self loops
    if(prev && prev->seed == seed) return prev;
    
    Node* cur = g.addNode(seed);

    if(prev){
	prev->addNext(read_idx, *prev_pos, cur);
	cur->addPrev(read_idx, cur_pos, prev);
    }
    *prev_pos = cur_pos;
    return cur;
}

inline double SeedFactory::getScoreFromDPTable(const int n, const int k,
					       const DPCell* dp){
    int q = access2d(k+1, n, k);
    double score = fabs(dp[q].min);
    if(score < dp[q].max) return dp[q].max;
    else return score;
}

void SeedFactory::atWork(int x){
    unique_lock<mutex> lock(door, defer_lock);
    while(true){
	lock.lock();
	while(!done && jobs.empty()){
	    trumpet.wait(lock);
	}
	if(!jobs.empty()){
	    Read r = move(jobs.front());
	    jobs.pop();
	    lock.unlock();
	    getAndSaveSubseqSeeds(r);
	}else{
	    return;
	}
    }
}

void SeedFactory::getAndSaveSubseqSeeds(const Read &r){
    int len = r.seq.length();
    
    int i;
    char cur[n+1];
    kmer seed;
    double score;
    Node *prev = nullptr;
    size_t prev_pos;

    //calculate an extra column, can skip next position if score at
    //[n+1][k] does not reach threshold; otherwise does not need recalculation
    //if backtrack from [n+1][k] does not use first char
    DPCell dp[(n+2)*(k+1)];
	
    for(i=0; i<len-n; i+=1){
	r.seq.copy(cur, n+1, i);
	fillDPTable(cur, n+1, k, table, dp);
	//printf("called at %d\n", i);

	//get seed from pos i
	score = getScoreFromDPTable(n, k, dp);
	if(score >= threshold){
	    backtrackDPTable(cur, n, k, dp, &seed);
	    //add to graph
	    prev = storeSeedWithPosInGraph(seed, r.idx, i,
					   &prev_pos, prev, graph);
	}

	score = getScoreFromDPTable(n+1, k, dp);
	if(score >= threshold){
	    if(!backtrackDPTable(cur, n+1, k, dp, &seed)){//first char not used
		++i; //skip recalculation of next position
		//add to graph
		prev = storeSeedWithPosInGraph(seed, r.idx, i,
					       &prev_pos, prev, graph);
	    }
	}else{//score less than threshold even with extra column, actual score can only be lower
	    ++i;
	}
    }

    //handle last seed, either it's never calculated or the previous
    //iteration did not work (i.e., score >= threshold but used 1st char)
    if(i == len - n){
	r.seq.copy(cur, n, i);
	fillDPTable(cur, n, k, table, dp);
	//printf("called at %d\n", i);
	score = getScoreFromDPTable(n, k, dp);
	if(score >= threshold){
	    backtrackDPTable(cur, n, k, dp, &seed);
	    prev = storeSeedWithPosInGraph(seed, r.idx, i,
					   &prev_pos, prev, graph);
	}
    }
}

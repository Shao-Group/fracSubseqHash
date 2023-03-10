/*
  Given a fasta read file, for each read, generate subseq seeds and save them
  in a file. The seeds are generated with parameters n, k, and a set of
  random tables. The given table is loaded if it exists, otherwise a new set
  of random tables is generated and saved with the provided name. (See the 
  manuscript for more info on the algorithm of generating subseq seeds).

  Only seeds with scores at least the threshold are kept.
  
  The seeds are generated in parallel with NUMTHREADS threads.

  Seeds for each read is stored in a separate file. The files are meant to be
  loaded to generate a seed graph where nodes are seeds and two seeds are
  connected if they are obtained from consecutive windows (ignoring windows
  where no seed pass the threshold) on some read.

  By: Ke@PSU
  Last edited: 03/10/2023
*/

#include "util.h"
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
    const char* output_dir;
    const int dir_len;
    
    queue<Read> jobs;
    vector<thread> minions;
    bool done;
    mutex door;
    condition_variable trumpet;
    
    void getAndSaveSubseqSeeds(const Read &r){
	vector<Seed> seeds_list;
	getSubseqSeedsThreshold(r.seq, n, k, table, threshold, seeds_list);

	char output_filename[200];
	sprintf(output_filename, "%.*s/%zu.subseqseed",
		dir_len, output_dir, r.idx);

	saveSubseqSeeds(output_filename, seeds_list);
    }
    
    void atWork(int x){
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

public:
    SeedFactory(const int n, const int k, const RandTableCell* table,
		const double threshold, const char* output_dir, const int dir_len):
	n(n), k(k), table(table), threshold(threshold),
	output_dir(output_dir), dir_len(dir_len), done(false){

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
	printf("usage: genSubseqSeeds.out readFile n k randTableFile\n");
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


    //output directory
    char output_dir[200];
    int dir_len = strstr(argv[1], ".efa") - argv[1];
    int tablename_st = strlen(table_filename) - 1;
    for(; tablename_st>=0; --tablename_st){
	if(table_filename[tablename_st] == '/'){
	    break;
	}
    }
    ++tablename_st;
    dir_len = sprintf(output_dir, "%.*s-seeds-%s-t%f",
		      dir_len, argv[1],
		      table_filename+tablename_st,
		      THRESHOLDFACTOR);

    mkdir(output_dir, 0744);

    //input reads and process
    ifstream fin(argv[1], ifstream::in);

    {
	SeedFactory factory(n, k, table, threshold, output_dir, dir_len);   
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
    }

    printf("%s %s %d %d %f %s done\n", argv[0], argv[1], n, k, threshold, argv[4]);
    
    return 0;
}

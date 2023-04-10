/*
  Given a set of seed files (readable by loadSubseqSeedsPos), 
  output pairs of reads with the number of unique seeds they share.
  To avoid reporting transitive overlapping pairs, for each seed, 
  reads containing it are sorted in reverse order according to the 
  position of the seed. Only adjacent pairs in this order are counted. 
  
  By: Ke@PSU
  Last edited: 04/07/2023
*/

#include "util.h"
#include <sys/stat.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

struct Occurrence{
    int read_id;
    unsigned int pos;

    Occurrence(): read_id(0), pos(0){}
    Occurrence(const int id, const unsigned int pos): read_id(id), pos(pos){}
    Occurrence(const Occurrence& o): read_id(o.read_id), pos(o.pos){}
    bool operator < (const Occurrence& x) const{
	return pos > x.pos;
    }
};

void loadSubseqSeedsPos(const char* filename, const int read_id,
			map<kmer, vector<Occurrence> > &all_seeds){
    FILE* fin = fopen(filename, "rb");
    Seed s;
    while(fread(&s, sizeof(s), 1, fin) == 1){
	auto result = all_seeds.emplace(s.v, 1);
	if(result.second == false){
	    //if(result.first->second.back().read_id < read_id){
	    result.first->second.emplace_back(read_id, s.pos);
	    //}
	}else{
	    Occurrence* x = &(result.first->second[0]);
	    x->read_id = read_id;
	    x->pos = s.pos;
	}
    }
    if(ferror(fin)){
	fprintf(stderr, "Error reading %s\n", filename);
    }
    fclose(fin);
}


int main(int argc, const char * argv[])    
{   
    if(argc != 3){
	printf("usage: overlapBySeedsPos.out seedsDir numFiles\n");
	return 1;
    }

    int n = atoi(argv[2]);
    //int threshold = atoi(argv[3]);


    char filename[500];
    int i = strlen(argv[1]);
    memcpy(filename, argv[1], i);
    if(filename[i-1] != '/'){
	filename[i] = '/';
	++i;
    }
    
    map<kmer, vector<Occurrence> > all_seeds;
    int j;
    
    struct stat test_file;
    for(j=1; j<=n; j+=1){
	sprintf(filename+i, "%d.subseqseed", j);
	if(stat(filename, &test_file) != 0){//seed file does not exist
	    fprintf(stderr, "Stopped, cannot find file %d.subseqseed\n", j);
	    break;
	}
	loadSubseqSeedsPos(filename, j, all_seeds);
    }

    sprintf(filename+i, "overlapPos-n%d.all-pair", n);

    Table share_ct(n);
    Table share_ct_rev(n);

    int a, b, c;
    
    for(auto& seed : all_seeds){
	c = seed.second.size();
	if(c > 1){
	    sort(seed.second.begin(), seed.second.end());
	    a = seed.second[0].read_id;
	    for(i=1; i<c; ++i){
		b = seed.second[i].read_id;
		if(a < b) ++ share_ct.access(a, b);
		else if (b < a) ++ share_ct_rev.access(b, a);
		a = b;
		// if(share_ct.access(a, b) == threshold){
		// fprintf(fout, "%d %d\n", a, b);
		// }
	    }
	}
    }

    share_ct.saveNoneZeroEntries(filename);
    share_ct_rev.saveNoneZeroEntries(filename, "a", true);

    //sort the output
    char cmd[2000];
    //sprintf(cmd, "sort -k1g,2 -k2g,3 -o %s %s", filename, filename);
    sprintf(cmd, "bash -c 'sort -k1g,2 -k2g,3 -o %s{,}'", filename);
    return system(cmd);
}

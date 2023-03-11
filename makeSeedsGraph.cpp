/*
  Load a set of seed files, one for each read, as vertices of the graph.
  Two seeds a and b are connected by a directed edge if they are adjacent
  seeds (in this order) on some read.

  Each read is therefore represented by a path in the resulting graph.

  After all reads are processed, remove nodes (seeds) that only appear
  in one read.

  Output the graph in dot format.
  
  By: Ke@PSU
  Last edited: 03/10/2023
*/

#include "util.h"
#include "SeedsGraph.hpp"
#include <sys/stat.h>
#include <iostream>
#include <fstream>

using namespace std;

typedef SeedsGraph<kmer> Graph;
typedef Graph::Node Node;

string kmerToString(const kmer& x, unsigned int k, char* buf){
    decode(x, k, buf);
    return string(buf);
}

struct ReadPath{
    size_t read_idx;
    Node *head, *tail;

    ReadPath(const size_t read_idx):read_idx(read_idx),
				    head(nullptr), tail(nullptr){};
};

inline Node* storeSeedWithPosInGraph(
    kmer seed, const size_t read_idx, const size_t cur_pos,
    size_t* prev_pos, Node* prev, Graph& g){
    //avoid self loops -- already done at seed generation
    //if(prev && prev->seed == seed) return prev;
    
    Node* cur = g.addNode(seed);

    if(prev){
	prev->addNext(read_idx, *prev_pos, cur);
	cur->addPrev(read_idx, cur_pos, prev);
    }
    *prev_pos = cur_pos;
    return cur;
}

void loadSubseqSeeds(const char* filename, const size_t read_idx,
		     Graph& g, ReadPath& p){
    FILE* fin = fopen(filename, "rb");
    Seed s;
    size_t prev_pos;
    Node* prev=nullptr;

    if(fread(&s, sizeof(s), 1, fin) == 1){
	//first node
	prev = g.addNode(s.v);
	prev_pos = s.pos;
	p.head = p.tail = prev;

	//following nodes
	while(fread(&s, sizeof(s), 1, fin) == 1){
	    prev = storeSeedWithPosInGraph(s.v, read_idx, s.pos,
					   &prev_pos, prev, g);
	    p.tail = prev;
	}
    }
    if(ferror(fin)){
	fprintf(stderr, "Error reading %s/n", filename);
    }
    fclose(fin);
}

void saveGraphToDot(const char* filename, const unsigned int k,
		    const Graph& g, const vector<ReadPath>& paths){
    ofstream fout(filename, ofstream::out);
    char buf[k+1];
    buf[k] = '\0';
    
    fout << "digraph{" << endl;
    g.printNodesInDot(fout, kmerToString, k, buf);
    g.printEdgesInDot(fout);
    
    //add head and tail nodes for each path
    for(const ReadPath& p : paths){
	fout << "st" << p.read_idx << " [label=\"Read " << p.read_idx
	     << " head\"];" << endl;
	fout << "ed" << p.read_idx << " [label=\"Read " << p.read_idx
	     << " tail\"];" << endl;
	//add edge from head to first node
	fout << "st" << p.read_idx << " -> n" << p.head->id << ";" << endl;

	//add edge from last node to tail
	fout << "n" << p.tail->id << " -> ed" << p.read_idx << ";" << endl;
    }
    
    fout << "} //end of graph" << endl;
}

int main(int argc, const char * argv[])
{
    if(argc != 4){
	printf("usage: makeSeedsGraph.out seedsDir k numFiles\n");
	return 1;
    }

    unsigned int n = atoi(argv[3]);
    unsigned int k = atoi(argv[2]);

    char filename[500];
    unsigned int dir_len = strlen(argv[1]);
    memcpy(filename, argv[1], dir_len);
    if(filename[dir_len-1] != '/'){
	filename[dir_len] = '/';
	++dir_len;
    }
    
    Graph g;
    vector<ReadPath> paths;
    paths.reserve(n);
    size_t j;
    
    struct stat test_file;
    for(j=1; j<=n; j+=1){
	sprintf(filename+dir_len, "%zu.subseqseed", j);
	if(stat(filename, &test_file) != 0){//seed file does not exist
	    fprintf(stderr, "Stopped, cannot find file %zu.subseqseed\n", j);
	    break;
	}
	paths.emplace_back(j);
	loadSubseqSeeds(filename, j, g, paths.back());
    }

    //only keep reads that appear on multiple distinct reads
    g.removeUniqSeeds();

    //output
    sprintf(filename+dir_len, "overlap-n%d-graph.dot", n);
    
    saveGraphToDot(filename, k, g, paths);
    
    return 0;
}


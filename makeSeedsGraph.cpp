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

inline Node* storeSeedWithPosInGraph(
    kmer seed, const size_t read_idx, const size_t cur_pos,
    const size_t cur_span,
    size_t* prev_pos, Node* prev, Graph& g){
    //avoid self loops -- already done at seed generation
    //if(prev && prev->seed == seed) return prev;
    
    Node* cur = g.addNode(seed);

    if(prev){
	prev->addNext(read_idx, *prev_pos, cur_span, cur);
	cur->addPrev(read_idx, cur_pos, cur_span, prev);
    }
    *prev_pos = cur_pos;
    return cur;
}

void loadSubseqSeeds(const char* filename, const size_t read_idx, Graph& g){
    FILE* fin = fopen(filename, "rb");
    Seed s;
    size_t prev_pos;
    Node* prev=nullptr;
    Node *head=nullptr, *tail=nullptr;

    if(fread(&s, sizeof(s), 1, fin) == 1){
	//first node
	prev = g.addNode(s.v);
	prev_pos = s.pos;
	head = tail = prev;

	//following nodes
	while(fread(&s, sizeof(s), 1, fin) == 1){
	    prev = storeSeedWithPosInGraph(s.v, read_idx, s.pos, s.span,
					   &prev_pos, prev, g);
	    tail = prev;
	}

	g.addReadPath(read_idx, head, tail);
    }
    if(ferror(fin)){
	fprintf(stderr, "Error reading %s/n", filename);
    }
    fclose(fin);
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
    
    Graph g(n);
    size_t j;

    //load all seeds
    struct stat test_file;
    for(j=1; j<=n; j+=1){
	sprintf(filename+dir_len, "%zu.subseqseed", j);
	if(stat(filename, &test_file) != 0){//seed file does not exist
	    fprintf(stderr, "Stopped, cannot find file %zu.subseqseed\n", j);
	    break;
	}
	loadSubseqSeeds(filename, j, g);
    }

    //only keep reads that appear on multiple distinct reads
    g.removeUniqSeeds();

    //output to dot file
    sprintf(filename+dir_len, "overlap-n%d-graph.dot", n);
    char buf[k+1];
    buf[k] = '\0';
    g.saveGraphToDot(filename, kmerToString, k, buf);

    //save graph to binary file
    sprintf(filename+dir_len, "overlap-n%d.graph", n);
    g.saveGraph(filename);

    //test save and load graph produce an identical copy
    /*
    Graph g2;
    g2.loadGraph(filename);
    sprintf(filename+dir_len, "overlap-n%d-graph.dot.cp", n);
    g2.saveGraphToDot(filename, kmerToString, k, buf);
    */
    
    return 0;
}


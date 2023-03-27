/*
  Read in a previously generated seed graph
  and output in dot format (for format changes).

  By: Ke@PSU
  Last edited: 03/27/2023
*/

#include "util.h"
#include "SeedsGraph.hpp"
#include <sys/stat.h>
#include <iostream>
#include <fstream>

using namespace std;

typedef SeedsGraph<kmer> Graph;

string kmerToString(const kmer& x, unsigned int k, char* buf){
    decode(x, k, buf);
    return string(buf);
}

int main(int argc, const char * argv[])
{
    if(argc != 3){
	printf("usage: reloadSeedsGraph.out graphFile k\n");
	return 1;
    }

    unsigned int k = atoi(argv[2]);

    Graph g;
    g.loadGraph(argv[1]);

    char filename[200];
    sprintf(filename, "%s-withloc.dot", argv[1]);

    char buf[k+1];
    buf[k] = '\0';
    g.saveGraphToDot(filename, kmerToString, k, buf);
    
    return 0;
}


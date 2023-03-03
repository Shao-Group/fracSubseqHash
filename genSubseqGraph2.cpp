/*
  Given a set of reads, extract their FracSubseqHash seeds as vertices.
  Two seeds a and b are connected by a directed edge if they are adjacent
  seeds (in this order) on some read.

  Each read is therefore represented by a path in the resulting graph.

  Output the graph in dot format.
  
  By: Ke@PSU
  Last edited: 02/25/2023
*/

#include "util.h"
#include "SeedGraph.hpp"
#include <sys/stat.h>
#include <iostream>
#include <fstream>

using namespace std;

#define EXPECTEDVALUE ((1lu<<30)+(1lu<<29))
#define THRESHOLDFACTOR 0.785

//from brewer color scheme spectral11
const char* COLORS[] = {"#910142", "#5e4fa2", "#f46d43",
	"#66c2a5", "#fee08b", "#e6f598", "#3288bd",
	"#d53e4f", "#abdda4", "#fdae61", "#ffffbf"};

typedef SeedGraph<kmer> Graph;
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

static inline Node* storeSeedWithPosInGraph(
    const kmer seed, const size_t read_idx, const size_t cur_pos,
    size_t* prev_pos, Node* prev, Graph& g, ReadPath& path){

    if(prev && prev->key == seed) return prev;
    Node* cur = g.addNode(seed);
    cur->addPath(read_idx, cur_pos, *prev_pos, prev);
    //update path record
    if(path.head == nullptr) path.head = cur;
    path.tail = cur;
    *prev_pos = cur_pos;
    return cur;
}

static inline double getScoreFromDPTable(const int n, const int k,
					 const DPCell* dp){
    int q = access2d(k+1, n, k);
    double score = fabs(dp[q].min);
    if(score < dp[q].max) return dp[q].max;
    else return score;
}
			      

void addToGraph(const string &read, const size_t read_idx, const int n,
		const int k, const RandTableCell *tp,
		const double threshold, Graph &g,
	        ReadPath& path){
    int len = read.length();
    int i;
    char cur[n+1];
    kmer seed;
    double score = threshold;
    Node *prev = nullptr;
    size_t prev_pos;
    
    //calculate an extra column, can skip next position if score at
    //[n+1][k] does not reach threshold; otherwise does not need recalculation
    //if backtrack from [n+1][k] does not use first char
    DPCell dp[(n+2)*(k+1)];

    for(i=0; i<len-n; i+=1){
	read.copy(cur, n+1, i);
	fillDPTable(cur, n+1, k, tp, dp);
	//printf("called at %d\n", i);

	//get seed from pos i
	score = getScoreFromDPTable(n, k, dp);
	if(score >= threshold){
	    backtrackDPTable(cur, n, k, dp, &seed);
	    //add to graph
	    prev = storeSeedWithPosInGraph(seed, read_idx, i,
					   &prev_pos, prev, g, path);
	}

	score = getScoreFromDPTable(n+1, k, dp);
	if(score >= threshold){
	    if(!backtrackDPTable(cur, n+1, k, dp, &seed)){//first char not used
		++i; //skip recalculation of next position
		//add to graph
		prev = storeSeedWithPosInGraph(seed, read_idx, i,
					       &prev_pos, prev, g, path);
	    }
	}else{//score less than threshold even with extra column, actual score can only be lower
	    ++i;
	}
    }

    //handle last seed, either it's never calculated or the previous
    //iteration did not work (i.e., score >= threshold but used 1st char)
    if(i == len - n){
	read.copy(cur, n, i);
	fillDPTable(cur, n, k, tp, dp);
	//printf("called at %d\n", i);
	score = getScoreFromDPTable(n, k, dp);
	if(score >= threshold){
	    backtrackDPTable(cur, n, k, dp, &seed);
	    prev = storeSeedWithPosInGraph(seed, read_idx, i,
					   &prev_pos, prev, g, path);
	}
    }
}

int main(int argc, const char * argv[])
{
    if(argc != 5){
	printf("usage: genSubseqGraph.out readFile n k randTableFile\n");
	return 1;
    }

    int n = atoi(argv[2]);
    int k = atoi(argv[3]);
    double threshold = THRESHOLDFACTOR * EXPECTEDVALUE * k;

    RandTableCell table[k*ALPHABETSIZE];
    const char* table_filename = argv[4];

    struct stat test_table;
    if(stat(table_filename, &test_table) == 0){//file exists
	loadRandTable(table_filename, k, table);
    }else{
	initRandTable(k, table);
	saveRandTable(table_filename, k, table);
    }

    ifstream fin(argv[1], ifstream::in);

    char output_filename[200];
    int output_len = strstr(argv[1], ".efa") - argv[1];
    int tablename_st = strlen(table_filename) - 1;
    for(; tablename_st>=0; --tablename_st){
	if(table_filename[tablename_st] == '/'){
	    break;
	}
    }
    ++tablename_st;
    output_len = sprintf(output_filename, "%.*s-%s-t%f.delaynode.dot",
			 output_len, argv[1],
			 table_filename+tablename_st,
			 THRESHOLDFACTOR);

    string read;
    size_t read_idx = 0;
    Graph g;
    vector<ReadPath> paths;
    
    while(fin.get() == '>'){
	fin >> read_idx;
	//skip the header
	fin.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(fin, read);
	//++ read_idx;
	ReadPath p(read_idx);
	addToGraph(read, read_idx, n, k, table, threshold, g, p);
	paths.push_back(p);
    }

    Node* cur;
    size_t cur_pos;
    char buf[k+1];
    buf[k]='\0';
    /*
    for(ReadPath p : paths){
	cout << "path " << p.read_idx <<endl;
	cur = p.head;
	cur_pos = 0;
	while(cur){
	    cout<< "# edges: " << cur->paths.size() << " " << cur->toString(kmerToString, k, buf) << endl; 
	    auto local = cur->getPathLowerBD(p.read_idx, cur_pos);
	    if(local != cur->paths.end()){
		cur = local->second.next;
		cur_pos = local->first.pos;
	    }else cur = nullptr;
	}
    }
    */

    ofstream fout(output_filename, ofstream::out);
    fout << "digraph{" << endl;
    fout << "graph[compound=true];" << endl;
    //g.printNodesInDot(fout, kmerToString, k, buf);

    size_t cur_color = 0;
    size_t cur_node_ct = 0;
    
    for(const ReadPath& p : paths){
	fout << "subgraph cluster_read" << p.read_idx << " {" << endl;
	fout << "edge [color=\"" << COLORS[cur_color] << "\"];" << endl;
	//add head and tail nodes for each path
	fout << "st" << p.read_idx << " [label=\"Read " << p.read_idx
	     << " head\"];" << endl;
	fout << "ed" << p.read_idx << " [label=\"Read " << p.read_idx
	     << " tail\"];" << endl;
	++ cur_color;
	cur = p.head;
	cur_pos = 0;
	if(cur->id2 == 0){//new node
	    cur->id2 = (++ cur_node_ct);
	    fout << cur->toString2(false, kmerToString, k, buf);
	}
	//add edge from head to first node
	fout << "st" << p.read_idx << " -> n" << cur->id2 << ";" << endl;
	
	while(cur){
	    auto local = cur->getPathLowerBD(p.read_idx, cur_pos);
	    //if(local != cur->paths.end()){
	    //local is still valid when cur == p.tail with
	    //first being the location (on read) of p.tail
	    if(cur != p.tail){
		if(local->second.next->id2 == 0){
		    local->second.next->id2 = (++ cur_node_ct);
		    fout << local->second.next->toString2(false, kmerToString, k, buf);
		}
		fout << "n" << cur->id2 << " -> n" << local->second.next->id2
		     << " [taillabel=" << local->first.pos << "];" << endl;
		cur = local->second.next;
	    }else cur = nullptr;
	    cur_pos = local->first.pos;
	}

	//add edge from last node to tail
	if(p.tail->id2 == 0){
	    p.tail->id2 = (++ cur_node_ct);
	    fout << p.tail->toString2(false, kmerToString, k, buf);
	}
	fout << "n" << p.tail->id2 << " -> ed" << p.read_idx
	     << " [taillabel=" << cur_pos << "];" << endl;
	fout << "}; // end of read "<< p.read_idx << endl;
    }
    
    fout << "} //end of graph" << endl;
    
    return 0;
}

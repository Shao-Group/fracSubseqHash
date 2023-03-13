/*
  Utility functions for the FracSubseqHash project.

  By: Ke@PSU
  Last edited: 02/06/2023
*/

#ifndef _UTIL_H
#define _UTIL_H 1

#include <cstdio>
#include <cstring>
#include <map>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
//#include <queue>


/*
  Each k-mer is represented by a long unsigned int 
  with the encoding A-00, C-01, G-10, T-11.
*/
typedef __uint128_t kmer;

struct Seed{
    kmer v;
    unsigned int pos;
    //number of consecutive windows that all produce this seed
    unsigned int span;
    Seed(const kmer v, const unsigned int pos): v(v), pos(pos), span(1) {};
    Seed() {};
};


#define ALPHABETSIZE 4
const char ALPHABET[ALPHABETSIZE] = {'A', 'C', 'G', 'T'};

static inline int alphabetIndex(char c){
    return 3 & ((c>>2) ^ (c>>1));
}

/*
  Encode the string representation of a k-mer.
*/
kmer encode(const char* s, const int k);

/*
  Decode an k-mer into its string representation.
  If str is not null, it is used to store the resulting string; 
  otherwise a new char array is allocated.
*/
char* decode(const kmer enc, const int k, char* str);


/*
  Following definitions and functions are related to the dynamic programming 
  procedure for computing a smallest bmer out of a given kmer.
  See the manuscript for more details.
*/


/*
  A cell in the DP table
*/
typedef struct{
    double max, min;
    bool max_choose_pre, min_choose_pre; //for backtracking, cur at [i][j], true: [i-1][j-1], false: [i-1][j]
    bool max_from_max, min_from_max; //for backtracking, true if the value is obtained from max of the prev cell (determined by above two bools)
} DPCell;
//extern DPCell dp[n][k];

/*
  The set of random tables determining a total order on all bmers.
*/
typedef struct{
    double A;
    bool B1, B2; //true: +1, false: -1
} RandTableCell;
//extern RandTableCell table[k][ALPHABETSIZE];

static inline size_t access2d(const size_t row_len,
			      const size_t i, const size_t j){
    return row_len * i + j;
}

/*
  Initialize the random tables:
  --tp is RandTableCell[k][ALPHABETSIZE] flattened;
  --A[k][ALPHABETSIZE] take double values between 2^{30} and 2^{31};
  --B1[k][ALPHABETSIZE] and B2[k][ALPHABETSIZE] take values +/-1
  such that (B1[i][], B2[i][]) is a permutation of {(+1,+1), 
  (+1,-1), (-1,+1), (-1,-1)}.
 */
void initRandTable(const int k, RandTableCell* tp);

/*
  Save/load a set of random tables to file.
  --tp is RandTableCell[b][ALPHABETSIZE] flattened.
*/
void saveRandTable(const char* filename, const int k, const RandTableCell* tp);
void loadRandTable(const char* filename, const int k, RandTableCell* tp);
void printRandTable(const int k, const RandTableCell* tp);

/*
  Given a k-mer seed, calculate its score according to the random tables tp.
  --tp is RandTableCell[k][ALPHABETSIZE] flattened;
*/
double getSeedScore(const kmer seed, const int k, const RandTableCell* tp);

/*
  Given a char-representation of an n-mer s and a set of random tables,
  fill the dp table according to the total order defined by the random tables.
  --tp is RandTableCell[k][ALPHABETSIZE] flattened;
  --dpp is DPCell[n+1][k+1] flattened.
*/
void fillDPTable(const char* s, const int n, const int k,
		 const RandTableCell* tp, DPCell* dpp);

/*
  After filling the dp table by fillDPTable, backtrack from the given
  cell [n][k] to obtain the selected k-mer.
  Return true if the first char of the n-mer s is selected, false otherwise. 
  This is used for the heuristic speedup of calculating one additional
  column of the dp table.
  --dpp is at least DPCell[n+1][k+1] (maybe DPCell[n+2][k+1]) flattened.
*/
bool backtrackDPTable(const char* s, const int n, const int k,
		      const DPCell* dpp, kmer* result);

/*
  Same as above but also record position info for each char of the 
  selected k-mer, offset by st which is the starting index of the 
  first char of s.
  --pos is assumed to have length at least k
*/
bool backtrackDPTableWithPos(const char* s, const int n, const int k,
			     const DPCell* dpp, kmer* result,
			     const int st, int* pos);


/*
  Given a char-representation of a read and a set of random tables, 
  calculate and store its minSubseq seeds and their starting positions in a
  vector sorted by the positions:
  - only the first one is kept if consecutive windows produce the same seed
    (in the current implementation, "consecutive" is understood to not
    consider windows that do not produce a seed);
  - if a seed appears multiple times in different (non-consecutive) locations,
    it will also be stored multiple times in the vector. 
  Only seeds whose scores are at least the given threshold are generated/stored.
*/
void getSubseqSeedsThreshold(const std::string &read,
			     const int n, const int k,
			     const RandTableCell* tp, const double threshold,
			     std::vector<Seed>& seeds_list);

/*
  Save seeds of a read to file.
  The seeds are saved in ascending order with respect to their starting
  positions in the read. A seed can appear multiple times if it is selected
  from different (non-consecutive) positions in the read.
  Each seed is stored in binary format by fwrite, and can be read by fread.
*/
void saveSubseqSeeds(const char* filename,
		     const std::vector<Seed>& seeds_list);

/*
  Read seeds of a read from file, merge them into a map where the seed
  is key and value is a vector of read ids.
  The seed location info is currently unused.
  ********
  This method is assumed to be called in ascending order of read_id so
  that each vector (associated to a seed) is in sorted order.
*/
void loadSubseqSeeds(const char* filename, const int read_id,
		     std::map<kmer, std::vector<int> > &all_seeds);


/*
  An upper diagonal matrix without the main diagonal,
  valid indices for Table::access(i, j) are 1<=i<j<=n
*/
class Table{
    size_t n;
    unsigned int* arr;

public:
    Table(size_t n);
    ~Table();
    unsigned int& access(size_t i, size_t j);
    void saveNoneZeroEntries(const char* filename);
};

#endif // util.h

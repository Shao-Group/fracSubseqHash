#include "util.h"

kmer encode(const char* s, const int k){
    kmer enc = 0lu;
    int i;
    for(i=0; i<k; i+=1){
	enc = (enc << 2) | alphabetIndex(s[i]);
    }
    return enc;
}

char* decode(const kmer enc, const int k, char* str){
    if(str == NULL){
	str = (char*)malloc(sizeof *str *k);
    }
    kmer enc_copy = enc;
    int i;
    for(i=k-1; i>=0; i-=1){
	str[i] = ALPHABET[enc_copy & 3];
	enc_copy >>= 2;
    }
    return str;
}


void initRandTable(const int k, RandTableCell* tp){
    std::random_device rd;
    std::vector<int> pos;
    std::vector<int> possign;

    int i, j, q;
    for(i=0; i<ALPHABETSIZE; i+=1){
	possign.push_back(i);
    }

    unsigned seed;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution((int64_t)1<<30, (int64_t)1<<31);

    for(i=0; i<k; i+=1)
    {
	for(j=0; j<ALPHABETSIZE; j+=1){
	    tp[access2d(ALPHABETSIZE, i, j)].A = distribution(generator);
	}        


	seed = unsigned(std::chrono::system_clock::now().time_since_epoch().count());
	std::shuffle(possign.begin(), possign.end(), std::default_random_engine(seed));
	for(j=0; j<ALPHABETSIZE; j+=1)
	{
	    q = access2d(ALPHABETSIZE, i, j);
	    tp[q].B1 = possign[j] % 2;
	    tp[q].B2 = possign[j] / 2;
	}
    }
}

void saveRandTable(const char* filename, const int k, const RandTableCell* tp){
    FILE* fout = fopen(filename, "wb");
    fwrite(tp, sizeof(RandTableCell), k * ALPHABETSIZE, fout);
    fclose(fout);
}

void loadRandTable(const char* filename, const int k, RandTableCell* tp){
    FILE* fin = fopen(filename, "rb");
    size_t expected = k * ALPHABETSIZE;
    size_t ret = fread(tp, sizeof(RandTableCell), expected, fin);
    if(ret != expected){
	if(feof(fin)) fprintf(stderr, "Rand tables in %s are too small\n", filename);
	else if(ferror(fin)) fprintf(stderr, "Error reading %s\n", filename);
    }
    fclose(fin);
}

void printRandTable(const int k, const RandTableCell* tp){
    int i, j, q;
    for(i=0; i<k; i+=1){
	for(j=0; j<ALPHABETSIZE; j+=1){
	    printf("%.6lf ", tp[access2d(ALPHABETSIZE, i, j)].A);
	}
	printf("\n");
    }
    printf("\n");

    for(i=0; i<k; i+=1){
	for(j=0; j<ALPHABETSIZE; j+=1){
	    q = access2d(ALPHABETSIZE, i, j);
	    printf("%+d,%+d ", tp[q].B1?1:-1, tp[q].B2?1:-1);
	}
	printf("\n");
    }
    printf("\n");
}


double getSeedScore(const kmer seed, const int k, const RandTableCell* tp){
    double omega = 0;

    int i, cur, q;
    for(i=0; i<k; ++i){
	cur = (seed >> ((k-i-1)<<1)) & 3;
	q = access2d(ALPHABETSIZE, i, cur);
	omega = (tp[q].B1 ? omega : -omega) + (tp[q].B2 ? tp[q].A : -tp[q].A);
    }

    return fabs(omega);
}

void fillDPTable(const char* s, const int n, const int k,
		 const RandTableCell* tp, DPCell* dpp){
    int del = n-k, i, j, c, q;

    memset(dpp, 0, sizeof *dpp * (n+1) * (k+1));

    c = alphabetIndex(s[0]);
    q = access2d(k+1, 1, 1);
    dpp[q].min = dpp[q].max = tp[c].B2 ? tp[c].A : -tp[c].A; //tp[0][c]
    dpp[q].min_choose_pre = dpp[q].max_choose_pre = true;

    int minj, maxj, prev;
    double v1, v2;
    for(i=2; i<=n; ++i){
	minj = std::max(1, i-del);
	maxj = std::min(i, k);
	
	for(j=minj, q=access2d(k+1, i, minj), prev=access2d(k+1, i-1, minj-1);
	    j<=maxj; ++j, ++q, ++prev){
	    //dpp[i][j] = dpp[i-1][j]
	    if(i-1 < j){//[i-1][j] is not a meaningful cell
		dpp[q].min = 1e15;
		dpp[q].max = -1e15;
	    }else{
		memcpy(&dpp[q], &dpp[prev+1], sizeof *dpp);
		dpp[q].min_choose_pre = dpp[q].max_choose_pre = false;
		dpp[q].max_from_max = true;
		dpp[q].min_from_max = false;
	    }
	    
	    //compare with dpp[i-1][j-1]
	    c = access2d(ALPHABETSIZE, j-1, alphabetIndex(s[i-1]));
	    if(tp[c].B1){
		v1 = dpp[prev].min;
		v2 = dpp[prev].max;
	    }else{
		v1 = -dpp[prev].min;
		v2 = -dpp[prev].max;
	    }

	    if(tp[c].B2){
		v1 += tp[c].A;
		v2 += tp[c].A;
	    }else{
		v1 -= tp[c].A;
		v2 -= tp[c].A;
	    }

	    if(v1 < v2){
		if(v1 <= dpp[q].min){
		    dpp[q].min = v1;
		    dpp[q].min_choose_pre = true;
		    dpp[q].min_from_max = false;
		}
		if(v2 >= dpp[q].max){
		    dpp[q].max = v2;
		    dpp[q].max_choose_pre = true;
		    dpp[q].max_from_max = true;
		}
	    }else{
		if(v2 <= dpp[q].min){
		    dpp[q].min = v2;
		    dpp[q].min_choose_pre = true;
		    dpp[q].min_from_max = true;
		}
		if(v1 >= dpp[q].max){
		    dpp[q].max = v1;
		    dpp[q].max_choose_pre = true;
		    dpp[q].max_from_max = false;
		}
	    }
	}//end for j
    }//end for i    
}//end fillDPTable

bool backtrackDPTable(const char* s, const int n, const int k,
		      const DPCell* dpp, kmer* result){
    *result = 0;
    kmer c;
    int i = 0, cur = n;
    bool select, from_max;
    
    int q = access2d(k+1, n, k);
    double score = fabs(dpp[q].min);
    if(dpp[q].max > score){
	score = dpp[q].max;
	select = dpp[q].max_choose_pre;
	from_max = dpp[q].max_from_max;
    }else{
	select = dpp[q].min_choose_pre;
	from_max = dpp[q].min_from_max;
    }

    while(i < (k<<1)){
	if(select){
	    c = alphabetIndex(s[cur-1]);
	    c <<= i;
	    *result |= c;
	    i += 2;
	    q -= (k+2); //[i][j] to [i-1][j-1]
	}else{
	    q -= (k+1); //[i][j] to [i-1][j]
	}
	cur -= 1;
	
	if(from_max){
	    select = dpp[q].max_choose_pre;
	    from_max = dpp[q].max_from_max;
	}else{
	    select = dpp[q].min_choose_pre;
	    from_max = dpp[q].min_from_max;
	}
    }

    return (q == 0);
}

bool backtrackDPTableWithPos(const char* s, const int n, const int k,
			     const DPCell* dpp, kmer* result,
			     const int st, int* pos){
    *result = 0;
    kmer c;
    int i = 0, cur = n;
    bool select, from_max;
    
    int q = access2d(k+1, n, k);
    double score = fabs(dpp[q].min);
    if(dpp[q].max > score){
	score = dpp[q].max;
	select = dpp[q].max_choose_pre;
	from_max = dpp[q].max_from_max;
    }else{
	select = dpp[q].min_choose_pre;
	from_max = dpp[q].min_from_max;
    }

    while(i < (k<<1)){
	if(select){
	    c = alphabetIndex(s[cur-1]);
	    pos[k-1-(i>>1)] = st + cur - 1;
	    c <<= i;
	    *result |= c;
	    i += 2;
	    q -= (k+2); //[i][j] to [i-1][j-1]
	}else{
	    q -= (k+1); //[i][j] to [i-1][j]
	}
	cur -= 1;
	
	if(from_max){
	    select = dpp[q].max_choose_pre;
	    from_max = dpp[q].max_from_max;
	}else{
	    select = dpp[q].min_choose_pre;
	    from_max = dpp[q].min_from_max;
	}
    }

    return (q == 0);
}

static inline void storeSeedWithPosInVector(const kmer seed, const unsigned int pos,
					    std::vector<Seed>& seeds_list){
    //skip the same seed from consecutive positions
    if(seeds_list.size() > 0){
	Seed& s = seeds_list.back();
	if(s.v == seed){
	    ++ s.span;
	    return;
	}
    }
    seeds_list.emplace_back(seed, pos);
}

static inline double getScoreFromDPTable(const int n, const int k,
					 const DPCell* dp){
    int q = access2d(k+1, n, k);
    double score = fabs(dp[q].min);
    if(score < dp[q].max) return dp[q].max;
    else return score;
}

void getSubseqSeedsThreshold(const std::string &read,
			     const int n, const int k,
			     const RandTableCell* tp, const double threshold,
			     std::vector<Seed>& seeds_list){
    size_t len = read.length();
    unsigned int i;
    char cur[n+1];
    kmer seed;
    double score = threshold;
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
	    storeSeedWithPosInVector(seed, i, seeds_list);
	}

	score = getScoreFromDPTable(n+1, k, dp);
	if(score >= threshold){
	    if(!backtrackDPTable(cur, n+1, k, dp, &seed)){//first char not used
		++i; //skip recalculation of next position
		storeSeedWithPosInVector(seed, i, seeds_list);
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
	    storeSeedWithPosInVector(seed, i, seeds_list);
	}
    }
}

void saveSubseqSeeds(const char* filename,
		     const std::vector<Seed>& seeds_list){
    FILE* fout = fopen(filename, "wb");
    for(const Seed& s : seeds_list){
	fwrite(&s, sizeof(s), 1, fout);
    }
    fclose(fout);
}

void loadSubseqSeeds(const char* filename, const int read_id,
		     std::map<kmer, std::vector<int> > &all_seeds){
    FILE* fin = fopen(filename, "rb");
    size_t ret = 1;
    kmer seed;
    int ct;
    while(ret == 1){
	ret = fread(&seed, sizeof(kmer), 1, fin);
	fread(&ct, sizeof(int), 1, fin);
	fseek(fin, sizeof(int)*ct, SEEK_CUR); //skip the position info for now
	auto result = all_seeds.insert(std::pair<kmer, std::vector<int> >(seed, std::vector<int>(1, read_id)));
	if(result.second == false && result.first->second.back() < read_id){
	    result.first->second.push_back(read_id);
	}
    }
    if(ferror(fin)){
	fprintf(stderr, "Error reading %s\n", filename);
    }
    fclose(fin);
}

Table::Table(size_t n): n(n){
    size_t size = (n*(n-1))>>1;
    arr = new unsigned int[size];
    memset(arr, 0, sizeof *arr * size);
}

Table::~Table(){
    delete[] arr;
}

unsigned int& Table::access(size_t i, size_t j){
    return arr[((((n<<1)-i)*(i-1))>>1)+j-i-1];
}

void Table::saveNoneZeroEntries(const char* filename){
    FILE* fout = fopen(filename, "w");

    size_t i, j, k, size=(n*(n-1))>>1;
    for(k=0, i=1, j=2; k<size; ++k){
	if(arr[k] > 0){
	    fprintf(fout, "%zu %zu %u\n", i, j, arr[k]);
	}
	j += 1;
	if(j > n){
	    i += 1;
	    j = i + 1;
	}
    }
    
    fclose(fout);
}

/*
  Link all seeds of a read as a path. Multiple paths are merged into a graph.
  Only keep seeds that appear on multiple reads (paths).

  By: Ke@PSU
  Last edited: 03/12/2023
*/

#ifndef _SEEDSGRAPH_H
#define _SEEDSGRAPH_H 1

#include <map>
#include <string>
//#include <mutex>
#include <fstream>

template<class T>
class SeedsGraph{
public:
    struct Locus;
    struct Path;
    class Node;
    struct ReadPath;
    
private:
    std::map<T, Node> nodes;
    // only addNode is protected as other operations are not to be used
    // in parallel
    //std::mutex under_construction;
    std::vector<ReadPath> paths;

    /*
      Helper function for removeNode() and removeUniqSeeds().
    */
    void skipNode(Node* n);

    /*
      Helper functions for saveGraphToDot().
      Print all the nodes/edges/paths in the graph in the dot format to fout.
    */
    template<class... Args>
    void printNodesInDot(std::ofstream& fout,
			 std::string (*decode)(const T&, Args...),
			 Args... args) const;
    void printEdgesInDot(std::ofstream& fout) const;
    void printReadPathsInDot(std::ofstream& fout) const;
    
public:
    SeedsGraph() {};
    SeedsGraph(const int num_read_paths){
	paths.reserve(num_read_paths);
    };
    
    /*
      Getters.
    */
    size_t numNodes() const;
    Node* getNode(const T& key) const;

    /*
      Add a node for a given key into the graph, do nothing if such 
      a node already exists. In either case, return a pointer to 
      the node corresponding to the given key.
    */
    Node* addNode(T& key);

    /*
      Add a read represented by a head and a tail pointers to some nodes
      in this graph.
    */
    ReadPath& addReadPath(const size_t read_id,
			  Node* head, Node* tail);

    /*
      Remove a node n, n is assumed to be in nodes.
      All the adjacent in- and out-edges are sewed accordingly.
    */
    void removeNode(Node* n);

    /*
      Remove all nodes (seeds) that only appear in one read.
     */
    void removeUniqSeeds();

    /*
      Graph IO with the given filename.
      The decode function is used for transforming the key into a string.
    */
    template<class... Args>
    void saveGraphToDot(const char* filename,
			std::string (*decode)(const T&, Args...),
			Args... args) const; //to dot format
    void saveGraph(const char* filename) const; //to binary
    void loadGraph(const char* filename); //from binary
};



template<class T>
struct SeedsGraph<T>::Locus{//location info (on the read it originates from) of a seed
    size_t read_id;
    size_t pos;

    Locus(const size_t id, const size_t pos): read_id(id), pos(pos) {};
    Locus(): Locus(0, 0) {};
    
    /*
      Ordered first by read_id, then by pos.
    */
    bool operator < (const Locus& x) const;
    /*
      Equals if have the same read_idx and pos.
    */
    bool operator == (const Locus& x) const;
};



template<class T>
struct SeedsGraph<T>::Path{
    Node *prev, *next;
    Path(Node* p, Node* n): prev(p), next(n) {};
    Path(const Path& o) = delete;
};



template<class T>
class SeedsGraph<T>::Node{
    // only protect addPrev/addNext as other operations are not to
    // be used in parallel
    //std::mutex under_construction;
    
public:
    T seed;
    size_t id; //assigned in construction order
    std::map<Locus, Path> locations;
    size_t read_ct; // number of distinct reads that contain this seed
    
    Node(T& seed, size_t id):seed(std::move(seed)), id(id), read_ct(0) {};
    Node(Node&& other);
    /*
      Add an edge to this node; 
      Increment read_ct if the read_id does not appear in locations.
    */
    void addPrev(const size_t read_id, const size_t pos, Node* prev);
    void addNext(const size_t read_id, const size_t pos, Node* next);

    /*
      String representation of this node (without edge info) in dot format.
    */
    template<class... Args>
    std::string toString(std::string (*decode)(const T&, Args...),
			 Args... args) const;
    
    /*
      Node IO to binary file, helper functions to SeedsGraph::saveGraph.
    */
    void saveNode(std::ofstream& fout) const;
    void loadNode(std::ifstream& fin);
    
};

template<class T>
struct SeedsGraph<T>::ReadPath{
    size_t read_idx;
    Node *head, *tail;


    ReadPath(const size_t read_idx, Node* head, Node* tail):
	read_idx(read_idx), head(head), tail(tail) {};
    
    ReadPath(const size_t read_idx): ReadPath(read_idx, nullptr, nullptr) {};

    ReadPath(ReadPath&& o): read_idx(std::exchange(o.read_idx, 0)),
			    head(std::exchange(o.head, nullptr)),
			    tail(std::exchange(o.tail, nullptr)) {};

    /*
      ReadPath IO to binary file, helper functions to SeedsGraph::saveGraph.
    */
    void saveReadPath(std::ofstream& fout) const;
    void loadReadPath(std::ifstream& fin,
		      const std::map<size_t, Node*>& dict); 
};











/************* LOCUS *************/

template<class T>
bool SeedsGraph<T>::Locus::operator < (const Locus& x) const{
    if(read_id == x.read_id) return pos < x.pos;
    else return read_id < x.read_id;
}

template<class T>
bool SeedsGraph<T>::Locus::operator == (const Locus& x) const{
    return (read_id == x.read_id) && (pos == x.pos);
}


/************* NODE *************/

template<class T>
SeedsGraph<T>::Node::Node(Node&& other): seed(std::move(other.seed)),
					 id(std::exchange(other.id, 0)) {
    //const std::lock_guard<std::mutex> lock(other.under_construction);
    locations = std::move(other.locations);
    read_ct = std::exchange(other.read_ct, 0);
}

template<class T>
void SeedsGraph<T>::Node::addPrev(const size_t read_id, const size_t pos,
				  Node* prev){
    //const std::lock_guard<std::mutex> lock(under_construction);
    
    auto result = locations.emplace(std::piecewise_construct,
				    std::forward_as_tuple(read_id, pos),
				    std::forward_as_tuple(prev, nullptr));
    //new Locus inserted, since each read is processed in ascending pos
    //if this seed has been found on this same read before, it must
    //appear in the previous position in the map. Therefore, it is
    //sufficient to check the previous position (by decrementing the
    //iterator). If the read id is different from current id, then it's
    //the first appearance of this seed on this read.
    if(result.second) {
	if(result.first == locations.begin()) ++ read_ct;
	else{
	    -- result.first;
	    if(result.first->first.read_id < read_id) ++ read_ct;
	}
    }
    else result.first->second.prev = prev;
}

template<class T>
void SeedsGraph<T>::Node::addNext(const size_t read_id, const size_t pos,
				  Node* next){
    //const std::lock_guard<std::mutex> lock(under_construction);
    
    auto result = locations.emplace(std::piecewise_construct,
				    std::forward_as_tuple(read_id, pos),
				    std::forward_as_tuple(nullptr, next));
    if(result.second) {
	if(result.first == locations.begin()) ++ read_ct;
	else{
	    -- result.first;
	    if(result.first->first.read_id < read_id) ++ read_ct;
	}
    }
    else result.first->second.next = next;				
}


template<class T> template<class... Args>
std::string SeedsGraph<T>::Node::toString(
    std::string (*decode)(const T&, Args...),
    Args... args) const{

    return "n" + std::to_string(id) + " [label=\"" +
	decode(seed, args...) + "\"];";
}

template<class T>
void SeedsGraph<T>::Node::saveNode(std::ofstream& fout) const{
    fout.write(reinterpret_cast<const char*>(&seed), sizeof(seed));
    fout.write(reinterpret_cast<const char*>(&id), sizeof(id));
    fout.write(reinterpret_cast<const char*>(&read_ct), sizeof(read_ct));
    size_t num_locations = locations.size();
    fout.write(reinterpret_cast<const char*>(&num_locations), sizeof(num_locations));
    Node* x;
    decltype(id) nullid = 0;
    for(const auto& l : locations){
	fout.write(reinterpret_cast<const char*>(&(l.first)), sizeof(l.first));
	x = l.second.prev;
	fout.write(reinterpret_cast<const char*>(&(x?x->id:nullid)), sizeof(nullid));
	x = l.second.next;
	fout.write(reinterpret_cast<const char*>(&(x?x->id:nullid)), sizeof(nullid));
    }
}

template<class T>
void SeedsGraph<T>::Node::loadNode(std::ifstream& fin){
    //seed and id have been read to create this node
    //fin.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    //fin.read(reinterpret_cast<char*>(&id), sizeof(id));
    fin.read(reinterpret_cast<char*>(&read_ct), sizeof(read_ct));
    size_t i, num_locations;
    fin.read(reinterpret_cast<char*>(&num_locations), sizeof(num_locations));
    
    Locus l;
    size_t x, y;
    for(i=0; i<num_locations; ++i){
	fin.read(reinterpret_cast<char*>(&l), sizeof(l));
	fin.read(reinterpret_cast<char*>(&x), sizeof(x));
	fin.read(reinterpret_cast<char*>(&y), sizeof(y));
	locations.emplace_hint(locations.end(),
			       std::piecewise_construct,
			       std::forward_as_tuple(l.read_id, l.pos),
			       std::forward_as_tuple(reinterpret_cast<Node*>(x),
						     reinterpret_cast<Node*>(y)));
    }
}


/************* READPATH *************/
template<class T>
void SeedsGraph<T>::ReadPath::saveReadPath(std::ofstream& fout) const{
    if(head){//skip empty ReadPaths
	fout.write(reinterpret_cast<const char*>(&read_idx), sizeof(read_idx));
	fout.write(reinterpret_cast<const char*>(&(head->id)), sizeof(head->id));
	fout.write(reinterpret_cast<const char*>(&(tail->id)), sizeof(tail->id));
    }
}

template<class T>
void SeedsGraph<T>::ReadPath::loadReadPath(std::ifstream& fin,
					   const std::map<size_t, Node*>& dict){
    //id has been read for creating this ReadPath
    //fin.read(reinterpret_cast<char*>(&read_idx), sizeof(read_idx));
    size_t id;
    fin.read(reinterpret_cast<char*>(&id), sizeof(id));
    auto it = dict.find(id);
    //assert(it != dict.end())
    head = it->second;
    
    fin.read(reinterpret_cast<char*>(&id), sizeof(id));
    it = dict.find(id);
    //assert(it != dict.end())
    tail = it->second;
}


/************* SEEDSGRAPH *************/

template<class T>
size_t SeedsGraph<T>::numNodes() const{
    return nodes.size();
}

template<class T>
typename SeedsGraph<T>::Node* SeedsGraph<T>::getNode(const T& key) const{
    auto it = nodes.find(key);
    if(it != nodes.end()) return &(it->second);
    else return nullptr;
}

template<class T>
typename SeedsGraph<T>::Node* SeedsGraph<T>::addNode(T& key){
    //const std::lock_guard<std::mutex> lock(under_construction);

    auto it = nodes.lower_bound(key);
    if(it != nodes.end() && it->first == key) return &(it->second);
    else{
	size_t id = nodes.size() + 1;
	it = nodes.emplace_hint(it, key, Node(key, id));
	return &(it->second);
    }
}

template<class T>
inline typename SeedsGraph<T>::ReadPath&
SeedsGraph<T>::addReadPath(const size_t read_id,
			   Node* head, Node* tail){
    paths.emplace_back(read_id, head, tail);
    return paths.back();
}

template<class T>
void SeedsGraph<T>::skipNode(Node* n){
    //sewing all the paths through this node to skip this node
    for(auto& l : n->locations){
	Path& p = l.second;
	if(p.prev){
	    //search by current locus, should find the locus
	    //of the prev seed on this read
	    auto predecessor = p.prev->locations.lower_bound(l.first);
	    if(predecessor != p.prev->locations.begin()){//should always be true
		--predecessor;
		//assert(predecessor->second.next == n)
		predecessor->second.next = p.next;
	    }
	}
	if(p.next){
	    //search by current locus, should find the locus
	    //of the next seed on this read
	    auto successor = p.next->locations.upper_bound(l.first);
	    //assert(successor->second.prev == n)
	    successor->second.prev = p.prev;
	}
    }
}

template<class T>
void SeedsGraph<T>::removeNode(Node* n){
    skipNode(n);
    nodes.erase(n->seed);
}

template<class T>
void SeedsGraph<T>::removeUniqSeeds(){
    // move the head and tail pointers of each path to point to the first
    // and last (resp.) non-unique seeds in the path
    Node* cur;
    for(ReadPath& p : paths){
	cur = p.head;
	while(cur && cur->read_ct < 2){
	    //if read_ct is 1, the first (or the only, if there is no loop in the path)
	    //in locations must correspond to the cur path
	    p.head = cur->locations.begin()->second.next;
	    removeNode(cur);
	    cur = p.head;
	}

	if(cur){
	    cur = p.tail;
	    while(cur && cur->read_ct < 2){
		p.tail = cur->locations.rbegin()->second.prev;
		removeNode(cur);
		cur = p.tail;
	    }
	}else{
	    //entire path has been removed
	    p.tail = nullptr;
	}
    }

    //remove remaining unique nodes
    auto it = nodes.begin();
    while(it != nodes.end()){
	if(it->second.read_ct < 2){
	    skipNode(&(it->second));
	    it = nodes.erase(it);
	}else{
	    ++ it;
	}
    }
}

template<class T> template<class... Args>
void SeedsGraph<T>::saveGraphToDot(const char* filename,
				   std::string (*decode)(const T&, Args...),
				   Args... args) const{
    std::ofstream fout(filename, std::ofstream::out);
    
    fout << "digraph{" << std::endl;
    
    printNodesInDot(fout, decode, args...);
    printEdgesInDot(fout);
    printReadPathsInDot(fout);    
    
    fout << "} //end of graph" << std::endl;
}

template<class T> template<class... Args>
void SeedsGraph<T>::printNodesInDot(std::ofstream& fout,
				    std::string (*decode)(const T&, Args...),
				    Args... args) const{
    for(const auto& it : nodes){
	fout << it.second.toString(decode, args...) << std::endl;
    }
}

template<class T>
void SeedsGraph<T>::printEdgesInDot(std::ofstream& fout) const{
    for(const auto& it : nodes){
	const Node& cur = it.second;
	std::map<size_t, unsigned int> out_edges; // next_node.id -> weight
	for(const auto& loc_it : cur.locations){
	    if(loc_it.second.next){
		auto result = out_edges.emplace(loc_it.second.next->id, 1);
		if(!result.second){//key exists
		    ++ result.first->second;
		}
	    }
	}

	for(const auto& n : out_edges){
	    fout << "n" << cur.id << " -> n" << n.first
		 << " [label=\"" << n.second << "\"];" << std::endl; 
	}
    }
}

template<class T>
void SeedsGraph<T>::printReadPathsInDot(std::ofstream& fout) const{
    for(const ReadPath& p : paths){
	if(p.head){
	    fout << "st" << p.read_idx << " [label=\"Read " << p.read_idx
		 << " head\"];" << std::endl;
	    fout << "ed" << p.read_idx << " [label=\"Read " << p.read_idx
		 << " tail\"];" << std::endl;
	    //add edge from head to first node
	    fout << "st" << p.read_idx << " -> n" << p.head->id << ";" << std::endl;
	    
	    //add edge from last node to tail
	    fout << "n" << p.tail->id << " -> ed" << p.read_idx << ";" << std::endl;
	}else{
	    fout << "// read " << p.read_idx
		 << " has no overlapping seeds with others" << std::endl;
	}
    }
}


/*
  Save the graph to the given filename.
*/
template<class T>
void SeedsGraph<T>::saveGraph(const char* filename) const{
    std::ofstream fout(filename, std::ios_base::binary);
    size_t num_nodes = nodes.size();
    fout.write(reinterpret_cast<char*>(&num_nodes), sizeof(num_nodes));
    for(const auto& n : nodes){
	n.second.saveNode(fout);
    }
    for(const auto& p : paths){
	p.saveReadPath(fout);
    }
}

template<class T>
void SeedsGraph<T>::loadGraph(const char* filename){
    std::ifstream fin(filename, std::ios_base::binary);
    
    nodes.clear();
    size_t i, num_nodes = nodes.size();
    fin.read(reinterpret_cast<char*>(&num_nodes), sizeof(num_nodes));
    T seed;
    size_t id;

    std::map<size_t, Node*> dict; //restore the pointers by id
    dict.emplace(0, nullptr);
    
    for(i=0; i<num_nodes; ++i){
	fin.read(reinterpret_cast<char*>(&seed), sizeof(seed));
	fin.read(reinterpret_cast<char*>(&id), sizeof(id));
	auto it = nodes.emplace_hint(nodes.end(),
				     std::piecewise_construct,
				     std::forward_as_tuple(seed),
				     std::forward_as_tuple(seed, id));
	it->second.loadNode(fin);
	dict.emplace(id, &(it->second));
    }

    for(auto& it : nodes){
	for(auto& l : it.second.locations){
	    Path& p = l.second;
	    id = reinterpret_cast<size_t>(p.prev);
	    auto result = dict.find(id);
	    //assert(result != dict.end())
	    p.prev = result->second;

	    id = reinterpret_cast<size_t>(p.next);
	    result = dict.find(id);
	    //assert(result != dict.end())
	    p.next = result->second;
	}
    }

    paths.clear();
    while(fin.read(reinterpret_cast<char*>(&id), sizeof(id))
	  && fin.gcount() == sizeof(id)){
	paths.emplace_back(id);
	paths.back().loadReadPath(fin, dict);
    }
}

#endif // SeedsGraph.h

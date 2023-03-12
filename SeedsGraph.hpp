/*
  Link all seeds of a read as a path. Multiple paths are merged into a graph.
  Only keep seeds that appear on multiple reads (paths).

  By: Ke@PSU
  Last edited: 03/03/2023
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
    
private:
    std::map<T, Node> nodes;
    // only addNode is protected as other operations are not to be used
    // in parallel
    //std::mutex under_construction; 
    
public:
    /*
      Getters.
    */
    size_t numNodes() const;
    Node* getNode(const T& key) const;

    /*
      Adds a node for a given key into the graph, do nothing if such 
      a node already exists. In either case, return a pointer to 
      the node corresponding to the given key.
    */
    Node* addNode(T& key);

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
      Print all the nodes/edges in the graph in the dot format to fout.
      The decode function is used for transforming the key into a string.
    */
    template<class... Args>
    void printNodesInDot(std::ofstream& fout,
			 std::string (*decode)(const T&, Args...),
			 Args... args) const;
    void printEdgesInDot(std::ofstream& fout) const;

    /*
      Graph IO with the given filename.
    */
    void saveGraph(const char* filename);
    void loadGraph(const char* filename);
};



template<class T>
struct SeedsGraph<T>::Locus{//location info (on the read it originates from) of a seed
    const size_t read_id;
    const size_t pos;

    Locus(const size_t id, const size_t pos): read_id(id), pos(pos) {};
    
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
	size_t id = nodes.size();
	it = nodes.emplace_hint(it, key, Node(key, id));
	return &(it->second);
    }
}

template<class T>
void SeedsGraph<T>::removeNode(Node* n){
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

    nodes.erase(n->seed);
}

/*
  Remove all nodes (seeds) that only appear in one read.
*/
template<class T>
void SeedsGraph<T>::removeUniqSeeds(){
    for(auto& it : nodes){
	if(it.second.read_ct < 2){
	    removeNode(&it.second);
	}
    }
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
	std::map<Node*, unsigned int> out_edges;
	for(const auto& loc_it : cur.locations){
	    if(loc_it.second.next){
		auto result = out_edges.emplace(loc_it.second.next, 1);
		if(!result.second){//key exists
		    ++ result.first->second;
		}
	    }
	}

	for(const auto& n : out_edges){
	    fout << "n" << cur.id << " -> n" << n.first->id
		 << " [label=\"" << n.second << "\"];" << std::endl; 
	}
    }
}


/*
  Save the graph to the given filename.
*/
template<class T>
void SeedsGraph<T>::saveGraph(const char* filename){
    //TODO
}

template<class T>
void SeedsGraph<T>::loadGraph(const char* filename){
    //TODO
}

#endif // SeedsGraph.h

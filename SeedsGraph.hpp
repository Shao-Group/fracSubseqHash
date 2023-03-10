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
#include <mutex>
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
    std::mutex under_construction; 
    
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
      All the adjacent edges are also removed.
    */
    void removeNode(Node* n);

    /*
      Remove all nodes (seeds) that only appear in one read.
     */
    void removeUniqSeeds();

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
};



template<class T>
class SeedsGraph<T>::Node{
    // only protect addPrev/addNext as other operations are not to
    // be used in parallel
    std::mutex under_construction; 
    
public:
    const T seed;
    std::map<Locus, Path> locations;
    size_t read_ct; // number of distinct reads that contain this seed
    
    Node(T& seed):seed(std::move(seed)), read_ct(1) {};
    Node(Node&& other);
    /*
      Add an edge to this node; 
      Increment read_ct if the read_id does not appear in locations.
    */
    void addPrev(const size_t read_id, const size_t pos, Node* prev);
    void addNext(const size_t read_id, const size_t pos, Node* next);
    
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
SeedsGraph<T>::Node::Node(Node&& other): seed(std::move(other.seed)) {
    const std::lock_guard<std::mutex> lock(other.under_construction);
    locations = std::move(other.locations);
    read_ct = std::exchange(other.read_ct, 0);
}

template<class T>
void SeedsGraph<T>::Node::addPrev(const size_t read_id, const size_t pos,
				  Node* prev){
    const std::lock_guard<std::mutex> lock(under_construction);
    
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
    const std::lock_guard<std::mutex> lock(under_construction);
    
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
    const std::lock_guard<std::mutex> lock(under_construction);

    auto it = nodes.lower_bound(key);
    if(it != nodes.end() && it->first == key) return &(it->second);
    else{
	it = nodes.emplace_hint(it, key, Node(key));
	return &(it->second);
    }
}

/*
  Remove a node n, n is assumed to be in nodes.
  All the adjacent edges are also removed.
*/
template<class T>
void SeedsGraph<T>::removeNode(Node* n){
    //TODO
}

/*
  Remove all nodes (seeds) that only appear in one read.
*/
template<class T>
void SeedsGraph<T>::removeUniqSeeds(){
    //TODO
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

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

template<class T>
class SeedGraph{
public:
    struct Locus;
    struct Edge;
    class Node;
    
private:
    std::map<T, Node> nodes;
    
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
    Node* addNode(const T& key);

    /*
      Remove a node n, n is assumed to be in nodes.
      All the adjacent edges are also removed.
    */
    void removeNode(Node* n);
};

template<class T>
struct SeedGraph<T>::Locus{//location info (on the read it originates from) of an edge
    const size_t read_id;
    //starting positions on the read of the head/tail nodes (an edge is tail -> head)
    const unsigned int head_pos, tail_pos;

    Locus(const size_t id, const unsigned int headp, const unsigned int tailp):
	read_id(id), head_pso(headp), tail_pos(tailp) {};
    
    /*
      Ordered first by read_idx, then by tail_pos. (Equality of these two implies the same head_pos.)
    */
    bool operator < (const Locus& x) const;
    /*
      Equals if have the same read_idx and pos.
    */
    bool operator == (const Locus& x) const;
};

template<class T>
struct SeedGraph<T>::Edge{
    unsigned int weight; //number of distinct reads that contain this edge
    std::set<Locus> reads;
};

template<class T>
class SeedGraph<T>::Node{
    const T seed;
    size_t read_ct; // number of distinct reads that contain this seed
    std::map<Node*, Edge> in, out; //in- and out-edges
    
public:
    Node(const T& seed):key(seed) {};
    Node(Node&& other):key(std::move(other.key)),
		       read_ct(std::exchange(other.read_ct)),
		       in(std::move(other.in)),
		       out(std::move(other.out)) {};

    /*
      Adds a path prev->this-|. 
      If prev is not nullptr, then there should be a path in prev:
      pp->prev-|, which should be updated to pp->prev->this.
    */
    void addPath(const size_t read_idx, const size_t cur_pos,
		 const size_t prev_pos, Node* prev);

    /*
      Given the read_id and pos on read, return the link through this node
      which is a subpath of this read.
      --For the exact version, assume the seed at pos on the read is 
      this->key, otherwise return end().
      --For the lowerbd version, try to search for the smallest pos on read 
      such that the seed is this->key, and pos >= pos_lower_bd; if no
      such pos exists, return end().
    */
    auto getPathExact(const size_t read_idx, const size_t cur_pos) const;
    auto getPathLowerBD(const size_t read_idx,
			const size_t pos_lower_bd) const;
    /*
      Takes a function that transforms the key to a string.
      Returns a string representation of this node, 
      -- if long_fmt, include its id, key and all the in- and out-edges 
      in pairs representing local paths passing through this node.
      -- otherwise, only provide id and key in dot format.
    */
    template<class... Args>
    std::string toString(bool long_fmt,
	std::string (*decode)(const T&, Args...), Args... args) const;
    /*
      Use id2
    */
    template<class... Args>
    std::string toString2(bool long_fmt,
	std::string (*decode)(const T&, Args...), Args... args) const;
    
    /*
      Ordered by key.
    */
    bool operator < (const Node& x) const;
    /*
      Equals if have the same key, ie, each seed is represented by a 
      single node.
    */
    bool operator == (const Node& x) const;
};

#include "SeedsGraph.tpp"

#endif // SeedsGraph.h

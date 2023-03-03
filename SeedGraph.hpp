/*
  Link all seeds of a read as a path. Multiple paths are merged into a graph.

  By: Ke@PSU
  Last edited: 02/25/2023
*/

#ifndef _SEEDGRAPH_H
#define _SEEDGRAPH_H 1

#include <map>
#include <string>
#include <fstream>

template<class T>
class SeedGraph{
public:
    struct Locus;
    struct Link;
    class Node;
    
private:
    std::map<T, Node> nodes;
    
public:
    /*
      Getters.
    */
    size_t numNodes() const;
    Node* getNode(T key) const;

    /*
      Print all the nodes in the graph in the dot format to fout.
      The decode function is used for transforming the key into a string.
    */
    template<class... Args>
    void printNodesInDot(std::ofstream& fout,
	std::string (*decode)(const T&, Args...), Args... args) const;

    /*
      Adds a node for a given key into the graph, do nothing if such 
      a node already exists. In either case, return a pointer to 
      the node corresponding to the given key.
    */
    Node* addNode(T key);
};

template<class T>
struct SeedGraph<T>::Locus{
    size_t read_idx;
    size_t pos;

    Locus(const size_t read_idx, const size_t pos):
	read_idx(read_idx), pos(pos) {};

    /*
      Ordered first by read_idx, then by pos.
    */
    bool operator < (const Locus& x) const;
    /*
      Equals if have the same read_idx and pos.
    */
    bool operator == (const Locus& x) const;
};

template<class T>
struct SeedGraph<T>::Link{
    Node *prev, *next;

    Link(Node* prev):prev(prev), next(nullptr) {};
};

template<class T>
class SeedGraph<T>::Node{
public:
    size_t id; // assigned in construction order
    size_t id2; //assigned in order of output
    T key;
    /*
      Stores prev and next for each path passing through this node,
      indexed by the read_idx of that path and the position on the
      read from where the key stored in this node is obtained (as a seed).
    */
    std::map<Locus, Link> paths;

    Node(const size_t id, const T seed):id(id), id2(0), key(seed) {};
    Node(Node&& other):id(std::exchange(other.id, 0)),
		       id2(std::exchange(other.id2, 0)),
		       key(std::move(other.key)),
		       paths(std::move(other.paths)) {};

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

#include "SeedGraph.tpp"

#endif // SeedGraph.h

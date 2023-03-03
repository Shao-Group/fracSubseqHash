//-*-C-*- for emacs 

/***** LOCUS *****/
template<class T>
bool SeedGraph<T>::Locus::operator < (const Locus& x) const{
    if(read_idx == x.read_idx) return pos < x.pos;
    else return read_idx < x.read_idx;
}

template<class T>
bool SeedGraph<T>::Locus::operator == (const Locus& x) const{
    return (read_idx == x.read_idx) && (pos == x.pos);
}

/***** NODE *****/
template<class T>
bool SeedGraph<T>::Node::operator < (const Node& x) const{
    return key < x.key;
}

template<class T>
bool SeedGraph<T>::Node::operator == (const Node& x) const{
    return key == x.key;
}

template<class T>
void SeedGraph<T>::Node::addPath(const size_t read_idx, const size_t cur_pos,
				 const size_t prev_pos, Node* prev){
    this->paths.emplace(std::piecewise_construct,
			std::forward_as_tuple(read_idx, cur_pos),
			std::forward_as_tuple(prev));
    if(prev != nullptr){
	auto it = prev->paths.find(Locus(read_idx, prev_pos));
	it->second.next = this;
    }
}


template<class T>
auto SeedGraph<T>::Node::getPathExact(const size_t read_idx,
				      const size_t cur_pos) const{
    return paths.find(Locus(read_idx, cur_pos));
}

template<class T>
auto SeedGraph<T>::Node::getPathLowerBD(const size_t read_idx,
					const size_t pos_lower_bd) const{
    auto it = paths.lower_bound(Locus(read_idx, pos_lower_bd));
    if(it->first.read_idx == read_idx) return it;
    else return paths.end();
}


template<class T> template<class... Args>
std::string SeedGraph<T>::Node::toString(bool long_fmt,
    std::string (*decode)(const T&, Args...), Args... args) const{
    
    std::string result = "n" + std::to_string(id)
    + " [label=\"" + decode(key, args...) + "\"];\n";
    if(long_fmt){
	for(const auto& p : paths){
	    result += "read " + std::to_string(p.first.read_idx)
		+ ", pos " + std::to_string(p.first.pos) + "\nprev: "
		+ (p.second.prev? decode(p.second.prev->key, args...) : "--")
		+ "\nnext: "
		+ (p.second.next? decode(p.second.next->key, args...) : "--")
		+ "\n";
	}
    }
    return result;
}


template<class T> template<class... Args>
std::string SeedGraph<T>::Node::toString2(bool long_fmt,
    std::string (*decode)(const T&, Args...), Args... args) const{
    
    std::string result = "n" + std::to_string(id2)
    + " [label=\"" + decode(key, args...) + "\"];\n";
    if(long_fmt){
	for(const auto& p : paths){
	    result += "read " + std::to_string(p.first.read_idx)
		+ ", pos " + std::to_string(p.first.pos) + "\nprev: "
		+ (p.second.prev? decode(p.second.prev->key, args...) : "--")
		+ "\nnext: "
		+ (p.second.next? decode(p.second.next->key, args...) : "--")
		+ "\n";
	}
    }
    return result;
}

/***** SEEDGRAPH *****/
template<class T>
size_t SeedGraph<T>::numNodes() const{
    return nodes.size();
}

template<class T>
typename SeedGraph<T>::Node* SeedGraph<T>::getNode(T key) const{
    auto it = nodes.find(key);
    if(it != nodes.end()) return &(it->second);
    else return nullptr;
}

template<class T>
typename SeedGraph<T>::Node* SeedGraph<T>::addNode(T key){
    auto it = nodes.lower_bound(key);
    if(it == nodes.end() || it->first != key){//key does not exist in map
	size_t id = nodes.size();
	it = nodes.emplace_hint(it, key, Node(id, key));
    }
    return &(it->second);
}



template<class T> template<class... Args>
void SeedGraph<T>::printNodesInDot(std::ofstream& fout,
  std::string (*decode)(const T&, Args...), Args... args) const{
    for(const auto& it : nodes){
	fout << it.second.toString(false, decode, args...);
    }
}

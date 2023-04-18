#!/usr/bin/env python3

import sys
import os.path

dir = "sample-reads"
header_ext = "header-sorted"
found_ext = "all-pair"

argc = len(sys.argv)
if argc < 3 or argc > 4:
   example = "sim-e-coli-pb-le20k-nosub"
   print("Usage: makeStringGraph.py <sample-basename> <overlap-basename> [overlap-extension]")
   print(f'Example: makeStringGraph.py {example} {example}-n60-k42-t0/unionPos 10')
   print(f'Will use files {dir}/{example}.{header_ext} and {dir}/{example}-n60-k42-t0/unionPos.{found_ext}10')
   exit(1)


output_filename = f'{dir}/{sys.argv[2]}.graphml'
if argc == 4:
   found_ext += sys.argv[3]
   output_filename += sys.argv[3]

#check all needed files exist
header_name = f'{dir}/{sys.argv[1]}.{header_ext}'
if not os.path.isfile(header_name):
   print(f'Cannot read node file: {header_name}')
   exit(1)

found_name = f'{dir}/{sys.argv[2]}.{found_ext}'
if not os.path.isfile(found_name):
   print(f'Cannot read found all file: {found_name}')
   exit(1)



#generate graph-tool graph
import itertools
from graph_tool.all import *
g = Graph()
g.set_fast_edge_removal(True) #add edge becomes O(1) with O(E) additional space

#create nodes from the sorted header file
try:
   with open(header_name) as f:
      vmap = {int(read):int(idx) for read,idx in
              zip([line.split()[0][1:] for line in f.readlines()],
                  itertools.count(0,1))}         
except OSError as e:
   print(e.strerror)
   exit(1)

g.add_vertex(len(vmap))
g.vertex_properties["id"] = g.new_vertex_property("int")
# dict preserves input order since 3.7, otherwise use
# g.vp.id.a = list(map(lambda kv: kv[0], sorted(vmap.items(), key=lambda kv: kv[1])))
g.vp.id.a = list(vmap.keys())

#create all edges from the found_all file
try:
   with open(found_name) as f:
      g.add_edge_list(list(map(lambda x:(vmap[int(x[0])],
                                         vmap[int(x[1])],
                                         int(x[2])),
                               (line.split() for line in f.readlines()))),
                      eprops=[("weight", "int")])
except OSError as e:
   print(e.strerror)
   exit(1)

g.save(output_filename)

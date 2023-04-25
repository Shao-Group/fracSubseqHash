#!/usr/bin/env python3

import sys
import os.path
import itertools
import numpy as np
from graph_tool.all import *


def main(argc, argv):
   dir = "sample-reads"
   header_ext = "header-sorted"
   irr_ext = "irr-edges-directed"
   true_ext = "truepairs-directed"
   found_ext = "all-pair"
   
   if argc < 3 or argc > 4:
      example = "sim-e-coli-pb-le20k-nosub"
      print("Usage: makeGraphFromOverlap.py <sample-basename> <overlap-basename> [overlap-extension]")
      print(f'Example: makeGraphFromOverlap.py {example} {example}-n60-k42-t0/unionPos 10')
      print(f'Will use files {dir}/{example}.{{{header_ext},{true_ext},{irr_ext}}} and {dir}/{example}-n60-k42-t0/unionPos.{found_ext}10')
      exit(1)


   output_filename = f'{dir}/{argv[2]}.graphml'
   if argc == 4:
      found_ext += argv[3]
      output_filename += argv[3]

   #check all needed files exist
   header_name = f'{dir}/{argv[1]}.{header_ext}'
   if not os.path.isfile(header_name):
      print(f'Cannot read node file: {header_name}')
      exit(1)

   true_name = f'{dir}/{argv[1]}.{true_ext}'
   if not os.path.isfile(true_name):
      print(f'Cannot read true pair file: {true_name}')
      exit(1)

   irr_name = f'{dir}/{argv[1]}.{irr_ext}'
   if not os.path.isfile(irr_name):
      print(f'Cannot read irr pair file: {irr_name}')
      exit(1)

   found_name = f'{dir}/{argv[2]}.{found_ext}'
   if not os.path.isfile(found_name):
      print(f'Cannot read found all file: {found_name}')
      exit(1)


   #generate graph-tool graph
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


   # 0 - wrong, 1 - correct, 2 - irreducible
   g.ep["type"] = g.new_edge_property("int8_t")
   try:
      with open(true_name) as f:
         for line in f.readlines():
            l = line.split()
            s = vmap[int(l[0])]
            t = vmap[int(l[1])]
            e = g.edge(s, t)
            if e is None:
               e = g.add_edge(s, t)
            g.ep.type[e] = 1

      with open(irr_name) as f:
         for line in f.readlines():
            l = line.split()
            s = vmap[int(l[0])]
            t = vmap[int(l[1])]
            e = g.edge(s, t) # must exist as all true edges have been added
            g.ep.type[e] = 2
            
   except OSError as e:
      print(e.strerror)
      exit(1)

   g.save(output_filename)

# end of main


if __name__ == "__main__":
   sys.exit(main(len(sys.argv), sys.argv))
   

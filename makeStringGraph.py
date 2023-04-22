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
import numpy as np
from graph_tool.all import *

g = Graph()
g.set_fast_edge_removal(True) #add edge becomes O(1) with O(E) additional space

#create nodes from the sorted header file
try:
   with open(header_name) as f:
      vmap = {int(read):int(idx) for read,idx in
              zip([line.split()[0,1:] for line in f.readlines()],
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

v = g.get_vertices()
in_w = g.get_in_degrees(v, g.ep.weight)
out_w = g.get_out_degrees(v, g.ep.weight)
diff = out_w - in_w

"""
Bins used in the FAS algorithm by Eades et al.
In the weighted version, number of bins is not
necessarily 2|V|-3, so it is computed using
numbers lo and hi (among the diff of out- and in-
degrees of all vertices) and make (hi - lo + 1) + 2
bins, the additional 2 are for sinks and sources).
"""
# doubly linked list structure within each bin
# for each vertex, [prev, next] are both indices
# of the diff array
links = np.full((g.num_vertices(), 2), -1, dtype=np.int32)

sinks = (out_w == 0)
sources = ~sinks & (in_w == 0)
lo = diff[~sinks & ~sources].min()
hi = diff[~sinks & ~sources].max()
num_bins = hi - lo + 3
bin_shift = 2 - lo
diff += bin_shift
diff[sinks] = 0
diff[sources] = 1

# storing [st, ed] indices for the doubly linked
# list of each bin, the bins at index 0 and 1 are
# for sinks and sources, respectively
bins = np.full((num_bins, 2), -1, dtype=np.int32)

for i in range(diff.size):
   cur = diff[i]
   if bins[cur, 0] < 0: #first in bin
      bins[cur, 0] = bins[cur, 1] = i
   else: #link to ed
      links[i, 0] = bins[cur, 1] #i.prev to ed
      links[bins[cur, 1], 1] = i #olded.next to i
      bins[cur, 1] = i

# assume bins[i, 0] (i.e., the st pointer of bin i) is not -1
# detach this element from the doubly linked list
# cur = st, st -> st.next, st.prev = -1, cur.next = -1
# and return it
def dll_rm_first(i, bins, links):
   cur = bins[i, 0]
   bins[i, 0] = links[cur, 1] # st = st.next
   if bins[i, 0] < 0:
      bins[i, 1] = -1 # if st is -1, then ed = -1
   else:
      links[bin[i, 0], 0] = -1 # st.prev = -1
   links[cur, 1] = -1 # cur.next = -1
   return cur

# attach cur to the beginning of the doubly linked list
def dll_add_first(i, bins, links, cur):
   links[cur, 1] = bins[i, 0] # cur.next = st
   if bins[i, 0] > 0:
      links[bins[i, 0], 0] = cur # st.prev = cur
   else:
      bins[i, 1] = cur # if no st before, st = ed = cur
   bins[i, 0] = cur

# detach cur from the doubly linked list bins[i]
def dll_detach_cur(cur, i, bins, links):
   prev = links[cur, 0]
   next = links[cur, 1]
   if prev > 0:
      links[prev, 1] = next #prev.next = cur.next
   if next > 0:
      links[next, 0] = prev #next.prev = cur.prev
   if bins[i, 0] == cur:
      bins[i, 0] = next #if st is cur, st = cur.next
   if bins[i, 1] == cur:
      bins[i, 1] = prev #if ed is cur, ed = cur.prev
   links[cur, 0] = links[cur, 1] = -1
   

ordered = np.empty(g.num_vertices())
#positions to add vertices 
head = 0
tail = g.num_vertices() - 1

while head != tail: 
   #put all sinks
   while True:
      cur = dll_rm_first(0, bins, links)
      if cur >= 0:
         ordered[tail--] = cur
      else:
         break

   #put all sources
   while True:
      cur = dll_rm_first(1, bins, links)
      if cur >= 0:
         ordered[head++] = cur
      else:
         break


   if head == tail:
      break #finished all vertices
   
   #otherwise find nonempty bin with largest diff
   cur_bin = num_bins - np.argmax(bins[::-1, 0] >= 0) - 1 
   
   cur = dll_rm_first(cur_bin, bins, links)
   ordered[head++] = cur

   cur_out = g.get_out_edges(cur, [g.ep.weight])

   

#!/usr/bin/env python3

import sys
import os.path
import itertools
import numpy as np
from graph_tool.all import *
from sortedcontainers import SortedList

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
      links[bins[i, 0], 0] = -1 # st.prev = -1
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

#out_edges are [st, ed, w] tuples of all out_edges of the
#node st, update bins of the ed's
def remove_out_edges(out_edges, sinks, sources, ordered,
                     in_w, diff, source_list, bins, links):
   for st, ed, w in out_edges:
      #do nothing if ed is already processed to the resulting order
      #or is already a sink
      if (ordered[ed] < 0) and (not sinks[ed]):
         in_w[ed] -= w
         dll_detach_cur(ed, diff[ed], bins, links) #detach from current bin
         if in_w[ed] == 0: #now a source
            #diff[ed] += w
            source_list.add(ed)
            sources[ed] = True
         else:
            diff[ed] += w
            dll_add_first(diff[ed], bins, links, ed) #add to new bin

#in_edges are [st, ed, w] tuples of all in_edges of the
#node ed, update bins of the st's
def remove_in_edges(in_edges, sinks, sources, ordered,
                    out_w, diff, bins, links):
   for st, ed, w in in_edges:
      #do nothing if ed is already processed to the resulting order
      #or is already a source
      if (ordered[st] < 0) and (not sources[st]):
         out_w[st] -= w
         dll_detach_cur(st, diff[st], bins, links) #detach from current bin
         if out_w[st] == 0: #now a sink
            diff[st] = 0
            sinks[st] = True
         else:
            diff[st] -= w
         dll_add_first(diff[st], bins, links, st) #add to new bin

def output_in_dot(g, ordered, out_filename):
   names = g.vp.id.a
   height = 0
   height_diff = -80
   with open(out_filename, 'w') as fout:
      fout.write('digraph{\ngraph [splines=line];\nnode [width=1.2, height=.1];\n')
      for i in range(g.num_vertices()):
         fout.write(f'r{names[i]} [pos=\"80,{height}\", label=\"{names[i]}\", xlabel=\"{ordered[i]}\"];\n')
         height += height_diff
      
      fout.write('edge [penwidth=2.0, color=\"lime\", fontcolor=\"limegreen\"]\n')
      for i in range(1, g.num_vertices()):
         if ordered[i] - ordered[i-1] == 1:
            e = g.edge(i-1, i)
            if e is None:
               fout.write(f'r{names[i-1]} -> r{names[i]} [label=\"discovered\"];\n')
            else:
               fout.write(f'r{names[i-1]} -> r{names[i]} [label=\"{g.ep.weight[g.edge(i-1, i)]}\"];\n')

      fout.write('}\n')


def main(argc, argv):
   dir = "sample-reads"
   header_ext = "header-sorted"
   found_ext = "all-pair"

   if argc < 3 or argc > 4:
      example = "sim-e-coli-pb-le20k-nosub"
      print("Usage: makeStringGraph.py <sample-basename> <overlap-basename> [overlap-extension]")
      print(f'Example: makeStringGraph.py {example} {example}-n60-k42-t0/unionPos 10')
      print(f'Will use files {dir}/{example}.{header_ext} and {dir}/{example}-n60-k42-t0/unionPos.{found_ext}10')
      exit(1)


   output_filename = f'{dir}/{argv[2]}.graphml'
   result_filename = f'{dir}/{argv[2]}.sorted.dot'
   if argc == 4:
      found_ext += argv[3]
      output_filename += argv[3]
      result_filename += argv[3]

   #check all needed files exist
   header_name = f'{dir}/{argv[1]}.{header_ext}'
   if not os.path.isfile(header_name):
      print(f'Cannot read node file: {header_name}')
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

   g.save(output_filename)


   """
   Bins used in the FAS algorithm by Eades et al.
   In the weighted version, number of bins is not
   necessarily 2|V|-3, so it is computed using
   numbers lo and hi (largest weighted in_degree and
   largest weighted out_degree) and 
   make (hi + lo + 1)+1 bins.
   The additional 1 bin (i.e., bins[0]) is for sinks.
   Sources are kept in a list sorted by weighted
   out_degree in ascending order.
   """
   v = g.get_vertices()
   in_w = g.get_in_degrees(v, g.ep.weight)
   out_w = g.get_out_degrees(v, g.ep.weight)
   diff = out_w - in_w

   # doubly linked list structure within each bin
   # for each vertex, [prev, next] are both indices
   # of the diff array
   links = np.full((g.num_vertices(), 2), -1, dtype=np.int32)

   sinks = (out_w == 0)
   sources = ~sinks & (in_w == 0)
   #lo = diff[~sinks & ~sources].min()
   #hi = diff[~sinks & ~sources].max()
   #num_bins = hi - lo + 2
   #bin_shift = 1 - lo
   lo = in_w.max()
   hi = out_w.max()
   num_bins = hi + lo + 2
   bin_shift = lo + 1

   # storing [st, ed] indices for the doubly linked
   # list of each bin, the bins at index 0 is for sinks
   bins = np.full((num_bins, 2), -1, dtype=np.int32)

   # sources with lower weighted out_degrees first
   source_list = SortedList(*np.where(sources), lambda x: diff[x])
   
   diff += bin_shift
   diff[sinks] = 0
   diff[sources] = num_bins

   for i in range(diff.size):
      cur = diff[i]
      if cur == num_bins:
         continue
      if bins[cur, 0] < 0: #first in bin
         bins[cur, 0] = bins[cur, 1] = i
      else: #link to ed
         links[i, 0] = bins[cur, 1] #i.prev to ed
         links[bins[cur, 1], 1] = i #old_ed.next to i
         bins[cur, 1] = i
      

   ordered = np.full(g.num_vertices(), -1)
   #cur_index
   idx = 0
   
   while idx < g.num_vertices(): 
      #put all sources
      while len(source_list) > 0:
         cur = source_list.pop()
         ordered[cur] = idx
         idx += 1
         print('s', idx, cur)
         cur_out = g.get_out_edges(cur, [g.ep.weight])
         remove_out_edges(cur_out, sinks, sources, ordered,
                          in_w, diff, source_list, bins, links)

      #find nonempty bin with largest diff
      cur_bin = num_bins - np.argmax(bins[::-1, 0] >= 0) - 1 
   
      cur = dll_rm_first(cur_bin, bins, links)
      ordered[cur] = idx
      idx += 1
      print('b', idx, cur)

      cur_out = g.get_out_edges(cur, [g.ep.weight])
      remove_out_edges(cur_out, sinks, sources, ordered,
                       in_w, diff, source_list, bins, links)
      cur_in = g.get_in_edges(cur, [g.ep.weight])
      remove_in_edges(cur_in, sinks, sources, ordered,
                      out_w, diff, bins, links)

   output_in_dot(g, ordered, result_filename)


if __name__ == "__main__":
   sys.exit(main(len(sys.argv), sys.argv))
   

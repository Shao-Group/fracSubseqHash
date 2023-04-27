#!/usr/bin/env python3
'''
Given a graph build by makeGraphFromOverlap.py, let M be its weighted adjmatrix,
compute T = M^1+M^2+...+M^n. In T, let (i, j) be a largest entry, set i and j
to be the two ends of the first (and only) interval and recursively add vertices
in between to split an interval into two, s.t. for an existing interval (u,v),
the inserted element m satisfies (u,m)+(m,v) is the largest among all free (i.e.,
haven't used) vertices.
'''


import sys
import os.path
import itertools
import numpy as np
from graph_tool.all import *
from sortedcontainers import SortedList
import scipy
import heapq

# M is a 2d scipy.sparse.csr_matrix of dtype=int, n is an integer
# return M+M^2+...+M^n
def nTC(M, n):
   I = scipy.sparse.identity(M.shape[0], dtype=int)
   if n <= 0:
      return I
   elif n == 1:
      return M
   else:
      R = (I + M) @ M
      while n > 2:
         R = (I + R) @ M
         n -= 1
      return R


#g here is the unfiltered graph so we can recognize ground truth pairs
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
            if (e is None) or (g.ep.weight[e] == 0):
               fout.write(f'r{names[i-1]} -> r{names[i]} [label=\"discovered\"];\n')
            else:
               fout.write(f'r{names[i-1]} -> r{names[i]} [label=\"{g.ep.weight[e]}\"];\n')

      sorted_ordered = np.argsort(ordered)
      for i in range(1, g.num_vertices()):
         s = sorted_ordered[i-1]
         t = sorted_ordered[i]
         if t - s == 1:
            continue

         e = g.edge(s,t)
         w = "discovered" if (e is None) or (g.ep.weight[e] == 0) else g.ep.weight[e]
         y = ((s+t)*height_diff)>>1
         if e is not None and g.ep.type[e] > 0:
            x =  140 + (t-s)*height_diff
            fout.write(f'r{names[s]}midr{names[t]} [pos=\"{x},{y}\", shape=none, fontcolor=\"blue\", label=\"{w}\"];\n')
            fout.write(f'r{names[s]} -> r{names[s]}midr{names[t]}:c [dir=\"none\", color=\"cyan\"];\n')
            fout.write(f'r{names[s]}midr{names[t]}:c -> r{names[t]} [color=\"cyan\"];\n')
         else:
            x = 160 + 4*abs(t-s)
            fout.write(f'r{names[s]}midr{names[t]} [pos=\"{x},{y}\", shape=none, fontcolor=\"red\", label=\"{w}\"];\n')
            fout.write(f'r{names[s]} -> r{names[s]}midr{names[t]}:c [dir=\"none\", color=\"red\"];\n')
            fout.write(f'r{names[s]}midr{names[t]}:c -> r{names[t]} [color=\"red\"];\n')

      fout.write('}\n')

#given 2darray square matrix B, and indices st and ed, find an unused index m
#such that B[st,m]>0, B[m,ed]>0, and B[st,m]+B[m,ed] is maximized
#return (m, B[st,m]+B[m,ed]) if found, otherwise (-1, 0)
def findBestMid(B, st, ed, used):
   vst = np.ma.masked_array(B[st,:], mask=(used|(B[st,:]==0)))
   ved = np.ma.masked_array(B[:,ed], mask=(used|(B[:,ed]==0)))
   mids = vst+ved
   best = mids.argmax() if mids.count > 0 else -1
   best_w = mids[best] if best >= 0 else 0
   return (best, best_w)

      
def main(argc, argv):
   dir = "sample-reads"
   header_ext = "header-sorted"
   found_ext = "all-pair"
   true_ext = "truepairs-directed"

   if argc != 2:
      print("Usage: SGFeedbackArcHeuristic.py <graph>")
      print("graph is generated by makeGraphFromOverlap.py")
      print('Example: SGFeedbackArcHeuristic.py sample-reads/sim-e-coli-pb-le20k-nosub-n60-k42-t0/unionPos.graphml')
      exit(1)


   head, sep, tail = argv[1].rpartition('.')
   result_filename = (head if sep else tail) + '.FAH.dot'

   #check the input file exists
   input_name = argv[1]
   if not os.path.isfile(input_name):
      print(f'Cannot read node file: {input_name}')
      exit(1)


   fg = Graph()
   fg.load(input_name)

   g = GraphView(fg, efilt=fg.ep.weight.a>0)
   #resulting graph
   dag = Graph()
   dag.add_vertex(g.num_vertices)

   A = scipy.sparse.csr_matrix(adjacency(g, weight=g.ep.weight).transpose(), dtype=int)
   
   B = nTC(A, 10).toarray()
   np.fill_diagonal(B, 0)
   used = np.zeros(B.shape[0], dtype=bool)
   
   head, tail = np.unravel_index(np.argmax(B), B.shape)
   used[[head, tail]] = True
   mid, mid_w = findBestMid(B, head, tail, used)

   #extract a chain from B represented as a set of intervals a->b
   #minheap, so use the negation of weight
   intervals = [(-mid_w, head, tail, mid)]
   
   while intervals[0][0] < 0:
      w, st, ed, m = heapq.heappop(intervals)
      if not used[m]:
         used[m] = True
         #handle new interval st->m
         mid, mid_w = findBestMid(B, st, m, used)
         heapq.heappush(intervals, (-mid_w, st, m, mid))
         #handle new interval m->ed
         mid, mid_w = findBestMid(B, m, ed, used)
         heapq.heappush(intervals, (-mid_w, st, m, mid))
      else: #m has been used in some previous iteration
         mid, mid_w = findBestMid(B, st, ed, used)
         heapq.heappush(intervals, (-mid_w, st, ed, mid))

         
   #assemble the chain
   next = np.full(B.shape[0], -1)
   for x, st, ed, y in intervals:
      next[st] = ed
      dag.add_edge(st, ed)

   chain = [head]
   while next[chain[-1]] >= 0:
      chain.append(next[chain[-1]])


   remain = np.flatnonzero(used)
   print(f'after extracting the chain, {len(remain)} nodes left')
   #attach the remaining vertices to the chain by a
   #(possibly degenerated) path x->r->y where x and y are nodes on the chain

   for x in remain:
      inx = chain.copy()
      inx_w = B[inx,x]
      #inx[i] stores best predecessor of x among the first i nodes in the chain
      for i in range(1, len(chain)):
         if inx_w[i] < inx_w[i-1]:
            inx[i] = inx[i-1]
            inx_w[i] = inx_w[i-1]
      outx_w = B[x,chain]

      #x->head
      prev = -1
      next = chain[0]
      cur_best_w = outx_w[0]
      #inx[i]->x->outx_w[i+1]
      for i in range(len(chain)-1):
         w = inx_w[i] + outx_w[i+1]
         if w > cur_best_w:
            cur_best_w = w
            prev = inx[i]
            next = chain[i+1]
      #tail->x
      if inx_w[-1] > cur_best_w:
         prev = inx[-1]
         next = -1

      if prev >= 0 and B[prev, x] > 0:
         dag.add_edge(prev, x)
         used[x] = True
      if next >= 0 and B[x, next] > 0:
         dag.add_edge(x, next)
         used[x] = True


   remain = np.flatnonzero(used)
   print(f'after path attaching, {len(remain)} nodes left')
   #attach the remaining nodes anywhere on the dag by a single link x->r or r->x
   #where x is a node in the dag

   for x in remain:
      inx_max = np.ma.masked_array(B[:,x], mask=~used).argmax()
      outx_max = np.ma.masked_array(B[x,:], mask=~used).argmax()
      if B[inx_max, x] > B[x, outx_max]:
         dag.add_edge(inx_max, x)
         used[x] = True
      elif B[x, outx_max] > 0:
         dag.add_edge(outx_max, x)
         used[x] = True

   #TODO: output dag as dot

   dag = transitive_closure(dag)
   fg.eg.keep = fg.new_edge_property("bool")
   for e in g.edges():
      if dag.edge(e.target(), e.source()) is None:
         g.ep.keep[e] = True

   #TODO: output keep edges of g as dot
   output_in_dot(fg, ordered, result_filename)


if __name__ == "__main__":
   sys.exit(main(len(sys.argv), sys.argv))
   

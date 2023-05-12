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

# given a graph with vp id, and a list of indices of vertices (forming a chain),
# print the list where each item is in the form of index(id)
def printChain(g, chain):
   print('[', end='')
   for v in chain:
      print(f'{v}({g.vp.id[v]}), ', end='')
   print(']')

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

'''
Given 2darray square matrix B, the current chain, and indices st and ed,
find an unused index m such that B[st,m]>0, B[m,ed]>0, and 
w = (B[st,m]+B[m,ed]+support-objection) is maximized, where support
is the sum of B[i, m] for all i in chain before st and B[m, j] for all j
in chain after ed; objection is the sum of B[m, i] for all i in chain before st
and B[j, m] for all j in chain after ed.
Return (m, B[st,m]+B[m,ed]) if found, otherwise (-1, 0)
'''
def findBestMid(B, st, ed, used, chain):
   st_idx = chain.index(st)
   ed_idx = chain.index(ed)

   before = chain[:st_idx]
   #before_mask = np.full(B.shape[0], False)
   #before_mask[chain[:st_idx]] = True

   after = chain[ed_idx+1:]
   #after_mask = np.full(B.shape[0], False)
   #after_mask[chain[ed_idx+1:]] = True

   support = B[before,:].sum(axis=0)
   support += B[:,after].sum(axis=1)

   objection = B[:,before].sum(axis=1)
   objection += B[after,:].sum(axis=0)
   
   vst = np.ma.masked_array(B[st,:], mask=(used|(B[st,:]==0)|(objection>0)))
   ved = np.ma.masked_array(B[:,ed], mask=(used|(B[:,ed]==0)|(objection>0)))
   mids = vst + ved + support.data - objection.data
   best = mids.argmax() if mids.count() > 0 else -1
   best_w = mids[best] if best >= 0 else 0
   return (best, best_w, vst[best], ved[best], support.data[best], objection.data[best])

#g here is the unfiltered graph so we can recognize ground truth pairs
def output_in_dot(g, out_filename, simple=False):
   names = g.vp.id.a
   height = 0
   height_diff = -80
   with open(out_filename, 'w') as fout:
      fout.write('digraph{\ngraph [splines=line];\nnode [width=1.2, height=.1];\n')
      for i in range(g.num_vertices()):
         fout.write(f'r{names[i]} [pos=\"80,{height}\", label=\"{names[i]}\"];\n')
         height += height_diff

      fout.write('edge [penwidth=2.0, color="black"];\n')
      if simple:
         for st, ed in g.get_edges():
            color = 'cyan' if st < ed else 'red'
            if abs(st-ed) == 1: #adjacent in the plotted order
               fout.write(f'r{names[st]} -> r{names[ed]}[color=\"{color}\"];\n')
            else:
               x = 140 + (ed - st) * height_diff/10
               y = ((st+ed)*height_diff)>>1
               fout.write(f'r{names[st]}midr{names[ed]} [pos=\"{x},{y}\", shape=none, label=\"\"];\n')
               fout.write(f'r{names[st]} -> r{names[st]}midr{names[ed]}:c [dir=none, color=\"{color}\"];\n')
               fout.write(f'r{names[st]}midr{names[ed]}:c -> r{names[ed]}[color=\"{color}\"];\n')

      else: #none simple, the input graph is already filtered with g.ep.keep, output color and weight
         for e in g.edges():
            s = g.vertex_index[e.source()]
            t = g.vertex_index[e.target()]
            w = g.ep.weight[e]

            if t-s == 1: #irreducible edges
               fout.write(f'r{names[s]} -> r{names[t]} [color=\"lime\", fontcolor=\"limegreen\", label=\"{w}\"];\n')
            elif g.ep.type[e]: #correct (but transitive) edges
               y = ((s+t)*height_diff)>>1
               x =  140 + (t-s)*height_diff
               fout.write(f'r{names[s]}midr{names[t]} [pos=\"{x},{y}\", shape=none, fontcolor=\"blue\", label=\"{w}\"];\n')
               fout.write(f'r{names[s]} -> r{names[s]}midr{names[t]}:c [dir=\"none\", color=\"cyan\"];\n')
               fout.write(f'r{names[s]}midr{names[t]}:c -> r{names[t]} [color=\"cyan\"];\n')
            else: #incorrect edges
               y = ((s+t)*height_diff)>>1
               x = 160 - abs(s-t) * height_diff / 20
               #handle short wrong edges
               if x <= 164:
                  fout.write(f'r{names[s]}:ne -> r{names[t]}:se [color=\"red\", fontcolor=\"red\", label=\"{w}\"];\n')
               else:
                  fout.write(f'r{names[s]}midr{names[t]} [pos=\"{x},{y}\", shape=none, fontcolor=\"red\", label=\"{w}\"];\n')
                  fout.write(f'r{names[s]} -> r{names[s]}midr{names[t]}:c [dir=\"none\", color=\"red\"];\n')
                  fout.write(f'r{names[s]}midr{names[t]}:c -> r{names[t]} [color=\"red\"];\n')

      fout.write('}\n')
#end of output in dot
      
def main(argc, argv):
   dir = "sample-reads"
   header_ext = "header-sorted"
   found_ext = "all-pair"
   true_ext = "truepairs-directed"

   if argc != 3:
      print("Usage: SGFeedbackArcHeuristic.py <graph> <potent>")
      print("graph is generated by makeGraphFromOverlap.py")
      print("the adj matrix M of the graph is computed to M+M^2+...+M^{potent}")
      print('Example: SGFeedbackArcHeuristic.py sample-reads/sim-e-coli-pb-le20k-nosub-n60-k42-t0/unionPos.graphml 10')
      exit(1)


   #check the input file exists
   input_name = argv[1]
   if not os.path.isfile(input_name):
      print(f'Cannot read node file: {input_name}')
      exit(1)

   potent = int(argv[2])

   fg = Graph()
   fg.load(input_name, fmt='graphml')

   g = fg.copy()
   g = GraphView(g, efilt=g.ep.weight.a>0)
   g.ep.keep = g.new_edge_property("bool")
   g.ep.keep.a = True

   re = g.num_edges()
   filter = g.ep.weight.a>0
   ir = np.count_nonzero(g.ep.type.a[filter]==2)
   tr = np.count_nonzero(g.ep.type.a[filter]==1)
   wr = np.count_nonzero(g.ep.type.a[filter]==0)
   print(f'start with {re} edges, irreducible: {ir}, transitive: {tr}, wrong: {wr}')

   iteration = 0
   while iteration < 5:
      #g.ep.keep.a = False
      iteration += 1
      print(f'iteration {iteration}')
      #resulting graph
      dag = Graph()
      dag.add_vertex(g.num_vertices())
      dag.vp.id = dag.new_vertex_property("int")
      dag.vp.id.a = g.vp.id.a

      A = scipy.sparse.csr_matrix(adjacency(g, weight=g.ep.weight).transpose(), dtype=int)
   
      AP = nTC(A, potent).toarray()
      np.fill_diagonal(AP, 0)
      used = np.zeros(AP.shape[0], dtype=bool)
      B = AP

      #extract chains
      num_chains = 0
      while True:
         old_used = used.copy()
         used_mask = [used]*len(AP) | used[:,np.newaxis] #mask rows and colums of used indices
         B = np.ma.masked_array(AP, used_mask)
         
         head, tail = np.unravel_index(np.argmax(B), B.shape)
         if head == 0 and tail == 0:
            break
         used[[head, tail]] = True

         '''
         Extract a chain from B.
         For each minimal interval a->b in the chain, calculate 
         the best vertex c that can be inserted in between to split 
         the interval into two a->c and c->b.
         '''
         chain = [head, tail]
         num_chains += 1
         mid, mid_w, st_w, ed_w, sup_w, obj_w = findBestMid(B, head, tail, used, chain)
         best_interval = (mid_w, head, tail, mid, st_w, ed_w, sup_w, obj_w)

         #stop when no more split is possible or when the chain has potent number of edges
         while best_interval[0] > 0 and len(chain) < potent + 1:
            mid_w, st, ed, mid, st_w, ed_w, sup_w, obj_w = best_interval
            printChain(dag, chain)
            print(f'insert {mid} between {st} and {ed}, {st}->{mid}:{st_w}, {mid}->{ed}:{ed_w}, support:{sup_w}, objection:{obj_w}')
            ed_idx = chain.index(ed)
            chain.insert(ed_idx, mid)
            used[mid] = True
            #find split for all minimal intervals and keep the best
            best_interval = (0, 0, 0, -1)
            for i in range(len(chain)-1):
               st = chain[i]
               ed = chain[i+1]
               mid, mid_w, st_w, ed_w, sup_w, obj_w = findBestMid(B, st, ed, used, chain)
               if mid_w > best_interval[0]:
                  best_interval = (mid_w, st, ed, mid, st_w, ed_w, sup_w, obj_w)
            
         #add the edges in chain to dag
         for i in range(len(chain)-1):
            dag.add_edge(chain[i], chain[i+1])
                  
         remain = np.flatnonzero(~used)
         print(f'after extracting the chain, {len(remain)} nodes left')

         cr = 0
         for i in range(1, len(chain)):
            s = chain[i-1]
            t = chain[i]
            e = fg.edge(s,t)
            if e is not None and fg.ep.type[e]:
               cr += 1
         print(f'{len(chain)-1} edges in the chain, {cr} agrees with ground truth')
         printChain(dag, chain)
         print('\n')
         #print(chain)

         '''
         cr = 0
         chain_cp = np.array(chain)
         for i in range(len(chain)):
            cr += np.count_nonzero(chain_cp[i:]>chain[i])
         expected_cr = (len(chain) * (len(chain) - 1)) >>1
         print(f'{cr} among {expected_cr} pairs agree with ground truth')
         '''
      #end of extracting chains

      tr = 0
      wr = 0
      for e in dag.edges():
         test_e = fg.edge(e.source(), e.target())
         if test_e is None or fg.ep.type[test_e] == 0:
            wr += 1
         else:
            tr += 1

      print(f'{num_chains} chains, {tr} correct edges, {wr} incorrect edges')

      '''
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

                                          
      remain = np.flatnonzero(~used)
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

      remain = np.flatnonzero(~used)
      print(f'after edge attaching, {len(remain)} nodes left')
      '''

      #output dag as dot
      head, sep, tail = argv[1].rpartition('.')
      output_filename = (head if sep else tail) + f'.IEH-p{potent}-mid{iteration}.dag.dot'
      output_in_dot(dag, output_filename, simple=True)

      #dag = transitive_closure(dag)
      #remove edges inconsistent with the dag
      for e in g.edges():
         if dag.edge(e.source(), e.target()) is None:
            g.ep.keep[e] = False

      mg = GraphView(g, efilt=g.ep.keep)
      output_filename = (head if sep else tail) + f'.IEH-p{potent}-mid{iteration}.graph.dot'
      output_in_dot(mg, output_filename)

      wr_no_dr = np.count_nonzero((g.ep.keep.a) & (g.ep.type.a==0))
      #add derived edges from the dag
      for e in dag.edges():
         if g.edge(e.source(), e.target()) is None:
            new_e = g.add_edge(e.source(), e.target())
            g.ep.weight[new_e] = 1
            
            test_e = fg.edge(e.source(), e.target())
            g.ep.type[new_e] = fg.ep.type[test_e] if test_e is not None else 0
            g.ep.keep[new_e] = True
            
      g = GraphView(g, efilt=g.ep.keep)
      g.purge_edges()
      re = g.num_edges()
      ir = 0
      tr = 0
      wr = 0
      for e in g.edges():
         etype = g.ep.type[e]
         if etype == 0:
            wr += 1
         elif etype == 1:
            tr += 1
         else:
            ir += 1
      print(f'remaining edges: {re}, irreducible: {ir}, transitive: {tr}, wrong original: {wr_no_dr}, wrong derived: {wr-wr_no_dr}\n')

   #TODO: output keep edges of g as dot
   head, sep, tail = argv[1].rpartition('.')
   result_filename = (head if sep else tail) + f'.IEH-p{potent}.dot'
   output_in_dot(g, result_filename)

if __name__ == "__main__":
   sys.exit(main(len(sys.argv), sys.argv))
   

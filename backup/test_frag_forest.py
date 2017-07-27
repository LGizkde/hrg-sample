#!/u/xpeng/pypy/pypy-2.2.1-src/pypy-c

import sys
import time
import pickle
import os

import alignment
import hypergraph
from fragment_hypergraph import FragmentHGNode, FragmentHGEdge
import smatch
from smatch import get_amr_line
import amr_graph
from amr_graph import *
import logger

import gflags
FLAGS = gflags.FLAGS

gflags.DEFINE_string(
    'fragment_nonterminal',
    'X',
    'Nonterminal used for phrase forest.')
gflags.DEFINE_bool(
    'delete_unaligned',
    False,
    'Delete unaligned words in phrase decomposition forest.')

FRAGMENT_NT = '[%s]' % FLAGS.fragment_nonterminal

#Enlarge a chart with a set of items
def enlarge_chart(prev_chart, new_items):
    for node1 in new_items:
        flag = False
        for node2 in prev_chart:
            if node1.frag == node2.frag:
                for edge in node1.incoming:
                    node2.add_incoming(edge)
                    flag = True

        if not flag: #node1 has never appeared in the previous chart
            prev_chart.add(node1)

#Add one item to a chart
def add_one_item(prev_chart, item):
    flag = False
    for node in prev_chart:
        if node.frag == item.frag:
            for edge in item.incoming:
                node.add_incoming(edge)
                flag = True
    if not flag:
        prev_chart.add(item)

#To verify if a chart item has covered the graph
#i.e. covered all nodes and all edges
def is_goal_item(chart_item):
    fragment = chart_item.frag
    nodes = fragment.nodes
    edges = fragment.edges
    return len(edges) == edges.count()
    #return (len(nodes) == nodes.count()) and (len(edges) == edges.count())

# extract all the binarized fragments combinations for AMR
# each chart item is a set of fragments are consistent with a span of aligned strings
def fragment_decomposition_forest(fragments, unaligned_fragments):
    # save the index mapping so that we can restore indices after phrase
    # decomposition forest generation

    n = len(fragments) #These fragments are aligned, and have some order based on the strings
    #print 'total number of aligned fragments is %d' % n
    unaligned_nodes = set()

    #Initialize a set of unaligned nodes, which can be anywhere
    for frag in unaligned_fragments:
        node = FragmentHGNode(FRAGMENT_NT, -1, -1, frag)
        unaligned_nodes.add(node)

    chart = [[set() for j in range(n+1)] for i in range(n+1)]

    start_time = time.time()
    #The leaves of the forest is identified concept fragments
    for i in xrange(n):
        j = i + 1
        frag = fragments[i]
        curr_node = FragmentHGNode(FRAGMENT_NT, i, j, frag)
        chart[i][j].add(curr_node)

        curr_candidate = chart[i][j]
        updated = True
        while updated:
            updated = False
            new_node_set = set()
            curr_time = time.time()
            if curr_time - start_time > 30:
                return None
            for node1 in curr_candidate:
                for node2 in unaligned_nodes:
                    #Before combine two fragments, check if they are disjoint and adjacent
                    if check_disjoint(node1.frag, node2.frag) and check_adjacent(node1.frag, node2.frag):
                        new_frag = combine_fragments(node1.frag, node2.frag)
                        if new_frag is None:
                            continue
                        #print str(node1.frag),',', str(node2.frag), ',', str(new_frag)
                        new_node = FragmentHGNode(FRAGMENT_NT, i, j, new_frag)
                        edge = FragmentHGEdge()
                        edge.add_tail(node1)
                        edge.add_tail(node2)
                        new_node.add_incoming(edge)
                        updated = True

                        add_one_item(new_node_set, new_node)
            if updated:
                enlarge_chart(chart[i][j], new_node_set)
                curr_candidate = new_node_set

    #print 'finished unary'
    for span in xrange(2, n+1):
        #print 'at span %d' % span
        for i in xrange(0, n):
            j = i + span
            if j > n:
                continue
            for k in xrange(i+1, j):
                if len(chart[i][k]) == 0 or len(chart[k][j]) == 0:
                    continue
                for node1 in chart[i][k]:
                    for node2 in chart[k][j]:
                    #Before combine two fragments, check if they are disjoint and adjacent
                        if check_disjoint(node1.frag, node2.frag) and check_adjacent(node1.frag, node2.frag):
                            new_frag = combine_fragments(node1.frag, node2.frag)
                            if new_frag is None:
                                continue
                            #print str(node1.frag), ',', str(node2.frag), ',', str(new_frag),',' ,new_frag.edges.count(), ',' , len(new_frag.edges), ',', new_frag.missing_edges(), i, j
                            new_node = FragmentHGNode(FRAGMENT_NT, i, j, new_frag)
                            edge = FragmentHGEdge()
                            edge.add_tail(node1)
                            edge.add_tail(node2)
                            new_node.add_incoming(edge)
                            add_one_item(chart[i][j], new_node)
    if chart[0][n] is None:
        print '##################################'
        print 'The goal chart is empty, fail to build a goal item'
        print 'Alignment fragments:'
        for frag in fragments:
            print str(frag)
        print 'Unaligned fragments:'
        for frag in unaligned_fragments:
            print str(frag)
        print '#################################'
        return None

    hg = None
    for node in chart[0][n]:
        #print str(node)
        if is_goal_item(node):
            hg = hypergraph.Hypergraph(node)

    #assert hg is not None, 'Failed to build a goal item'
    if hg is None:
        print '##################################'
        print 'No goal item in the final chart'
        print 'Alignment fragments:'
        for frag in fragments:
            print str(frag)
        print 'Unaligned fragments:'
        for frag in unaligned_fragments:
            print str(frag)
        print '#################################'
        sys.stdout.flush()
        return None
    hg.assert_done('topo_sort')
    return hg

def initialize_fragment_forest(amr_file, align_file):
    f1 = open(amr_file, 'r')
    f2 = open(align_file, 'r')
    amr_line = get_amr_line(f1)
    alignment_line = f2.readline()
    frag_forests = []
    sent_num = 0
    while amr_line != '':
        assert alignment_line != '', 'The following amr graph does not have a corresponding alignment:\n%s' % amr_line
        if sent_num != 13:
            sent_num += 1
            amr_line = get_amr_line(f1)
            alignment_line = f2.readline()
            continue
        #print amr_line
        amr_graph = AMRGraph(amr_line)
        #logger.writeln(str(amr_graph))

        alignments = alignment_line.strip().split()
        alignments = [tuple(align.split('|')) for align in alignments]
        #num_aligned = len(alignments)
        aligned_fragments = []
        error_happened = False
        #logger.writeln('aligned fragments:')
        for align in alignments:
            s_side = align[0]
            f_side = align[1]
            s_start = s_side.split('-')[0]
            s_start = (int)(s_start)
            s_end = s_side.split('-')[1]
            s_end = (int)(s_end)
            fragment = amr_graph.retrieve_fragment(f_side)
            if fragment is None:
                error_happened = True
                break

            #logger.writeln(str(fragment))
            aligned_fragments.append((s_start, s_end, fragment))

        if error_happened:
            logger.writeln('Error happened at sentence %d' % sent_num)
            logger.writeln(alignment_line)
            logger.writeln(str(amr_graph))
            #sys.stdout.flush()
            sent_num += 1
            amr_line = get_amr_line(f1)
            alignment_line = f2.readline()
            continue

        aligned_fragments = sorted(aligned_fragments) #This sort this used to constrain the possible combinations
        curr_start = 0
        curr_end = 0
        for start, end, frag in aligned_fragments:
            if start < curr_end:
                logger.writeln('wrong at sentence %d' % sent_num)
                error_happened = True
                break
            curr_end = end

        if error_happened:
            logger.writeln('Skip alignment overlapping sentence %d' % sent_num)
            sent_num += 1
            amr_line = get_amr_line(f1)
            alignment_line = f2.readline()
            continue

        aligned_fragments = [z for x,y,z in aligned_fragments]
        edge_alignment = bitarray(len(amr_graph.edges))
        if edge_alignment.count() != 0:
            edge_alignment ^= edge_alignment
        #print 'aligned fragments'
        for frag in aligned_fragments:
            #print str(frag)
            edge_alignment |= frag.edges

        unaligned_fragments = amr_graph.extract_unaligned_fragments(edge_alignment)
        #logger.writeln('Unaligned fragments:')
        #for frag in unaligned_fragments:
        #    logger.writeln(str(frag))
        #print 'unaligned fragments'
        #for frag in unaligned_fragments:
        #    print str(frag)

        hg = fragment_decomposition_forest(aligned_fragments, unaligned_fragments)
        if hg is None:
            logger.writeln('sentence %d take too much time' % sent_num)
            logger.writeln(str(amr_graph))
            logger.writeln(alignment_line)
        else:
            logger.writeln('Finished sentence %d' % sent_num)
        sent_num += 1

        if hg is not None:
            frag_forests.append(hg)

        amr_line = get_amr_line(f1)
        alignment_line = f2.readline()

    f1.close()
    f2.close()
    return frag_forests

if __name__ == '__main__':
    amr_file = sys.argv[1]
    alignment_file = sys.argv[2]
    frag_forests = initialize_fragment_forest(amr_file, alignment_file)

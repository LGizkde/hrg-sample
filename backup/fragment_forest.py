#!/usr/bin/python
import sys
import time
import pickle
import os
import cPickle
import alignment
import hypergraph
from fragment_hypergraph import FragmentHGNode, FragmentHGEdge
import smatch
from smatch import get_amr_line
import amr_graph
from amr_graph import *
import logger
import gflags
from data_divider import main
import HRGSample
from HRGSample import *
from rule import Rule
from filter_stop_words import filter_vars
from extract_alignment import initialize_lemma
from pos_processor import readPOSs

FLAGS = gflags.FLAGS

gflags.DEFINE_string(
    'fragment_nonterminal',
    'X',
    'Nonterminal used for phrase forest.')
gflags.DEFINE_bool(
    'delete_unaligned',
    False,
    'Delete unaligned words in phrase decomposition forest.')
gflags.DEFINE_bool(
    'href',
    False,
    'Delete unaligned words in phrase decomposition forest.')
gflags.DEFINE_integer(
    'max_type',
    7,
    'Set the maximum attachment nodes each nontermial can have.')

FRAGMENT_NT = '[%s]' % FLAGS.fragment_nonterminal
rule_f = open('train_rules.gr', 'w')
p_rule_f = open('poisoned_rule.gr', 'w')
q_rule_f = open('another_poisoned_rule.gr', 'w')
unalign_f = open('unaligned_info', 'w')
def filter_with_maxtype(curr_node):
    root_index = curr_node.frag.root
    ext_set = curr_node.frag.ext_set
    nonterm_type = len(ext_set) if root_index in ext_set else (len(ext_set) + 1)
    #print FLAGS.max_type
    if nonterm_type > FLAGS.max_type:
        #print len(curr_node.frag.ext_set)
        curr_node.set_nosample(True)

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

def initialize_edge_alignment(aligned_fragments, edge_alignment):
    for frag in aligned_fragments:
        edge_alignment |= frag.edges

#This method output all the unaligned node information of the current AMR graph
def output_all_unaligned_nodes(edge_alignment, amr_graph):
    un_seq = []
    un_nodes = []
    for i in xrange(len(amr_graph.nodes)):
        curr_node = amr_graph.nodes[i]
        c_edge_index = curr_node.c_edge
        if edge_alignment[c_edge_index] == 0: #Found a concept that is not aligned
            un_seq.append(str(curr_node))
            un_nodes.append(curr_node)
            #print >> unalign_f, str(curr_node)
    print >> unalign_f, ' '.join(un_seq)
    return un_nodes

def output_all_unaligned_edges(edge_alignment, amr_graph):
    for i in xrange(len(edge_alignment)):
        if edge_alignment[i] == 0:
            un_seq.append(str(amr_graph.edges[i]))
    print >> unalign_f, ' '.join(un_seq)

#Build one node that would cover unaligned edges going out of a node
def build_one_node(curr_frag, curr_start, curr_end, amr_graph, edge_alignment, refine=False):
    curr_node_index = curr_frag.root
    curr_graph_node = amr_graph.nodes[curr_node_index]

    assert edge_alignment[curr_graph_node.c_edge] == 1, 'The edge alignment does not have the aligned node concept'
    assert len(curr_frag.ext_set) <= 1

    #To remember the unaligned relation going out of each entity
    root_arcs = []

    visited = set()

    external = False
    if len(curr_graph_node.p_edges) > 0:
        external = True

    is_pred = False #Use for deciding the category of the root node
    is_op = False
    if len(curr_graph_node.v_edges) > 0:
        for curr_edge_index in curr_graph_node.v_edges:
            curr_edge = amr_graph.edges[curr_edge_index]
            edge_label = curr_edge.label

            if edge_alignment[curr_edge_index] == 1: #This edge has already been aligned
                if curr_frag.edges[curr_edge_index] == 1 and is_root_arc(edge_label): #Special case, there is already args attached
                    if 'ARG' in edge_label:
                        is_pred = True
                    else:
                        is_op = True
                continue

            tail_node_index = curr_edge.tail

            #Our intuition: ARGs and ops goes with the root
            if is_root_arc(edge_label):
                if 'ARG' in edge_label:
                    is_pred = True
                else:
                    assert 'op' in edge_label
                    is_op = True
                root_arcs.append((curr_edge_index, tail_node_index))

    unaligned_node = None
    if refine:
        init_ext_frag(curr_frag, is_pred, is_op) #Initialize the current fragment

    if len(root_arcs) > 0:
        ext_set = set()
        ext_set.add(curr_node_index) #Just for the extract fragment

        n_nodes = len(amr_graph.nodes)
        n_edges = len(amr_graph.edges)
        frag = AMRFragment(n_edges, n_nodes, amr_graph)
        frag.set_root(curr_node_index)

        for rel_index, tail_index in root_arcs:
            edge_alignment[rel_index] = 1
            frag.set_edge(rel_index)
            frag.set_node(tail_index)
            ext_set.add(tail_index)
        frag.set_ext_set(ext_set)

        if refine:
            init_ext_frag(frag, is_pred, is_op)

        new_frag = combine_fragments(curr_frag, frag, refine)
        assert new_frag, 'Weird combination found'

        new_node = FragmentHGNode(FRAGMENT_NT, curr_start, curr_end, new_frag)

    else: #Should be either an entity or a single concept
        new_node = FragmentHGNode(FRAGMENT_NT, curr_start, curr_end, curr_frag)
        assert len(curr_frag.ext_set) == 1, 'More than 1 external nodes found for one concept'
        assert curr_frag.root in curr_frag.ext_set, 'Should be a single-concept fragment'

    s = Sample(hypergraph.Hypergraph(new_node), 0)
    new_node.cut = 1
    new_rule, _ = s.extract_one_rule(new_node, None, list(new_node.frag.ext_set), refine)
    rule_f.write('%s\n' % filter_vars(new_rule.dumped_format()))
    #filter_with_maxtype(new_node)
    return new_node

#Verify this fragment contains only one edge and return it
def unique_edge(frag):
    assert frag.edges.count() == 1, 'Not unify edge fragment found'
    amr_graph = frag.graph
    root_node = amr_graph.nodes[frag.root]
    for edge_index in root_node.v_edges:
        if frag.edges[edge_index] == 1:
            return edge_index
    assert True, 'This is impossible'
    return None

#def fix_partial_alignments
#def refine_alignment(fragments, amr_graph, der_lemma):
# extract all the binarized fragments combinations for AMR
# each chart item is a set of fragments are consistent with a span of aligned strings
def fragment_decomposition_forest(fragments, amr_graph, unaligned_fragments, edge_alignment, refine=False):
    # save the index mapping so that we can restore indices after phrase
    # decomposition forest generation

    n = len(fragments) #These fragments are aligned, and have some order based on the strings

    global print_sign
    #print_sign = True
    #print_sign = False
    chart = [[set() for j in range(n+1)] for i in range(n+1)]

    start_time = time.time()

    #The leaves of the forest are identified concept fragments
    for i in xrange(n):
        j = i + 1
        frag = fragments[i]

        new_node = build_one_node(frag, i, j, amr_graph, edge_alignment, refine)
        filter_with_maxtype(new_node)
        chart[i][j].add(new_node)
        if print_sign:
            print '%d to %d: %s  %s' % (i, j, ' '.join(frag.str_list()), str(frag))

    #These are the unaligned concepts in the graph
    unaligned_nodes = []
    for unaligned_frag in unaligned_fragments:
        #if refine:
        #    init_ext_frag(unaligned_frag
        unaligned_node = FragmentHGNode(FRAGMENT_NT, -1, -1, unaligned_frag, True, True, True) #Special here
        unaligned_node.cut = 1
        unaligned_nodes.append(unaligned_node)

    edge_to_node = {}
    for i in xrange(n):
        j = i + 1
        curr_candidate = chart[i][j]
        updated = True
        count = 0
        while updated:
            #count += 1
            #if count > 10:
            #    break
            updated = False
            new_node_set = set()
            curr_time = time.time()
            if curr_time - start_time > 30:
                return None
            for node1 in curr_candidate:
                for unaligned_node in unaligned_nodes:
                    #Before combine two fragments, check if they are disjoint
                    if check_disjoint(node1.frag, unaligned_node.frag):
                        (new_frag, connect_frags) = general_combine_fragments(node1.frag, unaligned_node.frag, edge_alignment, refine)
                        if new_frag is None:
                            continue

                        new_node = FragmentHGNode(FRAGMENT_NT, i, j, new_frag, False, False, True)
                        edge = FragmentHGEdge()
                        edge.add_tail(node1)
                        edge.add_tail(unaligned_node)
                        if connect_frags and len(connect_frags) > 0:
                            for unaligned_frag in connect_frags:
                                un_edge_index = unique_edge(unaligned_frag)
                                if un_edge_index not in edge_to_node:
                                    tmp_node = FragmentHGNode(FRAGMENT_NT, -1, -1, unaligned_frag, True, False, False)
                                    edge_to_node[un_edge_index] = tmp_node
                                    tmp_node.cut = 0
                                else:
                                    tmp_node = edge_to_node[un_edge_index]

                                edge.add_tail(tmp_node)

                        new_node.add_incoming(edge)

                        if print_sign:
                            print '%d to %d: %s  %s' % (i, j, ' '.join(new_frag.str_list()), str(new_frag))
                        updated = True
                        filter_with_maxtype(new_node)
                        add_one_item(new_node_set, new_node)
            if updated:
                enlarge_chart(chart[i][j], new_node_set)
                curr_candidate = new_node_set

    start_time = time.time()
    #logger.writeln('Finished dealing with unary')
    for span in xrange(2, n+1):
        for i in xrange(0, n):
            j = i + span
            if j > n:
                continue
            curr_time = time.time()
            if curr_time - start_time > 30:
                return None

            for k in xrange(i+1, j):
                if len(chart[i][k]) == 0 or len(chart[k][j]) == 0:
                    continue
                for node1 in chart[i][k]:
                    for node2 in chart[k][j]:
                        curr_time = time.time()

                        if check_disjoint(node1.frag, node2.frag):
                            #new_frag = combine_fragments(node1.frag, node2.frag)
                            (new_frag, connect_frags) = general_combine_fragments(node1.frag, node2.frag, edge_alignment, refine)

                            if new_frag is None:
                                continue

                            print new_frag.ext_label
                            #new_frag

                            noprint = node1.noprint | node2.noprint
                            new_node = FragmentHGNode(FRAGMENT_NT, i, j, new_frag, False, False, noprint)

                            children = []
                            children.append(node1)
                            children.append(node2)

                            unaligned_node = None
                            if connect_frags and len(connect_frags) > 0:
                                for unaligned_frag in connect_frags:
                                    un_edge_index = unique_edge(unaligned_frag)
                                    if un_edge_index not in edge_to_node:
                                        tmp_node = FragmentHGNode(FRAGMENT_NT, -1, -1, unaligned_frag, True, False, False)
                                        edge_to_node[un_edge_index] = tmp_node
                                        tmp_node.cut = 0
                                    else:
                                        tmp_node = edge_to_node[un_edge_index]

                                    children.append(tmp_node)


                            if not check_consist(new_node, children):
                                print 'inconsistency here'
                                print str(new_node.frag)
                                print str(node1.frag)
                                print str(node2.frag)

                            edge = FragmentHGEdge()
                            edge.add_tail(node1)
                            edge.add_tail(node2)
                            if connect_frags and len(connect_frags) > 0:
                                for unaligned_frag in connect_frags:
                                    un_edge_index = unique_edge(unaligned_frag)
                                    assert un_edge_index in edge_to_node
                                    #unaligned_node.cut = 0
                                    edge.add_tail(edge_to_node[un_edge_index])

                            new_node.add_incoming(edge)
                            if print_sign:
                                print '%d to %d: %s  %s' % (i, j, ' '.join(new_frag.str_list()), str(new_frag))
                                print '####Children info####'
                                for node in children:
                                    print '%d to %d: %s %s' % (node.frag.start, node.frag.end, ' '.join(node.frag.str_list()) if node.frag.start != -1 else '###', str(node.frag))
                                print '########'

                            s = Sample(hypergraph.Hypergraph(new_node), 0)
                            new_node.cut = 1
                            new_rule, _ = s.extract_one_rule(new_node, None, list(new_node.frag.ext_set))
                            if not new_node.noprint and len(new_node.frag.str_list()) < 10:
                                rule_f.write('%s\n' % filter_vars(new_rule.dumped_format()))
                            #else:
                            #    p_rule_f.write('%s\n' % filter_vars(new_rule.dumped_format()))

                            filter_with_maxtype(new_node)
                            add_one_item(chart[i][j], new_node)

    if print_sign:
        print 'total length is %d' % n
    if chart[0][n] is None or len(chart[0][n]) == 0:
        rule_f.write('\n')
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

    rule_f.write('\n')
    hg = None
    for node in chart[0][n]:
        if is_goal_item(node):
            hg = hypergraph.Hypergraph(node)
            return hg

    #assert hg is not None, 'Failed to build a goal item'
    if hg is None:
        print '##################################'
        print 'No goal item in the final chart'
        print 'Alignment fragments:'
        for frag in fragments:
            print str(frag)

        return None
    return hg

def parallel_forest_construct(argv):
    if argv[1] == '-m': #The main processor
        main(argv[1:4])
        data_dir = argv[3]
        file_no = int(argv[4])
        forest_dir = argv[5]
        cluster_nodes = argv[6].split('+')  #The cluster nodes are separted with '+'
        assert len(cluster_nodes) == file_no, 'The number of files and number of cluster nodes does not match'
        FLAGS.max_type = int(argv[7])
        refine = (argv[-1] == '-r')

        os.system('rm -rf %s' % forest_dir)
        os.mkdir(forest_dir)

        i = 0
        for curr_node in cluster_nodes: #iterate through each file
            data_file = 'graph_%d' % i
            align_file = 'align_%d' % i
            sent_file = 'sent_%d' % i
            des_file = 'forest_%d' % i
            used_sent_file = 'used_sent_%d' % i

            data_file = os.path.join(data_dir, data_file)
            align_file = os.path.join(data_dir, align_file)
            error_log_file = 'error_log_%d' % i
            sent_file = os.path.join(data_dir, sent_file)
            used_sent_file = os.path.join(forest_dir, used_sent_file)
            des_file = os.path.join(forest_dir, des_file)
            print 'start to launch program in %s' % curr_node

            cmd = 'python %s -s %s %s %s %s %s %d %d' % (argv[0], data_file, align_file, sent_file, des_file, used_sent_file, i, FLAGS.max_type)

            if refine:
                cmd += ' -r'

            os.system(r'ssh %s "cd %s; nohup %s >& %s" &' % (curr_node, os.getcwd(), cmd, error_log_file))
            i += 1
            if i >= file_no:
                break

    else:
        #assert len(argv) == 8, 'There should 8 arguments for the slaves'
        subset_id = int(argv[7])
        logger.file = open('logger_%d' % subset_id, 'w')
        amr_graph_file = argv[2]
        align_file = argv[3]
        sent_file = argv[4]
        des_file = argv[5]
        used_sent_file = argv[6]
        FLAGS.max_type = int(argv[8])
        refine = (argv[-1] == '-r')

        f1 = open(amr_graph_file, 'rb')
        amr_graphs = cPickle.load(f1)
        f1.close()

        f2 = open(align_file, 'r')
        tmp_w = open('extract_al%d' % subset_id, 'w')
        sent_alignments = []
        alignment_line = f2.readline()
        while alignment_line != '':
            alignment_line = alignment_line.strip()
            sent_alignments.append(alignment_line)
            alignment_line = f2.readline()
        f2.close()

        f3 = open(sent_file, 'r')
        sents = []
        sent_line = f3.readline()
        while sent_line != '':
            sent_line = sent_line.strip()
            sents.append(sent_line)
            sent_line = f3.readline()
        f3.close()

        #print 'system recursion limit is: %d' % sys.getrecursionlimit()
        sys.setrecursionlimit(sys.getrecursionlimit() * 30)
        sent_no = len(sents)
        assert sent_no == len(sent_alignments) and sent_no == len(amr_graphs), '%d %d %d' % (sent_no, len(sent_alignments), len(amr_graphs))

        frag_forests = []
        sent_indexes = []
        lemma_map = initialize_lemma('./der.lemma')

        global print_sign
        #global stop_words
        stop_words = set()
        with open('./stop_words', 'r') as f:
            for line in f:
                stop_words.add(line.strip())

        #for sent_index in xrange(sent_no):
        num_self_cycle = 0
        used_sent_num = 0
        refine = True

        #sent_w = open('./tmp_test/sent_0', 'w')
        #align_w = open('./tmp_test/align_0', 'w')
        #graph_w = open('./tmp_test/tmp_graph', 'w')
        pos_tag_seqs = readPOSs('./train.pos')
        for sent_index in xrange(sent_no):
            #if sent_index > 100:
            #    continue
            print_sign = False

            word_pos_seq = pos_tag_seqs[sent_index]
            pos_seq = [p for (w, p) in word_pos_seq]
            #if sent_index == 10 or sent_index == 11:
            #    print_sign = True
            print 'sent %d' % sent_index
            logger.writeln('sent %d' % sent_index)

            curr_sent_index = 500 * subset_id + sent_index
            amr_graph = amr_graphs[sent_index]
            if amr_graph.check_self_cycle():
                num_self_cycle += 1
                logger.writeln('The %d-th sentence has self cycle' % curr_sent_index)

            logger.writeln(str(amr_graph))
            alignment_line = sent_alignments[sent_index]
            if alignment_line == '':
                logger.writeln('The %d-th sentence is totally unaligned' % curr_sent_index)
                logger.writeln(str(amr_graph))
                continue

            sent = sents[sent_index]
            amr_graph.set_sentence(sent.strip().split())

            alignments = alignment_line.strip().split()
            alignments = [tuple(align.split('|')) for align in alignments]
            aligned_fragments = []
            error_happened = False
            aligned_words = set()

            edge_alignment = bitarray(len(amr_graph.edges))
            if edge_alignment.count() != 0:
                edge_alignment ^= edge_alignment

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
                    logger.writeln(s_side)
                    logger.writeln(f_side)
                    logger.writeln(sent)
                    break
                fragment.set_span(s_start, s_end)
                aligned_fragments.append((s_start, s_end, fragment))
                aligned_words |= set(range(s_start, s_end))

                edge_alignment |= fragment.edges

            if error_happened:
                logger.writeln('The %d-th sentence has wrong alignments' % curr_sent_index)
                logger.writeln(str(amr_graph))
                continue

            toks = sent.split()
            assert len(toks) == len(word_pos_seq)
            unaligned_toks = [(i, tok) for (i, tok) in enumerate(toks) if i not in aligned_words]
            (aligned, unaligned) = amr_graph.recall_unaligned_concepts(edge_alignment, unaligned_toks, lemma_map, stop_words, refine)

            for frag in aligned:
                aligned_fragments.append((frag.start, frag.end, frag))
                print '****Retrieved alignment****'
                print ' '.join(frag.str_list())
                print str(frag)
                print '********'

            if len(unaligned) > 10:
                print '******problematic sentence*****'
                print sents[sent_index]
                print amr_graphs[sent_index]
                print '****unaligned fragments****'
                for frag in unaligned:
                    print str(frag)
                    #print str(amr_graph.nodes[frag.root])
                print '********'
                continue

            aligned_fragments = sorted(aligned_fragments) #This sort this used to constrain the possible combinations

            curr_start = 0
            curr_end = 0
            tmp_w.write('%d\n' % sent_index)
            for start, end, frag in aligned_fragments:
                if start < curr_end:
                    error_happened = True
                    break

                #curr_hg = Hgraph()
                #var_mapping = {}
                #visited_index = set()
                #frag.init_new_hgraph(curr_hg, frag.root, var_mapping, visited_index, None, None, None,
                #frag_str =
                attrs = []
                new_node = FragmentHGNode(FRAGMENT_NT, start, end, frag)
                s = Sample(hypergraph.Hypergraph(new_node), 0)
                new_node.cut = 1
                new_rule, _ = s.extract_one_rule(new_node, None, list(frag.ext_set))

                frag_str = str(new_rule.e)
                attrs.append(frag_str)
                span_str = '%d %d' % (start, end)
                attrs.append(span_str)
                w_str = ' '.join(frag.str_list())
                attrs.append(w_str)
                p_str = ' '.join(pos_seq[start:end])
                attrs.append(p_str)
                tmp_w.write('%s\n' % '####'.join(attrs))

                curr_end = end

            if error_happened:
                tmp_w.write('\n')
                continue

            tmp_w.write('\n')
            continue
            aligned_fragments = [z for x,y,z in aligned_fragments]

            print >> unalign_f, 'sent %d:' % sent_index
            unaligned_words = set(range(len(toks))) - aligned_words
            un_seq = []
            for pos in unaligned_words:
                un_seq.append(toks[pos])
            print >> unalign_f, ' '.join(un_seq)

            hg = fragment_decomposition_forest(aligned_fragments, amr_graph, unaligned, edge_alignment, refine)
            if amr_graph.check_self_cycle():
                logger.writeln(str(amr_graph))
            elif hg is not None:
                frag_forests.append(hg)
                sent_indexes.append(curr_sent_index)
                used_sent_num += 1
            else:
                logger.writeln('Failed to build a goal for %d' % curr_sent_index)

            if len(unaligned) > 6:
                print '******problematic sentence*****'
                print sents[sent_index]
                print amr_graphs[sent_index]
                print '****unaligned fragments****'
                for frag in unaligned:
                    print str(frag)
                    #print str(amr_graph.nodes[frag.root])
                print '********'
                continue
            logger.writeln('Finished sentence %d.' % curr_sent_index)

        tmp_w.close()
        logger.writeln('The final used sentence num is %d' % used_sent_num)
        f4 = open(des_file, 'wb')
        cPickle.dump(frag_forests, f4)
        f4.close()
        f5 = open(used_sent_file, 'wb')
        cPickle.dump(sent_indexes, f5)
        f5.close()
        logger.writeln('finished')
        rule_f.close()
        #os.system('./run_dump')

if __name__ == '__main__':
    amr_file = sys.argv[1]
    sentence_file = sys.argv[2]
    alignment_file = sys.argv[3]
    #frag_forests = initialize_fragment_forest(amr_file, sentence_file, alignment_file)
    parallel_forest_construct(sys.argv)

#!/usr/bin/env python2.7
import sys
import amr
from amr import *
import amr_parser
from amr_parser import *
import amr_fragment
from amr_fragment import *
from collections import deque
#from fragment_forest import init_ext_frag

class AMREdge(object):
    def __init__(self, label, graph, h_node, t_node = None):
        self.label = label
        self.head = h_node
        self.tail = t_node
        self.is_coref = False #Denoting the tail node is a reentrance
        self.graph = graph

    def __str__(self):
        return self.label

    def set_coref(self, val):
        self.is_coref = val

    def isLeaf():
        return self.tail is None

'''There is some difference between leaf variables and the constant'''
class AMRNode(object):
    def __init__(self, graph, is_const = False):
        self.v_edges = []   #edges for relations
        self.c_edge = None #edge for node variable
        self.p_edges = []
        self.is_const = is_const
        self.use_quote = False
        self.has_reentrance = False
        self.graph = graph

    def is_var_node(self):
        return not self.is_const

    def set_quote(self, val):
        self.use_quote = val

    def set_reenter(self, val):
        self.has_reentrance = val

    def edge_set(self):
        for edge_index in self.p_edges:
            yield edge_index
        yield self.c_edge
        for edge_index in self.v_edges:
            yield edge_index

    def is_leaf(self):
        return len(self.v_edges) == 0

    def add_incoming(self, edge):
        self.v_edges.append(edge)

    def set_const_edge(self, edge):
        self.c_edge = edge

    def add_parent_edge(self, edge):
        self.p_edges.append(edge)

    def node_str(self):
        node_l = self.node_label()
        return self.graph.dict[node_l] if node_l in self.graph.dict else node_l

    def node_label(self):
        return self.graph.edges[self.c_edge].label

    def __str__(self):
        result = self.node_label()
        if self.is_var_node():
            result += ('/' + self.graph.dict[result])
        if self.use_quote:
            result = '"%s"' % result
        return result

class AMRGraph(object):
    def __init__(self, line):
        vars, var_values, rel_links = from_AMR_line(line)

        self.dict = {}
        label_to_node = {}
        self.nodes = []   #used to record all the nodes
        self.node_dict = {}
        self.edges = []
        self.edge_dict = {}

        #print 'get here'
        for i in range(len(vars)):
            curr_var = vars[i]
            self.dict[curr_var] = var_values[i] #maintain a dict for all the variable values

            #Setting a constant edge for a node: the variable name
            if curr_var not in label_to_node.keys(): #Haven't created a node for the current variable yet
                curr_node = AMRNode(self)
                self.node_dict[curr_node] = len(self.nodes)
                self.nodes.append(curr_node)
                curr_node_index = self.node_dict[curr_node]

                if i == 0:
                    self.root = curr_node_index
                const_edge = AMREdge(curr_var, self, curr_node_index)
                self.edge_dict[const_edge] = len(self.edges)
                self.edges.append(const_edge)
                curr_edge_index = self.edge_dict[const_edge]

                curr_node.set_const_edge(curr_edge_index) #The const edge is set immediately after the initialization of a node
                label_to_node[curr_var] = curr_node_index

        for i in range(len(vars)):
            curr_var = vars[i]

            curr_node_index = label_to_node[curr_var]
            curr_node = self.nodes[curr_node_index]
            if curr_var in rel_links:
                for rel, linked_val, is_var, is_coref in rel_links[curr_var]:
                    if is_var:
                        assert linked_val in label_to_node.keys(), 'Current coref variable %s is not a node yet' % linked_val
                        tail_node_index = label_to_node[linked_val] #Find the existed linked node index
                        edge = AMREdge(rel, self, curr_node_index, tail_node_index)
                        if is_coref:  #The node for this variable has already been generated
                            edge.set_coref(True)

                        self.edge_dict[edge] = len(self.edges)
                        self.edges.append(edge)
                        curr_edge_index = self.edge_dict[edge]

                        #curr_node = self.nodes[curr_node_index]
                        curr_node.add_incoming(curr_edge_index)

                        tail_node = self.nodes[tail_node_index]
                        tail_node.add_parent_edge(curr_edge_index)
                    else:
                        tail_node = AMRNode(self, True)  #Add a flag that it is a const node
                        if linked_val[0] == "\"" and linked_val[-1] == "\"":
                            linked_val = linked_val[1:-1]
                            tail_node.set_quote(True)

                        if '/' in linked_val:  #Dangerous here, pruned during the read amr procedure
                            try:
                                assert False, 'This should not happen again'
                            except:
                                print >> sys.stderr, linked_val
                                linked_val = linked_val.replace('/', '@@@@')
                                print >> sys.stderr, linked_val

                        self.node_dict[tail_node] = len(self.nodes)
                        self.nodes.append(tail_node)
                        tail_node_index = self.node_dict[tail_node]

                        #if '/' in linked_val:
                        #    linked_val = '///%s' % linked_val  #here is a special trick for identifying constant that has '/' in it
                        tail_const = AMREdge(linked_val, self, tail_node_index)

                        self.edge_dict[tail_const] = len(self.edges)
                        self.edges.append(tail_const)
                        tail_edge_index = self.edge_dict[tail_const]

                        tail_node.set_const_edge(tail_edge_index)
                        edge = AMREdge(rel, self, curr_node_index, tail_node_index)

                        self.edge_dict[edge] = len(self.edges)
                        self.edges.append(edge)
                        curr_edge_index = self.edge_dict[edge]

                        curr_node.add_incoming(curr_edge_index)
                        tail_node.add_parent_edge(curr_edge_index)

    def set_sentence(self, s):
        self.sent = s

    def print_info(self):
        #print all nodes info
        print 'Nodes information:'
        print 'Number of nodes: %s' % len(self.nodes)
        for node in self.node_dict.keys():
            print str(node), ',', self.node_dict[node]

        #print all edges info
        print 'Edges information'
        print 'Number of edges: %s' % len(self.edges)
        for edge in self.edge_dict.keys():
            print edge.label, ',', self.edge_dict[edge]

    def check_self_cycle(self):
        visited_nodes = set()
        sequence = []
        root_node = self.nodes[self.root]
        if len(root_node.p_edges) > 0:
            return True

        stack = [(self.root, 0)]
        while stack:
            curr_node_index, depth = stack.pop()
            if curr_node_index in visited_nodes:
                continue
            if depth >= len(sequence):
                sequence.append(curr_node_index)
            else:
                sequence[depth] = curr_node_index

            visited_nodes.add(curr_node_index)
            curr_node = self.nodes[curr_node_index]
            if len(curr_node.v_edges) > 0:
                for edge_index in reversed(curr_node.v_edges):  #depth first search
                    curr_edge = self.edges[edge_index]
                    child_index = curr_edge.tail
                    if child_index in sequence[:depth+1]:
                        return True
                    stack.append((child_index, depth+1))
        return False

    #Given a set of fragments, find the way they connect in the graph
    def derive_gld_rel(self, amr_fragments, frag_map, curr_alignment):
        stack = [(self.root, None, None, False)]
        root_to_frag = {}
        for frag in amr_fragments:
            assert frag.root not in root_to_frag
            root_to_frag[frag.root] = frag

        #tmp_edge_alignment =
        edge_alignment = bitarray(len(self.edges))
        if edge_alignment.count() != 0:
            edge_alignment ^= edge_alignment
        edge_alignment |= curr_alignment

        rels = []
        while len(stack) > 0:
            (curr_index, par_rel, par_id, is_approx) = stack.pop()

            curr_node = self.nodes[curr_index]
            #print curr_index, par_rel, par_id, is_approx

            if curr_index in root_to_frag:
                curr_frag = root_to_frag[curr_index]

                #edge_alignment |= curr_frag.edges
                curr_id = frag_map[curr_frag]
                if par_id is not None:
                    rels.append((par_id, par_rel, curr_id, is_approx))

                for edge_index in curr_node.v_edges:
                    if edge_alignment[edge_index] == 1:
                        continue
                    curr_edge = self.edges[edge_index]
                    edge_alignment[edge_index] = 1
                    stack.append((curr_edge.tail, str(curr_edge), curr_id, False))

                for other_index in curr_frag.ext_set:
                    if other_index == curr_index:
                        continue
                    curr_node = self.nodes[other_index]
                    for edge_index in curr_node.v_edges:
                        if edge_alignment[edge_index] == 1:
                            continue
                        edge_alignment[edge_index] = 1
                        curr_edge = self.edges[edge_index]
                        stack.append((curr_edge.tail, str(curr_edge), curr_id, True))

            else: #Found an unaligned node
                #edge_alignment[curr_node.c_edge] = 1
                for edge_index in curr_node.v_edges:
                    if edge_alignment[edge_index] == 1:
                        continue
                    curr_edge = self.edges[edge_index]
                    edge_alignment[edge_index] = 1
                    stack.append((curr_edge.tail, str(curr_edge), par_id, True))
        return rels



    #Return the order each node was visited
    def dfs(self):
        visited_nodes = set()
        sequence = []
        stack = [(self.root, None, 0, False)]
        while stack:
            curr_node_index, par_rel, depth, is_coref = stack.pop()
            sequence.append((curr_node_index, is_coref, par_rel, depth)) #push a tuple recording the depth search info
            if is_coref:
                continue

            visited_nodes.add(curr_node_index)
            curr_node = self.nodes[curr_node_index]
            if len(curr_node.v_edges) > 0:
                for edge_index in reversed(curr_node.v_edges):  #depth first search
                    curr_edge = self.edges[edge_index]
                    curr_rel = curr_edge.label
                    child_index = curr_edge.tail
                    if not curr_edge.is_coref:
                        stack.append((child_index, curr_rel, depth+1, False))
                    else:
                        stack.append((child_index, curr_rel, depth+1, True))
        return sequence

    #Try breadth-first search to find unaligned edges
    #For all edges coming out of the same node, add additional processing for args
    #When storing the fragments, be careful with the external nodes
    def extract_unaligned_fragments(self, edge_alignment):
        #print len(edge_alignment), edge_alignment.count(), edge_alignment
        n_nodes = len(self.nodes)
        n_edges = len(self.edges)

        unaligned_fragments = set()

        visited_nodes = set()
        stack = deque([self.root])
        while stack:
            curr_node_index = stack.popleft()
            #print str(curr_node)
            args_list = [] #Additional storage of args
            ops_list = [] #Additional storage of ops
            visited = curr_node_index in visited_nodes
            if visited:
                continue
            visited_nodes.add(curr_node_index)
            curr_node = self.nodes[curr_node_index]
            c_edge_index = curr_node.c_edge
            #c_edge_index = self.edge_dict[const_edge]
            #node_index = self.node_dict[curr_node]
            frag = None
            ext_set = set()

            external = False
            if len(curr_node.p_edges) > 0:
                external = True
            #Assume there is unaligned constant edge
            if edge_alignment[c_edge_index] == 0:
                frag = AMRFragment(n_edges, n_nodes, self)
                frag.set_node(curr_node_index)
                frag.set_edge(c_edge_index)
                frag.set_root(curr_node_index)
                #ext_set = set([curr_node_index])
                #ext_set.add(curr_node_index)
                #frag.set_ext_set(ext_set)
                #unaligned_fragments.add(frag)
            else:
                external = True

            if len(curr_node.v_edges) > 0:
                for curr_edge_index in curr_node.v_edges:
                    #edge_index = self.edge_dict[edge]
                    curr_edge = self.edges[curr_edge_index]
                    if edge_alignment[curr_edge_index] == 1: #This edge has already been aligned
                        external = True
                        if curr_edge.tail not in visited_nodes:
                            stack.append(curr_edge.tail)
                        continue
                    #curr_rel = curr_edge.label
                    tail_node_index = curr_edge.tail

                    args_list.append((curr_edge_index, tail_node_index))

                    if not tail_node_index in visited_nodes:
                        stack.append(tail_node_index)

            if external:
                ext_set.add(curr_node_index)

            if len(args_list) != 0:  #There are arguments going out of the current node
                if frag is None:
                    frag = AMRFragment(n_edges, n_nodes, self)
                    frag.set_root(curr_node_index)
                    #ext_set = set()
                #ext_set.add(curr_node_index)
                for rel_index, tail_index in args_list:
                    frag.set_edge(rel_index)
                    frag.set_node(tail_index)
                    ext_set.add(tail_index)

            if frag is not None:
                frag.set_ext_set(ext_set)
                unaligned_fragments.add(frag)

        return unaligned_fragments

    def recall_unaligned_concepts(self, edge_alignment, unaligned_toks, lemma_map, stop_words, refine=False):
        n_nodes = len(self.nodes)
        n_edges = len(self.edges)

        unaligned_fragments = set()
        aligned_fragments = set()

        visited_nodes = set()
        stack = deque([self.root])
        while stack:
            curr_node_index = stack.popleft()

            visited = curr_node_index in visited_nodes
            if visited:
                continue

            visited_nodes.add(curr_node_index)
            curr_node = self.nodes[curr_node_index]
            c_edge_index = curr_node.c_edge

            frag = None
            ext_set = set()

            external = False
            if len(curr_node.p_edges) > 0:
                external = True

            #Try to retrieve the concepts
            root_arcs = []

            is_pred = False #Use for deciding the category of the root node
            is_op = False

            if edge_alignment[c_edge_index] == 0:
                frag = AMRFragment(n_edges, n_nodes, self)
                frag.set_node(curr_node_index)
                frag.set_edge(c_edge_index)
                frag.set_root(curr_node_index)
                edge_alignment[c_edge_index] = 1

                concept_l = str(curr_node)
                concept_l = concept_label(concept_l)

                #Find the arguments structure for the current concept
                if len(curr_node.v_edges) > 0:
                    for curr_edge_index in curr_node.v_edges:

                        curr_edge = self.edges[curr_edge_index]
                        if edge_alignment[curr_edge_index] == 1: #This edge has already been aligned
                            external = True
                            if curr_edge.tail not in visited_nodes:
                                stack.append(curr_edge.tail)
                            continue

                        tail_node_index = curr_edge.tail
                        edge_label = str(curr_edge)
                        if is_root_arc(edge_label):
                            if 'ARG' in edge_label:
                                is_pred = True
                            else:
                                assert 'op' in edge_label
                                is_op = True

                            root_arcs.append((curr_edge_index, tail_node_index))
                        else:
                            external = True

                        #if not tail_node_index in visited_nodes:
                        #    stack.append(tail_node_index)

                if external:
                    ext_set.add(curr_node_index)

                if len(root_arcs) != 0 and False:  #There are arguments going out of the current node
                    for rel_index, tail_index in root_arcs:
                        edge_alignment[rel_index] = 1
                        frag.set_edge(rel_index)
                        frag.set_node(tail_index)
                        ext_set.add(tail_index)

                (pos, word) = match_word(concept_l, unaligned_toks, lemma_map, stop_words)
                frag.set_ext_set(ext_set)
                if refine:
                    init_ext_frag(frag, is_pred, is_op)

                if pos:
                    frag.set_span(pos, pos+1)
                    aligned_fragments.add(frag)
                else:
                    unaligned_fragments.add(frag)

            if len(curr_node.v_edges) > 0:
                for curr_edge_index in curr_node.v_edges:
                    curr_edge = self.edges[curr_edge_index]
                    if curr_edge.tail not in visited_nodes:
                        stack.append(curr_edge.tail)

        return (aligned_fragments, unaligned_fragments)

    #The concepts are all unary or split at the last stage
    def retrieve_fragment(self, integer_reps):
        n_nodes = len(self.nodes)
        n_edges = len(self.edges)

        frag = AMRFragment(n_edges, n_nodes, self)

        integer_concepts = sorted(integer_reps.split('+')) #newly added by xpeng
        integer_concepts = [x.split('.') for x in integer_concepts] #Each concept is a list of integers memorizing a path from the root to the concept
        integer_concepts = [[(int)(x) for x in y] for y in integer_concepts]

        n_concepts = len(integer_concepts)
        n_identified_concepts = 0

        #curr_node = self.root
        curr_num = 0
        lengths = [len(x) for x in integer_concepts]
        max_len = max(lengths)

        n_nodes = len(lengths)
        ext_set = set()
        curr_node_index = self.retrieve_first_concept(integer_concepts[0])
        if curr_node_index != self.root:
            ext_set.add(curr_node_index)

        #index = self.node_dict[curr_node]
        frag.set_root(curr_node_index) #The structure of the smallest grained fragment should be rooted structure
        #frag.set_node(index)
        curr_node = self.nodes[curr_node_index]

        c_edge_index = curr_node.c_edge
        #edge_index = self.edge_dict[edge]
        frag.set_edge(c_edge_index)

        if curr_node_index == self.root:
            #try:
            #    assert len(curr_node.p_edges) == 0, 'The root node has some parent nodes'
            #except AssertionError:
            #    #print str(self)
            #    #logger.writeln(str(self))
            #    #sys.exit(-1)
            #    return None

            length_two = len([None for x in lengths if (x == 2)])
            if len(curr_node.v_edges) > length_two:
                ext_set.add(curr_node_index)

        if n_nodes == 1:
            frag.set_ext_set(ext_set)
            return frag

        par_node_index = curr_node_index
        curr_node = None
        curr_depth = len(integer_concepts[0])

        #Starting from the node identified in the first step, retrieve the rest of the nodes and the relations that connect them
        for i in xrange(1, n_nodes):
            curr_len = lengths[i]
            try:
                assert len(integer_concepts[i]) == curr_depth+1, 'A relation has across more than 1 layer in just one step'
            except:
                print integer_concepts
                print integer_concepts[i]
                print curr_depth
                #sys.exit(-1)
                #logger.writeln(integer_concepts)
                return None
            (curr_edge_index, curr_node_index) = self.retrieve_one_concept(integer_concepts[i][curr_depth], par_node_index)
            if curr_edge_index == -1:
                return None
            #edge_index = self.edge_dict[curr_edge]
            frag.set_edge(curr_edge_index)
            #curr_node_index = self.node_dict[curr_node]
            curr_node = self.nodes[curr_node_index]
            if len(curr_node.p_edges) > 1:
                ext_set.add(curr_node_index)
            const_edge_index = curr_node.c_edge
            frag.set_node(curr_node_index)
            frag.set_edge(const_edge_index)

            if curr_len < max_len: #Not the leaf of the tree fragment
                if i < n_nodes-1 and curr_len == lengths[i+1]:
                    if len(curr_node.v_edges) > 0:
                        ext_set.add(curr_node_index)
                else:
                    n_next_len = len([None for x in lengths if (x == curr_len+1)])
                    if len(curr_node.v_edges) > n_next_len:
                        ext_set.add(curr_node_index)
                    par_node_index = curr_node_index
                    curr_depth = len(integer_concepts[i])
            else:
                if len(curr_node.v_edges) > 0:
                    ext_set.add(curr_node_index)

        frag.set_ext_set(ext_set)
        return frag

    def retrieve_first_concept(self, i_path):
        if len(i_path) == 1:
            assert i_path[0] == 0
            return self.root

        curr_node_index = self.root
        curr_node = self.nodes[curr_node_index]
        for curr_depth in xrange(1, len(i_path)):
            v_edges = curr_node.v_edges
            num = 0
            curr_index = i_path[curr_depth]
            for i in xrange(len(v_edges)): #search for non-coref nodes
                curr_edge_index = v_edges[i]
                curr_edge = self.edges[curr_edge_index]
                if not curr_edge.is_coref:
                    num += 1
                if num == curr_index+1:
                    curr_node_index = curr_edge.tail
                    curr_node = self.nodes[curr_node_index]
                    break
        return curr_node_index

    def retrieve_one_concept(self, child_num, par_node_index):
        par_node = self.nodes[par_node_index]
        v_edges = par_node.v_edges
        curr_index = 0
        curr_edge_index = -1
        curr_node_index = -1
        for i in xrange(len(v_edges)): #search for non-coref nodes
            curr_edge_index = v_edges[i]
            curr_edge = self.edges[curr_edge_index]
            if not curr_edge.is_coref:
                curr_index += 1
            if child_num+1 == curr_index:
                curr_node_index = curr_edge.tail
                break
        return (curr_edge_index, curr_node_index)

    #Do a depth-first traversal of the graph, print the amr format
    #Especially careful with the re-entrance structure
    def __str__(self):
        s = ""
        node_sequence = self.dfs()
        assert node_sequence[0][2] == None, 'The parent relation of the root should be none, %s' % node_sequence[0][2]

        dep_rec = 0 #record the current depth
        for curr_node_index, is_coref, par_rel, depth in node_sequence:
            curr_node = self.nodes[curr_node_index]
            #if len(curr_node.p_edges) > 0: #Current node has parents
            #    print str(curr_node) + ' :'
            #    for edge_index in curr_node.p_edges:
            #        print str(self.edges[edge_index])
            curr_c_edge = self.edges[curr_node.c_edge]
            curr_node_label = str(curr_c_edge)
            #if curr_node_label == 'c2':
            #    for edge in curr_node.p_edges:
            #        print self.edges[edge].label
            #    for edge in curr_node.v_edges:
            #        print self.edges[edge].label
            #    #print self.edges[curr_node.c_edge].label
            #if curr_node_label == 'c':
            #    for edge_index in curr_node.v_edges:
            #        print str(self.edges[edge_index])
            if par_rel == None:   #There is no relation going into this node
                if not is_coref and curr_node_label in self.dict.keys():
                    s += "(%s / %s" % (curr_node_label, self.dict[curr_node_label])
                else:
                    s += "(%s" % curr_node_label
            else:
                if depth < dep_rec:  #If the current layer is smaller than the current depth, then the previous few variables have finished traverse, print out the corresponding ) as finish
                    s += "%s" % ((dep_rec- depth) * ')')
                dep_rec = depth

                if curr_node.is_leaf():
                    if not is_coref and curr_node_label in self.dict.keys(): #In this case, current node is variable and is visited for the first time. Leaf variable
                        s += "\n%s:%s (%s / %s)"  % (depth*"\t", par_rel, curr_node_label, self.dict[curr_node_label])
                    else:
                        if curr_node_label not in self.dict.keys() and curr_node.use_quote:
                            s += '\n%s:%s "%s"' % (depth*"\t", par_rel, curr_node_label)
                        else:
                            s += "\n%s:%s %s" % (depth*"\t", par_rel, curr_node_label)
                else:
                    if not is_coref and curr_node_label in self.dict.keys(): #In this case, current node is variable and is visited for the first time. Not leaf variable
                        s += "\n%s:%s (%s / %s"  % (depth*"\t", par_rel, curr_node_label, self.dict[curr_node_label])
                    else:
                        s += "\n%s:%s %s" % (depth*"\t", par_rel, curr_node_label)
        if dep_rec != 0:
            s += "%s" % (dep_rec * ')')
        return s

    def print_variables(self):
        s = ''
        for node in self.nodes:
            var = str(node)
            s += str(node)
            s += '/'
            s += self.dict[var]
            s += ' '
        print s

#unaligned_words are formatted as tuples of (position, word)
def match_word(label, unaligned_words, lemma_map, stop_words):
    for (pos, word) in unaligned_words:
        if word.lower() in stop_words:
            continue
        lem_w = None
        if word in lemma_map:
            lem_w = list(lemma_map[word])[0]

        if len(label) > 4:
            if word[:3] == label[:3] or (lem_w and lem_w[:3] == label[:3]):
                return (pos, word)
        else:
            if word == label or (lem_w and lem_w == label):
                return (pos, word)
    return (None, None)

def concept_label(label):
    concept_label = label
    if '/' in concept_label:
        concept_label = concept_label.split('/')[1].strip()
    elif concept_label[0] == '"':
        assert concept_label[-1] == '"', 'weird constant %s' % concept_label
        concept_label = concept_label[1:-1]
    return concept_label

def is_root_arc(edge_label):
    return (edge_label[:3] == 'ARG' and 'of' not in edge_label) or edge_label[:2] == 'op'


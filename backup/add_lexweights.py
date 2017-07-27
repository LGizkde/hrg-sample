#!/usr/grads/lib/pypy-pypy-64bit/pypy/translator/goal/pypy-c

from __future__ import print_function
from __future__ import division

import sys
import time
import pickle
import os
import random
from math import log, exp, factorial, lgamma

#import syspath
import logprob
import alignment
from phrase_forest import make_rule, phrase_decomposition_forest
import logger
import gflags
from lexical_weighter import LexicalWeighter
from common import INF, ZERO
from levels import mark_level
from monitor import memory, resident
from rule import Rule

FLAGS = gflags.FLAGS

PHRASE_NT = 'A'

gflags.DEFINE_string(
    'base',
    'poisson',
    'Base distribution')
gflags.DEFINE_float(
    'alpha',
    5.0,
    'Concentration parameter in Dirichlet process.')
gflags.DEFINE_float(
    'discount',
    0.5,
    'Discount parameter in Pitman-yor process.')
gflags.DEFINE_integer(
    'maxcut',
    7,
    'no cut sampling when target or source span larger than maxcut.')
gflags.DEFINE_boolean(
    'sample_cut_only',
    False,
    'Sample at only nodes already cut (wrong implementation).')
gflags.DEFINE_integer(
    'sample_level',
    None,
    'Sample only nodes with level <= sample_level."')
gflags.DEFINE_integer(
    'level_inc',
    None,
    'Sample each level for #level_inc# iterations."')
gflags.DEFINE_boolean(
    'double_count',
    False,
    'Use double counting."')
gflags.DEFINE_boolean(
    'variable_alpha',
    False,
    'concentration parameter is different for each length')
gflags.DEFINE_boolean(
    'correct_edge_sampling',
    False,
    '"correct" path sampling. Number of incoming nodes is considered')
gflags.DEFINE_boolean(
    'lhs_conditional',
    False,
    'normalized over lhs.')
gflags.DEFINE_boolean(
    'seed_random',
    False,
    'seed random number generation.')
gflags.DEFINE_boolean(
    'sample_cut',
    True,
    'sample cut point (set to no to sample only minimal rules)')
gflags.DEFINE_boolean(
    'sample_edge',
    True,
    'sample edge switch')
gflags.DEFINE_integer(
    'splits',
    2,
    'two-way nt split by default, but you can change it.')
gflags.DEFINE_string(
    'model',
    'PY',
    'model. choose from DP|PY')
gflags.DEFINE_integer(
    'split_iter',
    10,
    'do symbol split every #split_iter iterations')
gflags.DEFINE_boolean(
    'split_at_iter_0',
    True,
    'allow split to happen at iter 0')
gflags.DEFINE_boolean(
    'type',
    False,
    'use type-based sampling')
gflags.DEFINE_boolean(
    'refine',
    False,
    'use symbol refinement')
gflags.DEFINE_boolean(
    'check_index',
    False,
    'Check cut index for errors.')

def timed(l):
    prev = time.time()
    for i, x in enumerate(l, 1):
        if i % FLAGS.interval == 0:
            logger.writeln('%s (%s/sec)' %
                           (i, FLAGS.interval/(time.time()-prev)))
            prev = time.time()
        yield x

def lncr(n, r):
    "log choose r from n"
    return lgamma(n+1) - lgamma(r+1) - lgamma(n-r+1)

def discrete(l):
    s = sum(l)
    if s == 0.0:
        length =  len(l)
        l = [1/length] * length
    else:
        l = [x/s for x in l]
    s = 0.0
    r = random.random()
    for i, x in enumerate(l):
        s += x
        if r <= s:
            return i

def rule_size(rule):
    return  len(rule.f) + len(rule.e) - 2 * rule.arity + rule.scope()

def geometric(rule):
    p = 0.99

    l = rule_size(rule)
    if l == 0:
        return 1.0
    else:
        return (1-p)**(l-1) * p

def poisson(rule):
    mean = 2.0
    l = rule_size(rule)

    denominator = log(factorial(l))
    numerator = l * log(mean) - mean

    result = exp(numerator - denominator)

    return result if result > ZERO else ZERO

def weibull(rule):
    shape = 2.9
    scale = 2.3

    l = rule_size(rule)

    result = ((shape/scale)
            * (pow(float(l)/scale, shape-1))
            * exp(-pow(float(l)/scale, shape)))

    return result if result > ZERO else ZERO

def uniform(rule):
    return 1.0

# emulate a Dirichlet distribution with three events
# used for the ab_test
def abtest_base(base):
    return 1.0/3

#def uniform_base(rule):
#    vocab_size_e = 34245 + 1
#    vocab_size_f = 32492 + 1
#    e_prob = pow(vocab_size_f, -float(len(rule.f)))
#    f_prob = pow(vocab_size_e, -float(len(rule.e)))
#
#    result = e_prob * f_prob
#
#    return result if result > ZERO else ZERO

def uniform_base(rule):
    l = rule_size(rule)

    result = pow(2.0, -l)
    return result if result > ZERO else ZERO

base = poisson
rule_size_prob = poisson

def choose_m(n, c1, c2, sampler):
    p = 0.0
    weights = [0] * (n+1)
    for i in xrange(n):
        for x in c1:
            p += logprob.elog(sampler.posterior(x))
            sampler.count(x)
    weights[n] = p
    for i in xrange(n):
        for x in c1:
            sampler.discount(x)
            p -= logprob.elog(sampler.posterior(x))
        for x in c2:
            p += logprob.elog(sampler.posterior(x))
            sampler.count(x)
        weights[n-i-1] = p
    # discount all c2's
    for i in xrange(n):
        for x in c2:
            sampler.discount(x)
    for i in xrange(n+1):
        weights[i] += lncr(n, i)
    weights = [logprob.eexp(w) for w in weights]
    #print(weights)
    return discrete(weights)

class Sampler(object):
    def count(self, x):
        pass

    def discount(self, x):
        pass

    def posterior(self, x):
        return 1.0

    def choice_posterior(self, c):
        result = 1.0
        for x in c:
            result *= self.posterior(x)
            if result < ZERO:
                return ZERO
        return result

class NPSampler(Sampler):
    "NP for non-paramatric"
    def __init__(self):
        self.counts = {}
        self.rule_size_counts = {}  # key=rule size, value=counts
        self.rule_size_tables = {}  # key=rule size, value=estimated tables
        self.n = 0

    def update_rule_size_tables(self):
        self.rule_size_tables = {}

        # should be fixed for python3
        for rule, counts in self.counts.iteritems():
            l = rule_size(rule)
            estimated_number_of_tables = pow(counts, FLAGS.discount)
            self.rule_size_tables[l] =  (self.rule_size_tables.get(l, 0.0)
                                            + estimated_number_of_tables)

#    def pitman_yor_posterior(self, c):
#        result = 1.0
#        for x in c:
#            l = rule_size(x)
#            n_r = self.counts.get(x, 0.0)
#            T_r = pow(n_r, FLAGS.discount)
#            T_e = pow(self.n, FLAGS.discount)
#
#            result *= ((n_r - T_r*FLAGS.discount + (T_e + FLAGS.alpha)*base(x))
#                        / (self.n + FLAGS.alpha))
#
#            if result < ZERO:
#                return ZERO
#
#        return result

    '''
    This is used to calculate the posterior prob of one rule
    Before this procedure, should discount the rule counts
    '''
    def pitman_yor_posterior_rule_size(self, x):
        alpha = FLAGS.alpha
        if FLAGS.variable_alpha == True:
            alpha = FLAGS.alpha * rule_size_prob(x)

        l = rule_size(x)
        n_r = self.counts.get(x, 0.0)
        n_l = self.rule_size_counts.get(l, 0.0)
        T_r = pow(n_r, FLAGS.discount)
        T_l = self.rule_size_tables.get(l, 0.0)

        return (((n_r - T_r*FLAGS.discount + (T_l*FLAGS.discount + alpha)*base(x))
                * (rule_size_prob(x)))
                / (n_l + alpha))

    # Used for generating likelihood graph
    # Single Dirichlet process, no rule size
    def simple_dirichlet_posterior(self, x):
        n_r = self.counts.get(x, 0.0)
        return ((n_r + FLAGS.alpha*base(x))
                / (self.n + FLAGS.alpha))

    def simple_dirichlet_posterior_for_choice(self, c):
        result = 1.0

        alpha = FLAGS.alpha
        for x in c:
            alpha = FLAGS.alpha
            n_r = self.counts.get(x, 0.0)

            result *= ((n_r + alpha*base(x))
                        / (self.n + alpha))

            self.count(x)

        # To make it absolutely correct
        for x in c:
            self.discount(x)

        return result

    '''
    xpeng: add one to the rule count of x and the total rule count
    if Pitman-Yor, add one to the rule size count of rule_size(x)
    '''
    def count(self, x):
        #print('count %s' % x)
        self.counts[x] = self.counts.get(x, 0) + 1
        self.n += 1

        if FLAGS.model == 'PY':
            l = rule_size(x)
            self.rule_size_counts[l] = self.rule_size_counts.get(l, 0) + 1

    '''
    xpeng: extract rule count of x and total rule count by 1
    if Pitman-Yor, then extract rule size count of rule_size(x) by 1
    there should be deletion if count is 0
    '''
    def discount(self, x):
        #print('discount %s' % x)
        if FLAGS.model == 'PY':
            l = rule_size(x)
            try:
                c = self.rule_size_counts[l]
                if c == 1:
                    del self.rule_size_counts[l]
                else:
                    self.rule_size_counts[l] = c - 1
            except KeyError:
                print('Warning: rule size %d not seen before in discounting' % l,
                      file=sys.stderr)

        try:
            c = self.counts[x]
            if c == 1:
                del self.counts[x]
            else:
                self.counts[x] = c - 1
            self.n -= 1
        except KeyError:
            print('Warning: %s not seen before in discounting' % x,
                  file=sys.stderr)

    def double_count(self, node):
        rule = make_composed_rule(node)
        #print('count %s' % rule)
        c = self.counts.get(rule, 0)
        if c == 0:
            for child in node.incoming[node.edge].tail:
                if not child.cut:
                    self.double_count(child)
        self.counts[rule] = c + 1
        self.n += 1
        l = rule_size(rule)
        self.rule_size_counts[l] = self.rule_size_counts.get(l, 0) + 1

    '''
    xpeng:it seems there is something wrong in this method
    '''
    def double_discount(self, node):
        rule = make_composed_rule(node)
        l = rule_size(x)   #xpeng:????????
        logger.writeln(x)
        logger.writeln(l)
        logger.writeln(rule)
        #print('discount %s' % rule)
        c = self.counts.get(rule)
        assert c is not None, '%s not seen before in discounting' % rule
        c -= 1
        if c == 0:
            del self.counts[rule]
            del self.rule_size_counts[l]
            for child in node.incoming[node.edge].tail:
                if not child.cut:
                    self.double_discount(child)
        else:
            self.counts[rule] = c
        self.n -= 1
        self.rule_size_counts[l] -= 1

    # Now this is really wrong...
    def rule_size_likelihood(self):
        result = 0.0
        for rule, count in self.counts.iteritems():
            l = rule_size(rule)
            f_n_l = float(self.rule_size_counts[l])
            l_prob = rule_size_prob(rule)
            result += count * logprob.elog(float(count) * l_prob / f_n_l)
        return result

    # Only works for dirichlet process, no rule size
    def dp_likelihood(self):
        n = 0.0
        counts = {}
        result = 0.0
        alpha = FLAGS.alpha
        for r, count in self.counts.iteritems():
            for _ in range(count):
                n_r = float(counts.get(r, 0.0))
                result += logprob.elog((n_r+alpha*base(r))/(n+alpha))
                counts[r] = counts.get(r, 0.0) + 1.0
                n += 1.0
        return result

    def nsamples(self):
        return self.n

    def ntypes(self):
        return len(self.counts)

class NTSampler(Sampler):
    "NT for nonterminal"
    def __init__(self):
        self.samplers = {}

    def count(self, x):
        self.samplers.setdefault(x.lhs, NPSampler()).count(x)

    def discount(self, x):
        self.samplers.setdefault(x.lhs, NPSampler()).discount(x)

    def posterior(self, x):
        return self.samplers.setdefault(x.lhs, NPSampler()).posterior(x)

    def update_rule_size_tables(self):
        for sampler in self.samplers.itervalues():
            sampler.update_rule_size_tables()

    def nsamples(self):
        return sum(s.nsamples() for s in self.samplers.itervalues())

    def ntypes(self):
        return sum(s.ntypes() for s in self.samplers.itervalues())

    def likelihood(self):
        return sum(s.likelihood() for s in self.samplers.itervalues())


def child_symbols(nt):
    "Return range of symbol indices given parent symbol index"
    return range(nt*FLAGS.splits+1, (nt+1)*FLAGS.splits+1)

def parent_symbol(nt):
    return (nt - 1)//2

def children(node):
    edge = node.incoming[node.edge]
    result = []
    for n in edge.tail:
        if n.cut:
            result.append(n)
        else:
            result.extend(children(n))
    return result

def cut_nodes_under(node):
    "pre-order, self not included"
    for child in children(node):
        yield child
        for c in cut_nodes_under(child):
            yield c

def nodes_turned_on_under(node):
    result = []
    queue = [node]
    while len(queue) > 0:
        curr = queue.pop(0)
        result.append(curr)
        for child in curr.incoming[curr.edge].tail:
            queue.append(child)
    return result

def nodes_in_fragment_under(node):
    "return all nodes that are in the fragment marked by cut points, including self, pre-order"
    yield node
    if not node.cut:
        for n in node.incoming[node.edge].tail:
            for n1 in nodes_in_fragment_under(n):
                yield n1

'''
xpeng: return the non-terminal # of the node
1.if is symbol-refined, also return the node.nt attribute
2.else just PHRASE_NT
'''
def get_nt(node):
    if FLAGS.refine:
        return '[%s-%s]' % (PHRASE_NT, node.nt)
    else:
        return '[%s]' % PHRASE_NT

class TreeFile():
    def __init__(self, filename):
        self.f = open(filename, 'w')
        self.i = 1

    def write(self, line):
        self.f.write(line)

    def dump(self, sample):
        self.f.write('# %s\n' % self.i)
        self.f.write('%s\n' % sample.tree_str())
        self.i += 1

    def close(self):
        self.f.close()

def dump_trees(samples, filename):
    logger.writeln('dump trees')
    treefile = TreeFile(filename)
    for s in timed(samples):
        # call this before dumping rules for each sample!
        LEXICAL_WEIGHTER.compute_lexical_weights(s.a)
        treefile.dump(s)
    treefile.close()

'''
xpeng: find the nearest cut-node ancester
'''
def cut_parent(node):
    #print(node)
    assert hasattr(node, 'parent')
    p = node.parent
    while not p.cut:
        p = p.parent
    return p

'''
xpeng: set the parent attribute for each node under the sub-tree rooted at the current node
the current tree, not the forest
'''
def set_parents_under(node):
    for child in node.incoming[node.edge].tail:
        child.parent = node
        set_parents_under(child)

'''
xpeng: the Sample class is used to record one sample of sentences(hypergraph-structured)
'''
class Sample():
    '''
    xpeng: initialization
    1.init the hypergraph and alignment attribute
    2.for each node, init the level, cut, chosen-edge and parent info
    '''
    def __init__(self, hg, a):
        self.hg = hg  # hypergraph
        self.a = a  # alignment

        mark_level(self.hg) # xpeng:calculate the level attribute for each node in the hypergraph, in topological order
        '''
        This is to iterate all the nodes, which does not specify one tree
        '''
        for node in self.hg.nodes:
            node.cut = 1 # xpeng: each node is initiated as a cut-site
            #if node.fj - node.fi <= FLAGS.sample_max_phrase_len:
            #    if random.random() <= 0.9:
            #        node.cut = True
            #    else:
            #        node.cut = False
            node.edge = random.randrange(len(node.incoming)) # xpeng: randomly choose incoming edge for each node
            # this is wrong. there may be multiple parents trying to claim a child.
            #for c in node.incoming[node.edge].tail:
            #    c.parent = node
            node.nt = 0
            node.pnt = None  # parent nt (before splitting)
            # sampled in a particular iteration
            node.sampled = False

        set_parents_under(self.hg.root)

    '''
    xpeng: make a direct rule with the node and it's tail nodes
    '''
    def make_one_level_rule(self, node):
        mychildren = node.incoming[node.edge].tail
        rule = make_rule([node.fi, node.fj, node.ei, node.ej],
                         [[c.fi, c.fj, c.ei, c.ej] for c in mychildren],
                         self.a.fwords,
                         self.a.ewords,
                         get_nt(node),
                         [get_nt(c) for c in mychildren])
        return rule

    '''
    xpeng: make a composed rule with the node and all its cut-node descendants
    '''
    def make_composed_rule(self, node):
        mychildren = children(node)
        rule = make_rule([node.fi, node.fj, node.ei, node.ej],
                         [[c.fi, c.fj, c.ei, c.ej] for c in mychildren],
                         self.a.fwords,
                         self.a.ewords,
                         get_nt(node),
                         [get_nt(c) for c in mychildren])
        return rule

    '''
    return all the composed rules in the sub-tree rooted at the node
    return value: (cut node, rule)
    '''
    def composed_rules_under(self, node):
        mychildren = children(node)
        rule = make_rule([node.fi, node.fj, node.ei, node.ej],
                         [[c.fi, c.fj, c.ei, c.ej] for c in mychildren],
                         self.a.fwords,
                         self.a.ewords,
                         get_nt(node),
                         [get_nt(c) for c in mychildren])
        yield node, rule
        for child in mychildren:
            for n, rule in self.composed_rules_under(child):
                yield n, rule

    def tree_str(self):
        return self.tree_str_helper(self.hg.root)

    def tree_str_helper(self, node, indent=0):
        result = ''
        rule = self.make_composed_rule(node)
        rule.feats = [1.0]
        rule.feats.extend(LEXICAL_WEIGHTER.score_rule(self.a, rule))
        result += ' '*indent + str(rule) + '\n'
        for child in children(node):
            result += self.tree_str_helper(child, indent + 4)
        return result

    def init_from_file(self, f):
        line = f.readline()
        sentnum = line[2:-1]
        line = f.readline()





    def __str__(self):
        return self.str_helper_expand(self.hg.root)

    def str_helper(self, node, indent=0):
        result = ''
        rule = self.make_composed_rule(node)
        result += ' '*indent + str(rule) + ' ' +  str(node) + '\n'
        for child in children(node):
            result += self.str_helper(child, indent + 4)
        return result

    def str_helper_expand(self, node, indent=0):
        result = ''
        rule = self.make_one_level_rule(node)
        result += ' '*indent + str(rule) + ' ' +  str(node) + ' ' + ('cut: %s' % node.cut) + '\n'
        for child in node.incoming[node.edge].tail:
            result += self.str_helper_expand(child, indent + 4)
        return result

if __name__ == '__main__':
    gflags.DEFINE_integer(
        'interval',
        5000,
        'Print stat every #interval# sentences.')
    gflags.DEFINE_integer(
        'iter',
        1,
        'Number of sampling iterations.')
    gflags.DEFINE_integer(
        'dump_iter',
        1,
        'Dump trees every #dump_iter# iterations.')
    gflags.DEFINE_string(
        'dump',
        'dump',
        'Dump directory.')

    try:
        argv = FLAGS(sys.argv)  # parse flags
    except gflags.FlagsError as e:
        print('%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS))
        sys.exit(1)

    ffilename = argv[1]
    efilename = argv[2]
    afilename = argv[3]
    rfilename = argv[4]

    ffile = open(ffilename)
    efile = open(efilename)
    afile = open(afilename)
    rfile = open(rfilename)
    # The following reads all the alignments from the three files, keep a list of alignments for each sentence
    alignments = alignment.Alignment.reader_pharaoh(ffile, efile, afile)


    LEXICAL_WEIGHTER = LexicalWeighter()

    os.system('rm -rf %s' % FLAGS.dump)
    os.mkdir(FLAGS.dump)
    logger.file = open(os.path.join(FLAGS.dump, 'log'), 'w')
    #flagfile = open(os.path.join(FLAGS.dump, 'run.flag'), 'w')
    #flagfile.write(FLAGS.FlagsIntoString())
    #flagfile.close()

    SAMPLER = None

    # initialization
    #samples = []
    logger.writeln()
    logger.writeln('read samples')
    treefile = TreeFile('dump-test')
    for i, a in enumerate(timed(alignments), 1):
        hg, a = phrase_decomposition_forest(a)
        #samples.append(Sample(hg, a))
        s = Sample(hg, a)
        LEXICAL_WEIGHTER.compute_lexical_weights(s.a)
        treefile.write('# %s\n' % i)
        line = rfile.readline()
        treefile.write(line)
        line = rfile.readline()
        while line.strip() != '':
            indent = len(line) - len(line.lstrip(' '))
            line = line[indent: -1]
            r = Rule()
            r.fromstr(line)
            r.feats = [1.0]
            r.feats.extend(LEXICAL_WEIGHTER.score_rule(a, r))
            treefile.write(' '* indent + str(r) + '\n')
            line = rfile.readline()

        #treefile.dump(s)
    treefile.close()
    ffile.close()
    efile.close()
    afile.close()

    #dump_trees(samples, os.path.join(FLAGS.dump, 'iter-0000'))

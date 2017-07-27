#!/usr/bin/env python3

import sys
import time
import pickle
import os
import random

import logprob
import alignment
from phrase_forest import make_rule, phrase_decomposition_forest
import logger
import gflags
from lexical_weighter import LexicalWeighter
from levels import mark_level
FLAGS = gflags.FLAGS

gflags.DEFINE_float(
    'alpha',
    0.001,
    'Concentration parameter in Dirichlet process.')
gflags.DEFINE_integer(
    'maxcut',
    10,
    'No cut when target or source span larger than maxcut.')
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

#random.seed(0)

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

def exponential(r):
    return 2**(len(r.f))

def geometric(rule):
    p = 0.99
    # l = len(rule.f) + len(rule.e)
    # l = len(rule.f)
    l = len(rule.f)
    if l == 0:
        return 1.0
    else:
        return (1-p)**(l-1) * p

base = geometric

class Sampler(object):
    def __init__(self):
        self.counts = {}
        self.n = 0

    def choice(self, choices, i):
        for x in choices[i]:
            self.discount(x)

        result = discrete([self.posterior(c) for c in choices])

        for x in choices[result]:
            self.count(x)

        return result

    def posterior(self, c):
        result = 1.0
        for x in c:
            result *= ((self.counts.get(x, 0.0) + FLAGS.alpha*base(x))/
                       (self.n + FLAGS.alpha))
        return result

    def count(self, x):
        self.counts[x] = self.counts.get(x, 0) + 1
        self.n += 1

    def discount(self, x):
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

    def likelihood(self):
        result = 0.0
        for x, c in self.counts.items():
            result += c*logprob.elog(float(c)/self.n)
        return result

sampler = Sampler()

def init_sample(hg):
    mark_level(hg)
    for node in hg.nodes:
        node.cut = True
        #if node.fj - node.fi <= FLAGS.sample_max_phrase_len:
        #    if random.random() <= 0.9:
        #        node.cut = True
        #    else:
        #        node.cut = False
        node.edge = random.randrange(len(node.incoming))
        # this is wrong. there may be multiple parents trying to claim a child.
        #for c in node.incoming[node.edge].tail:
        #    c.parent = node

def sample(hg):
    queue = [(hg.root, None)]
    while len(queue) > 0:
        node, parent = queue.pop(0)
        if FLAGS.sample_level is None or node.level <= FLAGS.sample_level:
            sample_edge(node, parent)
            sample_cut(node, parent)
        if node.cut:
            parent = node
        if FLAGS.sample_cut_only:
            for child in children(node):
                queue.append((child, parent))
        else:
            for child in node.incoming[node.edge].tail:
                queue.append((child, parent))

def posterior(x):
    return 1

def choice(l, i):
    result = discrete([posterior(x) for x in l])
    return result

def children(node):
    edge = node.incoming[node.edge]
    result = []
    for n in edge.tail:
        if n.cut:
            result.append(n)
        else:
            result.extend(children(n))
    return result

def sample_edge(node, parent):
    if parent is None:  # root node
        parent = node
    if len(node.incoming) == 1:  # trivial case
        return
    rule_lists = []
    old = node.edge
    for i in range(len(node.incoming)):
        node.edge = i
        rule_lists.append([r for r in composed_rules_under(parent)])
    i = sampler.choice(rule_lists, old)
    node.edge = i

def sample_cut(node, parent):
    if parent is None:  # root node
        return
    if node.fj - node.fi > FLAGS.maxcut or node.ej - node.ei > FLAGS.maxcut:
        return
    rule_lists = []
    old = node.cut

    node.cut = 0
    rule_lists.append([make_composed_rule(parent)])
    node.cut = 1
    rule_lists.append([make_composed_rule(parent), make_composed_rule(node)])

    i = sampler.choice(rule_lists, old)
    node.cut = i

def make_composed_rule(node):
    rule = make_rule([node.fi, node.fj, node.ei, node.ej],
                     [[c.fi, c.fj, c.ei, c.ej] for c in children(node)],
                     a.fwords,
                     a.ewords)
    return rule

def composed_rules_under(node):
    mychildren = children(node)
    rule = make_rule([node.fi, node.fj, node.ei, node.ej],
                     [[c.fi, c.fj, c.ei, c.ej] for c in mychildren],
                     a.fwords,
                     a.ewords)
    yield rule
    for child in mychildren:
        for rule in composed_rules_under(child):
            yield rule

def tree_str(node, indent=0):
    result = ''
    rule = make_composed_rule(node)
    rule.feats = [1.0]
    rule.feats.extend(lexical_weighter.score_rule(a, rule))
    result += ' '*indent + str(rule) + '\n'
    for child in children(node):
        result += tree_str(child, indent + 4)
    return result

def dump_tree(hg, filename):
    f = open(filename, 'w')
    tree = tree_str(hg.root)
    #print(tree)
    f.write(tree)
    f.close()

if __name__ == '__main__':
    gflags.DEFINE_integer(
        'interval',
        100,
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
    ffile = open(ffilename)
    efile = open(efilename)
    afile = open(afilename)
    alignments = alignment.Alignment.reader_pharaoh(ffile, efile, afile)

    lexical_weighter = LexicalWeighter()

    os.system('rm -rf %s' % FLAGS.dump)
    os.mkdir(FLAGS.dump)
    logger.file = open(os.path.join(FLAGS.dump, 'log'), 'w')

    starttime = time.time()

    # initialization
    samples = []
    dumpdir = os.path.join(FLAGS.dump, 'iter_0000')
    os.mkdir(dumpdir)
    logger.writeln()
    logger.writeln('initialization')
    for i, a in enumerate(alignments, 1):
        if i % FLAGS.interval == 0:
            logger.writeln('%s sentences at %s secs/sent' %
                           (i, (time.time()-starttime)/i))
            starttime = time.time()
        hg, a = phrase_decomposition_forest(a)
        lexical_weighter.compute_lexical_weights(a)
        init_sample(hg)
        samples.append((hg, a))
        for rule in composed_rules_under(hg.root):
            sampler.count(rule)
        filename = os.path.join(dumpdir, 'tree_%s' % str(i).rjust(8, '0'))
        dump_tree(hg, filename)
        #for rule in composed_rules_under(hg.root):
        #    print(rule)
    ffile.close()
    efile.close()
    afile.close()

    logger.writeln('%s rules, %s rule types, loglikelihood: %s' %
                   (sampler.n, len(sampler.counts), sampler.likelihood()))

    # sampling
    iteration = 1
    iter_start = time.time()
    if FLAGS.level_inc is not None:
        FLAGS.sample_level = 1
    while iteration <= FLAGS.iter:
        logger.writeln()
        logger.writeln('iteration %s' % iteration)
        if FLAGS.sample_level is not None:
            logger.writeln('level <= %s' % FLAGS.sample_level)

        if FLAGS.level_inc is not None and iteration % FLAGS.level_inc == 0:
                FLAGS.sample_level += 1
        if iteration % FLAGS.dump_iter == 0:
            dumpdir = os.path.join(FLAGS.dump,
                                   'iter_%s' % str(iteration).rjust(4, '0'))
            os.mkdir(dumpdir)
        for i, s in enumerate(samples, 1):
            if i % FLAGS.interval == 0:
                logger.writeln('%s sentences at %s secs/sent' %
                               (i, (time.time()-starttime)/i))
                starttime = time.time()
            hg, a = s
            lexical_weighter.compute_lexical_weights(a)
            sample(hg)
            if iteration % FLAGS.dump_iter == 0:
                filename = os.path.join(dumpdir,
                                        'tree_%s' % str(i).rjust(8, '0'))
                dump_tree(hg, filename)
        
        logger.writeln('iteration time: %s sec' % (time.time() - iter_start))
        logger.writeln('%s rules, %s rule types, loglikelihood: %s' %
                       (sampler.n, len(sampler.counts), sampler.likelihood()))
        iter_start = time.time()
        iteration += 1

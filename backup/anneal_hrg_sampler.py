#!/usr/bin/python
from __future__ import print_function
from __future__ import division

import sys
import time
import pickle
import os
import random
import cPickle
from math import log, exp, factorial, lgamma
import copy
import socket

#import syspath
import logprob
import alignment
#from phrase_forest import make_rule, phrase_decomposition_forest
import logger
import gflags
#from lexical_weighter import LexicalWeighter
from common import INF, ZERO
from levels import mark_level
from monitor import memory, resident
from fragment_forest import *
from collections import deque
from HRGSample import *

FLAGS = gflags.FLAGS

PHRASE_NT = 'X'

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
gflags.DEFINE_boolean(
    'anneal',
    False,
    'Use simulated annealing')
gflags.DEFINE_boolean(
    'cease_update',
    False,
    'Cease using a annealing factor')
gflags.DEFINE_integer(
    'annealing_intervals',
    1,
    'How many iterations per step')
gflags.DEFINE_float(
    'annealing_factor',
    0.1,
    'Annealing factor for simulated annealing')
gflags.DEFINE_float(
    'annealing_stepsize',
    0.1,
    'Annealing step size')
gflags.DEFINE_float(
    'annealing_ceiling',
    1.0,
    'Annealing upper bound')
gflags.DEFINE_list(
    'nodes',
    'node82,node83,node84,node85',
    'paralleled nodes')
gflags.DEFINE_integer(
    'port',
    12345,
    'port number')
gflags.DEFINE_string(
    'host',
    'cycle2.cs.rochester.edu',
    'host machine')
gflags.DEFINE_boolean(
    'slave',
    False,
    'Slave node')
gflags.DEFINE_string(
    'current',
    'cycle2',
    'current child node')
gflags.DEFINE_integer(
    'currid',
    0,
    'current child id')


#base = poisson
#rule_size_prob = poisson
def update_sampler(samples, g_sampler):
    for sample in timed(samples):
        for n, rule in sample.composed_rules_under(sample.hg.root):
            g_sampler.count(rule)

'''
xpeng: init the sampler with the current settings of the samples
rule, total rule and rule_size count. Also initiate the type-based indexer
'''
def init_split(samples, split=True):
    #global SAMPLER
    logger.writeln('initialization. split=%s' % split)
    SAMPLER = init_sampler()
    for sample in timed(samples):
        sample.set_sampler(SAMPLER)
        if split:
            for node in sample.hg.nodes:
                node.pnt = node.nt
                node.nt = random.choice(child_symbols(node.pnt))
        total_count = 0
        #print type(sample.hg.root)
        for n, rule in sample.composed_rules_under(sample.hg.root):
            total_count += 1
            SAMPLER.count(rule)
        logger.writeln('Total number of rules: %d' % total_count)
    #for rule, c in SAMPLER.counts.items():
    #    logger.writeln('Rule hash: %d, Rule counts: %d' % (hash(rule), c))
    return SAMPLER

def child_symbols(nt):
    "Return range of symbol indices given parent symbol index"
    return range(nt*FLAGS.splits+1, (nt+1)*FLAGS.splits+1)

def parent_symbol(nt):
    return (nt - 1)//2

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

def init_sampler():
    if FLAGS.lhs_conditional:
        sampler = NTSampler()
    else:
        sampler = NPSampler()

    if FLAGS.model == 'PY':
        NPSampler.likelihood = NPSampler.rule_size_likelihood
        NPSampler.posterior = NPSampler.pitman_yor_posterior_rule_size
    elif FLAGS.model == 'DP':
        NPSampler.likelihood = NPSampler.dp_likelihood
        NPSampler.posterior = NPSampler.simple_dirichlet_posterior
    else:
        assert False, 'unsupported model'
    return sampler

class TreeFile():
    def __init__(self, filename):
        self.f = open(filename, 'w')
        self.i = 1

    def dump(self, sample):
        self.f.write('# %s\n' % self.i)
        self.f.write('%s\n' % sample.tree_str())
        self.i += 1

    def close(self):
        self.f.close()

def dump_trees(samples, filename):
    file_handle = open(filename, 'w')
    for sample in samples:
        sample.dump_hrg_rules(file_handle)
    file_handle.close()

    #logger.writeln('dump trees')
    #treefile = TreeFile(filename)
    #for s in timed(samples):
    #    # call this before dumping rules for each sample!
    #    LEXICAL_WEIGHTER.compute_lexical_weights(s.a)
    #    treefile.dump(s)
    #treefile.close()

def choose_k(n):
    return random.randrange(n+1)

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

def index(node):
    p = node.parent
    n = node
    result = []
    while True:
        for i, c in enumerate(p.incoming[p.edge].tail):
            if c is n:
                result.append(i)
        if p.cut:
            break
        n = p
        p = p.parent
    result.reverse()
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
    gflags.DEFINE_string(
        'sample_file',
        'tempfile',
        'sample file name')
    gflags.DEFINE_string(
        'forest_dir',
        'forest_dir',
        'Forest file directory.')
    gflags.DEFINE_string(
        'prefix',
        'forest',
        'Separate forest file prefix.')
    gflags.DEFINE_string(
        'file_indexes',
        '0',
        'Different forest file indexes, separted by +.')
    try:
        argv = FLAGS(sys.argv)  # parse flags
    except gflags.FlagsError as e:
        print('%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS))
        sys.exit(1)


    f = open(FLAGS.sample_file, 'rb')
    samples = cPickle.load(f)
    f.close()
    os.system('rm -rf %s' % FLAGS.dump) #remove the existing directory
    os.mkdir(FLAGS.dump)
    likelihood_file = FLAGS.dump + '/likelihood_file'
    likelihood_handle = open(os.path.join(FLAGS.dump, 'likelihood_file'), 'w')
    logger.file = open(os.path.join(FLAGS.dump, 'log'), 'w')
    flagfile = open(os.path.join(FLAGS.dump, 'run.flag'), 'w')
    flagfile.write(FLAGS.FlagsIntoString())
    flagfile.close()

    dump_trees(samples, os.path.join(FLAGS.dump, 'iter-0000'))

    SAMPLER = init_split(samples, False)
    start_likelihood = SAMPLER.likelihood()
    likelihood_handle.write('%lf\n' % start_likelihood)
    likelihood_handle.flush()

    logger.writeln('%s rules, %s rule types, loglikelihood: %s' %
                   (SAMPLER.nsamples(), SAMPLER.ntypes(), start_likelihood))

    iteration = 1
    ascent = True
    FLAGS.cease_update = False
    logger.writeln('Starting annealing factor: %.6f at iteration 0' % FLAGS.annealing_factor)
    iter_start = time.time()

    while iteration <= FLAGS.iter:
        logger.writeln()
        logger.writeln('iteration %s' % iteration)

        if FLAGS.model == 'PY' and FLAGS.discount > 0:
            SAMPLER.update_rule_size_tables()

        for s in timed(samples):
            s.sample()

        if FLAGS.anneal and (not FLAGS.cease_update) and iteration % FLAGS.annealing_intervals == 0:
            if ascent:
                if (FLAGS.annealing_factor + FLAGS.annealing_stepsize > FLAGS.annealing_ceiling) or abs(FLAGS.annealing_factor + FLAGS.annealing_stepsize - FLAGS.annealing_ceiling) < 1e-4: #Bigger or 'equal' to ceiling
                    FLAGS.annealing_factor = FLAGS.annealing_ceiling
                    logger.writeln('Current annealing factor: %.6f at iteration %d' % (FLAGS.annealing_factor, iteration))
                    ascent = False
                else:
                    FLAGS.annealing_factor += FLAGS.annealing_stepsize
                    FLAGS.annealing_factor = (float)('%.2f' % FLAGS.annealing_factor)
                    logger.writeln('Current annealing factor: %.6f at iteration %d' % (FLAGS.annealing_factor, iteration))
            else:
                assert (FLAGS.annealing_factor > 1.0) or abs(FLAGS.annealing_factor - 1.0) < 1e-4, 'at the start of assent, the annealing factor is smaller than 1.0'
                if (FLAGS.annealing_factor - FLAGS.annealing_stepsize < 1.0) or abs(FLAGS.annealing_factor - FLAGS.annealing_stepsize - 1.0) < 1e-4: #descent to a number smaller or equal to 1.0
                    FLAGS.annealing_factor = 1.0
                    FLAGS.cease_update = True #Mark the end of annealing factor
                    logger.writeln('Current annealing factor: %.6f at iteration %d' % (FLAGS.annealing_factor, iteration))
                else:
                    FLAGS.annealing_factor -= FLAGS.annealing_stepsize
                    FLAGS.annealing_factor = (float)('%.2f' % FLAGS.annealing_factor)
                    logger.writeln('Current annealing factor: %.6f at iteration %d' % (FLAGS.annealing_factor, iteration))

        iter_likelihood = SAMPLER.likelihood()
        likelihood_handle.write(str(iter_likelihood)+ '\n')
        likelihood_handle.flush()
        logger.writeln('iteration time: %s sec' % (time.time() - iter_start))

        logger.writeln('%s rules, %s rule types, loglikelihood: %s' %
                       (SAMPLER.nsamples(), SAMPLER.ntypes(), iter_likelihood))
        logger.writeln('memory: %s' % memory())
        logger.writeln('resident memory: %s' % resident())

        if iteration % FLAGS.dump_iter == 0:
            dump_trees(samples,
                       os.path.join(FLAGS.dump, 'iter-%s' % str(iteration).rjust(4, '0')))

        iter_start = time.time()
        iteration += 1

    likelihood_handle.close()

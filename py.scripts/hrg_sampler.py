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
from filter_stop_words import filter_vars

FLAGS = gflags.FLAGS

PHRASE_NT = 'X'

file_h = None
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
    12355,
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
#gflags.DEFINE_boolean(
#    'href',
#    False,
#    'use heuristic symbol refinement')

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
    #file_handle = open(filename, 'w')
    new_file = open(filename+ '.r', 'w') #New
    for sample in samples:
        #sample.dump_hrg_rules(file_handle)
        rules = sample.extract_derivation() #New
        #match = re.search('(:[^A\s\n]+)|(:A[^0-9]+)|(:A[0-9]+[\s]+)', s)
        if not rules:
            continue
        for rule in rules:
            new_file.write('%s\n' % filter_vars(rule.dumped_format()))
        new_file.write('\n')

    new_file.close()
    #file_handle.close()

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
        'data',
        'data',
        'Data directory.')
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

    #The work done by the master program
    #1. initialize a global counter
    #2. communicate to the children to get ready
    #3. receive signal from children
    if not FLAGS.slave:
        os.system('rm -rf %s' % FLAGS.dump) #remove the existing directory
        os.mkdir(FLAGS.dump)
        likelihood_file = FLAGS.dump + '/likelihood_file'
        likelihood_handle = open(os.path.join(FLAGS.dump, 'likelihood_file'), 'w')
        #logger.file = open(FLAGS.dump, 'w')
        logger.file = open(os.path.join(FLAGS.dump, 'log'), 'w')
        file_ids = FLAGS.file_indexes.split('+')
        global_sampler = init_sampler()
        for id in file_ids:
            #hg_file = os.path.join(FLAGS.data, 'forest_%s' % id)
            ##hg_file = './small_file'
            #f = open(hg_file, 'rb')
            #hypergraphs = cPickle.load(f)
            #f.close()

            #sent_no_file = os.path.join(FLAGS.data, 'used_sent_%s' % id)
            #f = open(sent_no_file, 'rb')
            #sent_nos = cPickle.load(f)
            #f.close()

            #samples = []
            #for (hg, sent_num) in zip(hypergraphs, sent_nos)[:20]:
            #    samples.append(Sample(hg, sent_num))

            #f = open('small_sample', 'wb')
            #cPickle.dump(samples[:20], f)
            #f.close()
            sample_file = os.path.join(FLAGS.data, 'sample_%s' % id)
            f = open(sample_file, 'rb')
            samples = cPickle.load(f)
            f.close()
            #sample_file = os.path.join(FLAGS.data, 'sample_%s' % id)
            #f = open(sample_file, 'wb')
            #cPickle.dump(samples, f)
            #f.close()

            update_sampler(samples, global_sampler)
            samples[:] = []

        dump_f = open(os.path.join(FLAGS.dump, 'sampler_file_iter0'), 'wb')
        cPickle.dump(global_sampler, dump_f)
        dump_f.close()

        start_likelihood = global_sampler.likelihood()
        likelihood_handle.write(str(start_likelihood)+'\n')
        likelihood_handle.flush()

        logger.writeln('%s rules, %s rule types, loglikelihood: %s' %
                       (global_sampler.nsamples(), global_sampler.ntypes(), start_likelihood))

        server = socket.socket()
        host = socket.gethostname()
        port = FLAGS.port
        server.bind((host, port))
        server.listen(5)
        sys.setrecursionlimit(sys.getrecursionlimit() * 10)

        # Now start the cluster programs
        assert len(FLAGS.nodes) == len(file_ids), 'The length of cluster nodes and the # of files do not match'

        count = 0
        for curr_node in FLAGS.nodes:
            cmd = ' '.join(sys.argv)
            #curr_node_index = curr_node + '_'+ str(count)
            options = '--slave --current=%s --currid=%d --host=%s' % (curr_node, count, host)
            error_log_file = 'error_log%d' % count
            os.system(r'ssh %s "cd %s; nohup %s %s >& %s" &' % (curr_node, os.getcwd(), cmd, options, error_log_file))
            count += 1

        # sampling
        iteration = 1
        n_parallel_nodes = len(FLAGS.nodes)
        while iteration <= FLAGS.iter:
            iter_start = time.time()
            client_list = []
            for i in xrange(n_parallel_nodes):
                client, addr = server.accept()
                client_list.append(client) #collect all the clients to be notified
            logger.writeln('Finished the %d-th iteration on the clusters' % iteration)

            #The next step is to merge the sampler
            init_sampler_file = os.path.join(FLAGS.dump, 'sampler_file_iter%d' % (iteration-1))
            f = open(init_sampler_file, 'rb')
            Global_sampler = cPickle.load(f)
            f.close()
            logger.writeln('The likelihood from last iteration is %s' % Global_sampler.likelihood())

            for i in xrange(n_parallel_nodes):
                curr_sampler_file = os.path.join(FLAGS.dump, 'sampler_file%d' % i)
                #logger.writeln('Ready to read file %s' % curr_sampler_file)
                f = open(curr_sampler_file, 'rb')
                another_sampler = cPickle.load(f)
                f.close()
                Global_sampler.add(another_sampler)

            curr_likelihood = Global_sampler.likelihood()
            likelihood_handle.write(str(curr_likelihood)+'\n')
            likelihood_handle.flush()
            logger.writeln('%s rules, %s rule types, loglikelihood: %s' %
                        (Global_sampler.nsamples(), Global_sampler.ntypes(), curr_likelihood))

            dump_sampler_file = os.path.join(FLAGS.dump, 'sampler_file_iter%d' % iteration)
            f = open(dump_sampler_file, 'wb')
            cPickle.dump(Global_sampler, f)
            f.close()
            last_sampler_file = os.path.join(FLAGS.dump, 'sampler_file_iter%d' % (iteration-1))
            os.system('rm -rf %s' % last_sampler_file)
            iteration += 1
            logger.writeln('iteration time: %s sec' % (time.time() - iter_start))
            for i in xrange(0, n_parallel_nodes):
                client_list[i].send('okay')  #Whatever message, just to notify the child node to continue
                client_list[i].close()
        server.close()
    else:
        logger.file = open(os.path.join(FLAGS.dump, 'child%d.log' % FLAGS.currid), 'w')
        #file_h = open('discarded_rules')
        #LEXICAL_WEIGHTER = LexicalWeighter()
        #The samples should be loaded for only once
        #node = FLAGS.current
        sample_file = os.path.join(FLAGS.data, 'sample_%d' % FLAGS.currid)
        #sample_file = './small_sample'
        f = open(sample_file, 'rb')
        samples = cPickle.load(f) #Load a subset of the whole data
        #samples = samples[:3]
        f.close()
        sys.setrecursionlimit(sys.getrecursionlimit() * 10)

        #Dump the initial gramamr
        dump_trees(samples, os.path.join(FLAGS.dump, 'iter-0000-%d' % FLAGS.currid))

        SAMPLER = init_sampler()
        # sampling at each cluster
        iteration = 1

        while iteration <= FLAGS.iter:
            iter_start = time.time()
            init_sampler_file = os.path.join(FLAGS.dump, 'sampler_file_iter%d' % (iteration-1))
            f = open(init_sampler_file, 'rb')
            SAMPLER = cPickle.load(f)
            f.close()
            logger.writeln('The likelihood from last iteration is %s' % SAMPLER.likelihood())
            init_sampler = copy.deepcopy(SAMPLER) #Make a copy of the initial sampler
            logger.writeln()
            logger.writeln('iteration %s' % iteration)

            if FLAGS.model == 'PY' and FLAGS.discount > 0:
                SAMPLER.update_rule_size_tables()

            for s in timed(samples):
                s.set_sampler(SAMPLER)
                s.sample()

            logger.writeln('Top ten rule counts')
            for (rule, count) in SAMPLER.top_stats():
                logger.writeln(str(rule)+ ':%d' % count)
            logger.writeln('The likelihood for current iteration is %s' % SAMPLER.likelihood())

            SAMPLER.substract(init_sampler)
            #Save the update made by the current node
            curr_sampler_file = os.path.join(FLAGS.dump, 'sampler_file%d' % FLAGS.currid)
            os.system('rm -rf %s' % curr_sampler_file)
            f = open(curr_sampler_file, 'wb')
            cPickle.dump(SAMPLER, f)
            f.close()
            #logger.writeln('Having written file %s' % curr_sampler_file)

            logger.writeln('iteration time: %s sec' % (time.time() - iter_start))
            logger.writeln('%s rules, %s rule types' %
                        (SAMPLER.nsamples(), SAMPLER.ntypes()))
            logger.writeln('memory: %s' % memory())
            logger.writeln('resident memory: %s' % resident())

            if iteration % FLAGS.dump_iter == 0:
                dump_trees(samples,
                       os.path.join(FLAGS.dump, 'iter-%s-%d' % (str(iteration).rjust(4, '0'), FLAGS.currid)))

            #Report to the main process the sub task has been done
            s = socket.socket()
            host = FLAGS.host
            port = FLAGS.port
            s.connect((host, port))
            s.recv(1024)
            s.close()
            iteration += 1
        logger.writeln('Finished job')

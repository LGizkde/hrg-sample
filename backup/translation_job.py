from time import time
import socket
import re
import sys

from grammar import Grammar
from rule import Rule
from rule import get_num_edges
import gflags
FLAGS = gflags.FLAGS
from decode import Decoder
import hypergraph
from hypergraph import Hypergraph
import logger
from consensus_training import NgramCounter

gflags.DEFINE_boolean(
    'output_hypergraph',
    False,
    'Output hypergraph files.')
gflags.DEFINE_boolean(
    'reverse_kbest_feature_score',
    False,
    'Output negative feature costs in Kbest output. Turn this to true when \
    the weight tuner assumes the decoder maximizes (as opposed to minimize) \
    model scores.')
gflags.DEFINE_boolean(
    'preprocess',
    False,
    'Test set is marked with $number and $date.')

pattern = re.compile(r'(\$number|\$date)\s*\{\s*(.*?)\s*\}')

def preprocess(line):
    """Input: a sentence with $number and $dates tokens.
       Output: the sentence with brackets removed and a map for special
               symbols"""
    return pattern.sub(r'\1', line), pattern.findall(line)

class TranslationJob(object):
    def __init__(self,
                 id,
                 line,
                 grammars,
                 features):
        self.id = id
        #if FLAGS.preprocess:
        #    self.line, self.special = preprocess(line)
        #else:
        #    self.line = line

        # added by freesunshine, for incorporating oov words

        items = line.strip().split('|||')
        self.line = items[0].strip()

        self.oov_idx = []
        self.ner_items = []

        for i in range(1,len(items)):
            if items[i].find('####') == -1: # no '++' is found, unaligned words
                self.oov_idx = [int(x) for x in items[i].strip().split()]
            else: # ner rules
                for ne in items[i].strip().split('++'):
                    try:
                        span, shrg = ne.strip().split('####')
                    except:
                        print ne, ne.strip().split('####')
                        sys.exit(-1)
                    span = [int(x) for x in span.strip().split('-')]
                    self.ner_items.append((span, shrg.strip().replace('##', '|||')))

        #print self.oov_idx
        #print self.ner_items

        self.grammars = grammars
        self.features = features
        self.out = None   # output
        self.kbest = None
        self.time = None
        self.ref = None  # reference, class RefCounter

        self.suffix = str(id).rjust(5, '0')


    def run(self):
        # update per-sentence grammars, if there's any
        for g in self.grammars:
            g.update(self.id)

        self.flog = open('%s/%s_%s' % (FLAGS.run_dir,
                                  'log',
                                  self.suffix),
                    'w')
        if FLAGS.show_time:
            self.flog.write('running on %s\n\n' % socket.gethostname())
            self.flog.flush()

        fwords = self.line.strip().split()


        # added by freesunshine, build the local grammar for oov words for each sentence
        rules = []
        if self.oov_idx is not None and len(self.oov_idx) > 0:
            #oov_weight = 8.0
            oov_weight = 0.0001
            for idx in self.oov_idx:
                fw = fwords[idx]
                ew = "."
                rule_str = "[A0-0] ||| %s ||| %s ||| %lf %lf %lf" %(fw, ew, oov_weight, oov_weight, oov_weight)
                rr = Rule()
                rr.fromstr(rule_str)
                rules.append(rr)

        if self.ner_items is not None and len(self.ner_items) > 0:
            for item in self.ner_items:
                concept_weight = 10.0
                st = item[0][0]
                ed = item[0][1]
                fw = ' '.join(fwords[st:ed])
                #concept_weight *= pow((ed-st), 2)
                ew = item[1]
                value = int(ew[2])

                #Here is the feature for difference of nonterminal type
                #concept_weight /= pow(1.4, value)

                #Here is the feature for the favor of longer spans
                #concept_weight *= pow(2, ed-st)

                #Here is the feature for the number of edges
                #concept_weight /= pow(2.0, get_num_edges(ew))
                #print >>sys.stder, ew, concept_weight
                #rule_str = "[A1-1] ||| %s ||| %s ||| " % (fw, ew)
                rule_str = "%s ||| " % ew
                #weight = 5
                if fw == ';':
                    rule_str += "%lf %lf %lf" % (concept_weight, concept_weight, concept_weight)
                else:
                    rule_str += "%lf %lf %lf" % (concept_weight, concept_weight, concept_weight)
                rr = Rule()
                #print rule_str
                rr.fromstr(rule_str)
                rules.append(rr)

        #print '===== local_gr ====='
        #for r in rules:
        #    print r

        local_gr = None
        if len(rules) > 0:
          local_gr = Grammar(FLAGS.rule_bin_size)
          local_gr.build(rules, self.grammars[0].features)

        if FLAGS.preprocess:
            self.fidx2replacement = {}
            j = 0
            for i, token in enumerate(fwords):
                if token in ('$number', '$date'):
                    self.fidx2replacement[i] = self.special[j][1]
                    j += 1

        self.flog.write('[%s][%s words] %s\n' %
                   (self.id, len(fwords), self.line))

        decoder = Decoder(fwords,
                          self.grammars,
                          self.features,
                          local_gr)

        begin_time = time()
        if FLAGS.decoding_method == 'agenda':
            item = decoder.decode()
        elif FLAGS.decoding_method == 'cyk':
            item = decoder.decode_cyk()
        elif FLAGS.decoding_method == 'earley':
            item = decoder.decode_earley()
        else:
            assert False, '"%s" not valid decoding option' \
                    % FLAGS.decoding_method
        self.time = time() - begin_time

        if item is None:
            self.out = '[decoder failed to build a goal item]'
        else:
            ttt, succ = item
            item = ttt
            hg = Hypergraph(item)
            hg.set_semiring(hypergraph.SHORTEST_PATH)
            hg.set_functions(lambda x: x.cost, None, None)
            hg.topo_sort()
            self.kbest = hg.root.best_paths()
            #output_tokens = self.kbest[0].translation[:]

            #if FLAGS.preprocess:
            #    for i in range(len(output_tokens)):
            #        if output_tokens[i] in ('$number', '$date'):
            #            fidx = self.kbest[0].composed_rule.we2f[i]
            #            if fidx is not None:
            #                output_tokens[i] = self.fidx2replacement[fidx]

            # @freesunshine target side string output
            #self.out = ' '.join(output_tokens[FLAGS.lm_order-1:
            #                                  1-FLAGS.lm_order])

            self.flog.write('Decuction Tree:\n%s\n' % self.kbest[0].tree_str())
            #self.out = str(self.kbest[0].translation)
            #if succ:
            self.out = self.kbest[0].translation.to_amr_format()[0]
            #else:
            #    self.out = self.kbest[0].translation.toAMR()
            lines = [x.strip() for x in self.out.split('\n')]
            self.out = "".join(lines)

            self.hg = hg
            if FLAGS.output_hypergraph:
                self.write_hypergraph()

        self.flog.write('%s\n' % self.out)
        self.flog.write('\n')
        #if item is not None:
        #    self.flog.write(self.kbest[0].tree_str())
        #    self.flog.write('\n')
        #    self.flog.write(hg.stats())
        #    self.flog.write('\n')
        self.flog.write(decoder.agenda_stats())
        self.flog.write('\n')
        self.flog.write(decoder.chart.stats())
        self.flog.write('\n')
        for dotchart in decoder.dotcharts:
            self.flog.write(dotchart.stats())
            self.flog.write('\n')

        if FLAGS.show_time:
            timeline = '{:<35}{:>15.2f}\n'.format('[time]:', self.time)
            self.flog.write(timeline)
        self.write_output_file()
        if FLAGS.output_kbest:
            self.write_kbest_to_file()
        self.flog.close()

    def write_output_file(self):
        fout = open('%s/%s_%s' % (FLAGS.run_dir,
                                  FLAGS.output,
                                  self.suffix),
                    'w')
        fout.write('%s\n' % self.out)
        fout.write('\n') #Added by xpeng
        fout.close()

    def write_kbest_to_file(self):
        fkbest = open('%s/%s_%s' % (FLAGS.run_dir,
                                    FLAGS.kbest_output,
                                    self.suffix),
                      'w')
        if self.kbest is None:
            # TODO: generating an empty kbest list for a particular
            # weight may not be good for MERT, because MERT doesn't get
            # to know this set of weight doesn't work
            self.flog.write('failed to generate kbest\n')
            self.flog.write('\n')
        else:
            for path in self.kbest.iter_top(FLAGS.kbest_k):
                #out = ' '.join(x for x in path.translation[FLAGS.lm_order-1:
                #                                           1-FLAGS.lm_order])
                #out = path.translation.to_amr_format()
                out = "".join([x.strip() for x in path.translation.to_amr_format()[0].split('\n')])

                if FLAGS.reverse_kbest_feature_score:
                    fcosts = ' '.join(str(-x) for x in path.fcosts)
                else:
                    fcosts = ' '.join(str(x) for x in path.fcosts)

                # zmert requires 0-based indices for input sentences
                fkbest.write('%s\n' % ' ||| '.join([str(self.id - 1),
                                                    out,
                                                    fcosts]))
        fkbest.close()

    def write_hypergraph(self):
        #ngram_counter = NgramCounter(FLAGS.lm_order)
        #ngram_counter.mark_ngrams(self.hg)
        fname = '%s/%s_%s' % (FLAGS.run_dir,
                              'hg',
                              self.suffix)
        if FLAGS.reverse_kbest_feature_score:
            for edge in self.hg.edges():
                edge.fcosts = [-fcost for fcost in edge.fcosts]
        self.hg.serialize(fname)
        self.hg.show()
        if FLAGS.reverse_kbest_feature_score:
            for edge in self.hg.edges():
                edge.fcosts = [-fcost for fcost in edge.fcosts]

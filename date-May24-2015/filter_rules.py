#!/usr/bin/python
import sys
from rule import get_num_edges
from rule import retrieve_edges
def get_f_nums(f_str):
    fields = f_str.strip().split()
    n_nonterm = 0
    for i in xrange(len(fields)):
        if fields[i][0] == '[' and fields[i][-1] == ']':
            n_nonterm += 1
    if fields[-1].strip() == '.':
        #if len(fields)-n_nonterm-1 == 0:
        #    print f_str
        return (n_nonterm, len(fields)-n_nonterm-1)
    return (n_nonterm, len(fields)-n_nonterm)

def get_f_words(f_str):
    fields = f_str.strip().split()
    nonterms = []
    for i in xrange(len(fields)):
        if (not (fields[i][0] == '[' and fields[i][-1] == ']')) and fields[i] != '.' and fields[i] != '..':
            nonterms.append(fields[i].strip())
    return nonterms

def get_e_nums(e_str):
    edges = retrieve_edges(e_str)
    n_nonterm = 0
    for i in xrange(len(edges)):
        if edges[i][0] == 'A' and '$' in edges[i]:
            n_nonterm += 1
    return (n_nonterm, len(edges)-n_nonterm)

def all_stop_words(word_list, stop_words):
    number = 0
    for word in word_list:
        if word.lower() not in stop_words:
            #print word
            number += 1
            if number >= 2:
                return False
    return True

def main(argv):
    lines = []
    with open(argv[1], 'r') as f:
        lines = f.readlines()
    discarded_rules = open(argv[2], 'w')
    saved_rules = open(argv[3], 'w')
    stop_f = open(argv[4], 'r')
    stop_lines = stop_f.readlines()
    stop_words = set()
    for line in stop_lines:
        stop_words.add(line.strip())

    for line in lines:
        fields = [x.strip() for x in line.strip().split('|||')]
        (f_non, f_term) = get_f_nums(fields[1])
        (e_non, e_term) = get_e_nums(fields[2])

        f_words = get_f_words(fields[1])
        #print f_words

        if f_non == 0:
            discarded_rules.write(line)
            continue

        assert f_non == e_non, 'nonterm number does not match for %s:%d %d' % (line, f_non, e_non)

        #if (f_term == 0 and e_term != 0) or f_non >= 4 or f_non == 0:
        if (f_term == 0 and e_term >= 2) or f_non >= 5:
        #if ((f_term == 0 and e_term != 0) and f_non >= 4) or f_non >= 5:
            discarded_rules.write(line)
        elif e_term == 0:
            if not all_stop_words(f_words, stop_words):
                discarded_rules.write(line)
            else:
                saved_rules.write(line)
        else:
            saved_rules.write(line)

        #if e_term != 0:
        #    saved_rules.write(line)
        #if f_non >= 4:
        #    discarded_rules.write(line)
        #else:
        #    saved_rules.write(line)

if __name__ == '__main__':
    main(sys.argv)

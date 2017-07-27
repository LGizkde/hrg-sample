#!/usr/bin/python
import sys
from rule import get_num_edges
from rule import retrieve_edges
non_map_words = set(['was', 'were', 'am', 'is', 'are', 'be', 'a', 'an', 'the', ',', '.', '..', '...', ':', '(', ')', '@-@'])
def get_f_nums(f_str):
    fields = f_str.strip().split()
    n_nonterm = 0
    n_term = 0
    for i in xrange(len(fields)):
        if fields[i][0] == '[' and fields[i][-1] == ']':
            n_nonterm += 1
        elif fields[i].strip() not in non_map_words:
            n_term += 1

    return (n_nonterm, n_term)

def get_f_words(f_str):
    fields = f_str.strip().split()
    nonterms = []
    for i in xrange(len(fields)):
        if (not (fields[i][0] == '[' and fields[i][-1] == ']')) and fields[i] != '.' and fields[i] != '..':
            nonterms.append(fields[i].strip())
    return nonterms

def is_nonterm(tok):
    is_non = tok[0] == '[' and tok[-1] == ']'
    if is_non:
        assert tok[1] == 'A'
        return True
    return False

def to_lower(f_str):
    fields = f_str.strip().split()
    for (i, tok) in enumerate(fields):
        if not is_nonterm(tok):
            fields[i] = tok.strip().lower()
    return ' '.join(fields)

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
    with open(argv[1], 'r') as f:
        with open(argv[2], 'w') as wf:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                fields = line.split(' ||| ')
                fields[1] = to_lower(fields[1])

                if len(fields[2].split('_')) >= 2:
                    continue

                print>>wf, ' ||| '.join(fields)
            wf.close()
            f.close()

if __name__ == '__main__':
    main(sys.argv)

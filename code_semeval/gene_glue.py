#!/usr/bin/python
import sys
def combinations(prefix, left):
    if left == 0:
        return [prefix]
    res = []
    for i in xrange(2):
        res.extend(combinations(prefix+ '%d' % i, left - 1))
    return res

def suffixes(ntype):
    if ntype == 0:
        return ['1']
    return combinations('', ntype)

def hyperedge_rep(ntype, suffix, typ):
    edge_str = '. :A%d-%s$%d' % (ntype, suffix, typ)
    for i in xrange(ntype-1):
        edge_str += ' .*%d' % i
    return '(%s)' % edge_str

def pure_edge_rep(ntype, suffix):
    edge_str = '. :A%d-%s$' % (ntype, suffix)
    for i in xrange(ntype-1):
        edge_str += ' .*%d' % i
    return '(%s)' % edge_str

def main(argv):
    #print suffixes(7)
    numNonterm = int(argv[1])
    f = open(argv[2], 'w')
    f1 = open(argv[3], 'w')

    f1.write('[S] ||| [S] [X] ||| (.) ||| 1\n')
    f1.write('[S] ||| [X] ||| (.) ||| 1\n')

    for i in xrange(numNonterm+1):
        for nonterm_suffix in suffixes(i):
            lhs = 'A%d-%s' % (i, nonterm_suffix)
            rhs1_type1 = '[%s,0] [A0-0,1]' % lhs
            rhs1_type2 = '[A0-0,0] [%s,1]' % lhs
            rhs2_type1 = hyperedge_rep(i, nonterm_suffix, 0)
            rhs2_type2 = hyperedge_rep(i, nonterm_suffix, 1)
            f.write('[%s] ||| %s ||| %s ||| 1\n' % (lhs, rhs1_type1, rhs2_type1))
            f.write('[%s] ||| %s ||| %s ||| 1\n' % (lhs, rhs1_type2, rhs2_type2))
            f1.write('[X] ||| [%s] ||| %s ||| 1\n' % (lhs, pure_edge_rep(i, nonterm_suffix)))
    f.close()
    f1.close()

if __name__ == '__main__':
    main(sys.argv)

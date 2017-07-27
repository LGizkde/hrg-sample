#!/usr/bin/python
import sys
from rule import Rule
def add_unlexed_weight(orig_file, lexed_file):
    f = open(orig_file, 'r')
    #badrules = open('badrules', 'w')
    des_f = open(lexed_file, 'w')
    for i, line in enumerate(f):
        if line.strip() == '' or line.strip()[0] == '#':
            des_f.write(line)
            continue
        line = line.strip()
        try:
            rule = Rule()
            rule.fromstr(line)
        except AssertionError:
            print 'something wrong with: %s' % line
            continue
        des_f.write('%s\n' % rule.dumped_format())
    f.close()
    des_f.close()

if __name__ == '__main__':
    orig_file = sys.argv[1]
    lexed_file = sys.argv[2]
    add_unlexed_weight(orig_file, lexed_file)

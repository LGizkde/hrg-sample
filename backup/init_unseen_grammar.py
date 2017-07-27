#!/usr/bin/python
import sys
from extract_unseen import init_set
if __name__ == '__main__':
    f1 = open(sys.argv[1], 'r')
    unseen_set = init_set(f1)
    f2 = open(sys.argv[2], 'w')
    for unk_word in unseen_set:
        f2.write('[A1] ||| %s ||| (. :%s ) ||| 1.0 0.0 0.0\n' % (unk_word, unk_word))
    f2.close()
    f1.close()


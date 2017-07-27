#!/usr/bin/python
import sys
if __name__ == '__main__':
    snt_file = sys.argv[1]
    voc_file = sys.argv[2]
    snt_f = open(snt_file, 'r')
    snt = snt_f.readline()
    voc_dict = set()
    while snt:
        for word in snt.strip().split():
            voc_dict.add(word)
        snt = snt_f.readline()

    snt_f.close()
    voc_f = open(voc_file, 'w')
    for word in voc_dict:
        voc_f.write('%s\n' % word)
    voc_f.close()


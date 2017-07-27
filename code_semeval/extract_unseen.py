#!/usr/bin/python
import sys
def init_set(file_f):
    curr_s = set()
    line = file_f.readline()
    while line:
        curr_s.add(line.strip())
        line = file_f.readline()
    return curr_s

if __name__ == '__main__':
    f1 = open(sys.argv[1], 'r')
    f2 = open(sys.argv[2], 'r')
    dict_set = init_set(f1)
    new_set = init_set(f2)
    f3 = open(sys.argv[3], 'w')
    for word in new_set:
        if word not in dict_set:
            f3.write('%s\n' % word)
    f3.close()
    f2.close()
    f1.close()

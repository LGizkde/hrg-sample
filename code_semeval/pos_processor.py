#!/usr/bin/python
def readPOSs(pos_file):
    pos_seqs = []
    with open(pos_file, 'r') as f:
        for line in f:
            if line.strip() == '':
                continue
            word_tags = [('/'.join(word_tag.split('/')[:-1]), word_tag.split('/')[-1]) for word_tag in line.strip().split()]
            pos_seqs.append(word_tags)
        f.close()
    return pos_seqs

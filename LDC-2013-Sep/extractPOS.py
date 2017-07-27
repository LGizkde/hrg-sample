#!/usr/bin/python
import sys
import re
def get_all_word_tags(line):
    word_tags = []
    re_word_tag = re.compile('\([^\(\)]*\)')
    position = 0
    while position < len(line):
        match = re_word_tag.search(line, position)
        if not match:
            break
        parts = match.group(0)[1:-1].split()
        if parts[0] == '-RRB-':
            parts[0] = ')'
            parts[1] = ')'
        elif parts[0] == '-LRB-':
            parts[0] = '('
            parts[1] = '('
        word_tags.append(tuple(parts))
        position = match.end()
    return word_tags

def extract_pos(lines):
    seqs = []
    new_seq = []
    for line in lines:
        if line.strip() == '':
            seqs.append(new_seq)
            new_seq = []
            continue
        word_tags_inline = get_all_word_tags(line)
        new_seq.extend(word_tags_inline)
    return seqs

def extract_snt_pos(seqs, snts):
    word_tag_seqs = []
    new_seq = []
    for seq, snt in zip(seqs, snts):
        w_in_seq = [w for (p, w) in seq]
        #assert ''.join(snt.strip().split()) == ''.join(w_in_seq), '\n%s %s' % (snt, ' '.join(w_in_seq))
        #assert ''.join(snt.strip().split()) == ''.join(w_in_seq), '\n%s\n%s' % (''.join(snt.strip().split()), ''.join(w_in_seq))
        assert len(''.join(snt.strip().split())) == len(''.join(w_in_seq)), '\n%s\n%s' % (''.join(snt.strip().split()), ''.join(w_in_seq))
        words = snt.strip().split()
        i = 0
        curr_s = ''
        skip = False
        in_paren = False
        for (pos, word) in seq:
            if skip:
                skip = False
                #print 'current'
                continue

            if words[i] == ')':
                assert in_paren
                in_paren = False

            if len(words[i]) <= len(word):
                if in_paren:
                    new_seq.append((words[i], 'STOP'))
                else:
                    new_seq.append((words[i], pos))

            if words[i] == '(':
                assert not in_paren
                in_paren = True

            curr_s += words[i]
            i += 1
            if len(curr_s) == len(word):
                #assert curr_s.lower() == word.lower(), '%s %s: %s' % (curr_s, word, snt)
                curr_s = ''
                continue
            while curr_s != word:
                if len(curr_s) > len(word):
                    skip = True
                    #print curr_s, 'here'
                    new_seq.append((curr_s, 'STOP'))
                    break
                try:
                    assert i < len(words), 'inconsistency for sentence: %s' % snt
                except:
                    print curr_s
                    print word
                    print seq, snt
                    sys.exit(-1)

                if words[i] == ')':
                    assert in_paren
                    in_paren = False

                if in_paren:
                    new_seq.append((words[i], 'STOP'))
                else:
                    new_seq.append((words[i], pos))

                if words[i] == '(':
                    assert not in_paren
                    in_paren = True

                curr_s += words[i]
                i += 1
            curr_s = ''
        word_tag_seqs.append(new_seq)
        new_seq = []
    return word_tag_seqs


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as f:
        lines = f.readlines()
        seqs = extract_pos(lines)
        #print seqs
        #for seq in seqs:
        #    word_tags = ['%s/%s' % (word ,pos) for (pos, word) in seq]
        #    print ' '.join(word_tags)
        with open(sys.argv[2], 'r') as f2:
            snts = f2.readlines()
            word_tag_seqs = extract_snt_pos(seqs, snts)
            for seq in word_tag_seqs:
                word_tags = ['%s/%s' % (word ,pos) for (word, pos) in seq]
                print ' '.join(word_tags)
        f.close()

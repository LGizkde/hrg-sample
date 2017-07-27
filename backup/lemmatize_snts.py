#!/usr/bin/python
import sys
from extract_alignment import initialize_lemma
from pos_processor import readPOSs

#Use the lemma dict to find the lemma of the word
def get_lemma(word, pos, infl_lemma, trans_lemma):
    if pos[0] in 'NV' and word in infl_lemma: #Noun or verb
        return list(infl_lemma[word])[0]
    elif pos[0] == 'V' and word in trans_lemma: #Verb
        return list(trans_lemma[word])[0]
    return word

#To lemmatize a sentence
def lemma_snt(sentence, pos_seq, infl_lemma, trans_lemma):
    #assert len(sentence.strip().split()) == len(pos_seq)
    words = [word for (word, pos) in pos_seq]
    assert sentence.strip() == ' '.join(words), '\n%s : %s' % (sentence, ' '.join(words))
    poss = [pos for (word, pos) in pos_seq]
    return ' '.join([get_lemma(word, pos, infl_lemma, trans_lemma) for (word, pos) in zip(sentence.strip().split(), poss)])

def main(argv):
    orig_file = argv[1]
    des_file = argv[2]
    pos_file = argv[5]

    infl_lemma = initialize_lemma(argv[3])
    trans_lemma = initialize_lemma(argv[4])
    word_pos_seqs = readPOSs(pos_file)

    with open(orig_file, 'r') as f:
        with open(des_file, 'w') as wf:
            for (i, line) in enumerate(f):
                wf.write('%s\n' % lemma_snt(line, word_pos_seqs[i], infl_lemma, trans_lemma))

            wf.close()
            f.close()

if __name__ == '__main__':
    main(sys.argv)

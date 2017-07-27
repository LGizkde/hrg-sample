#!/usr/bin/python
def get_num_edges(rule_str):
    re_edgelabel = re.compile(':[^\s\)]*')
    position = 0
    mapping = []
    num_edges = 0
    while position < len(rule_str):
        match = re_edgelabel.search(rule_str, position)
        if not match:
            break
        num_edges += 1
        position = match.end()
    return num_edges

def retrieve_edges(
def reform_edge(s):
    re_edgelabel = re.compile(':[^\s\)]*')
    position = 0
    new_s = ''
    mapping = []
    while position < len(s):
        match = re_edgelabel.search(s, position)
        if not match:
            new_s = new_s + s[position : ]
            break
        new_s = new_s + s[position : match.start()]
        token = match.group(0)
        if token[0] == ':' and token[1] == '[' and token[-1] == ']': # :[A0,1]  or  :[S]
            sym, idx = token[2:-1].split(',')
            new_s = new_s + ':%s$%s' %(sym, idx)
            mapping.append(int(idx))
        elif token[0] == ':' and token.find('$') != -1:  # :A0$1  or  :S$
            sym, idx = token[1:].split('$')
            new_s = new_s + token
            if idx != '':
                mapping.append(int(idx))
            else:
                mapping.append(0)
        else:  # :want
            new_s = new_s + token
        position = match.end()
    #print new_s
    #print mapping
    return (new_s, mapping)

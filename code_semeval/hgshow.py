#!/usr/bin/env python3
import sys

from decode import Item, Deduction
from hypergraph import Deserializer, Node, Edge

class SNode(Node):
    def dot_label(self, detailed):
        return str(self.id)

class SEdge(Deduction):
    def dot_label(self):
        r = str(self.rule)
        return r[:r.rfind('|||')]

if __name__ == '__main__':
    deserializer = Deserializer(SNode, SEdge)
    hg = deserializer.deserialize(sys.argv[1])
    hg.show()

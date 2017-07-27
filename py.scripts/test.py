#!/usr/bin/python
import sys
from util import log
from os import system
from rule import get_num_edges
if __name__ == '__main__':
    #status = system(r'ssh node85 "cat /proc/meminfo; ls ./"')
    #if status == 0:
    #    print 'okay'
    num_edges = get_num_edges('(. :a/after  :op1 (. :n/now ) :quant (. :b/between  :mod (. :n/next ) :op1 .*0  :op2 (. :quant (. :3 ) :t/temporal-quantity  :unit (. :y/year ))))')
    print num_edges

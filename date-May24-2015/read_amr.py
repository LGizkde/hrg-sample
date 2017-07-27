#!/usr/bin/env python2.6
import sys
import os
import time
import random
import amr
from smatch import get_amr_line
if __name__ == "__main__":
    amr_file = sys.argv[1]
    f = open(amr_file, 'r')
    wf = open(sys.argv[2], 'w')

    while True:
        curr_amr = get_amr_line(f)
        if not curr_amr or curr_amr.strip() == "":
            break
        wf.write(curr_amr+ '\n')

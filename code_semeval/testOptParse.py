#!/usr/bin/env python
import sys, math
if __name__ == "__main__":
    import optparse
    optparser = optparse.OptionParser()
    optparser.add_option('-w', '--weight', nargs=2, dest='weightfiles', help="lexical weight tables")
    optparser.add_option('-W', '--words', nargs=2, dest='words', help="parallel text files (words)")
    optparser.add_option('-P', '--pharaoh', dest='pharaoh', action='store_true', default=False, help="input is Pharaoh-style alignment (requires -W)")
    optparser.add_option('-r', '--ratio', dest='ratiofile', help="likelihood ratio file")
    (opts,args) = optparser.parse_args()
    sys.stdout.write(str(opts.weightfiles) + '\n')
    sys.stdout.write(str(opts.words)+ '\n')
    sys.stdout.write(str(opts.ratiofile)+ '\n')
    sys.stdout.write(str(args[0])+ '\n')
    #sys.stdout.write("okay\n")

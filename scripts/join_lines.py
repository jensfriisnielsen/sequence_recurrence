#!/usr/bin/env python2.7

import sys
import fileinput

prev_f1 = None
elements = set()
for line in fileinput.input():
    sline = line.strip()
    fields = sline.split()
    f1 = fields[0]
    f2 = fields[1]
    if prev_f1 is None:
        prev_f1 = f1

    if f1 == prev_f1:
        elements.add(f2)
    else:
        sys.stdout.write( "%s" % (prev_f1) )
        for e in sorted(list(elements)):
            sys.stdout.write( "\t%s" % (e) )
        sys.stdout.write( "\n" )
        elements = set()
        elements.add(f2)
        
    prev_f1 = f1

sys.stdout.write( "%s" % (prev_f1) )
for e in sorted(list(elements)):
    sys.stdout.write( "\t%s" % (e) )
sys.stdout.write( "\n" )

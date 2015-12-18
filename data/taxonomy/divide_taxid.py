#!/usr/bin/env python2.7

import fileinput
import sys
import string

transtab = string.maketrans(" ","_")

#fn_division = '/home/projects/pr_46500/references/taxonomy/division.dmp'
#fn_nodes = '/home/projects/pr_46500/references/taxonomy/nodes.dmp'

fn_division = sys.argv[1]
fn_nodes = sys.argv[2]

fout_handles = []

divisions = {}
f_division = open(fn_division)
for line in f_division:
    sline = line.rstrip("\n")
    fields = map(lambda x: x.strip(), sline.split("|"))
    id = int(fields[0])
    abbr = fields[1]
    name = fields[2]
    divisions[id] = name
    foutname = name.translate(transtab) + ".taxid.lst"
    fout = open(foutname,'w')
    fout_handles.append(fout)
f_division.close()
    
f_nodes = open(fn_nodes)
for line in f_nodes:
    sline = line.rstrip("\n")
    fields = map(lambda x: x.strip(), sline.split("|"))
    taxid = int(fields[0])
    divid = int(fields[4])
    fout = fout_handles[divid]
    fout.write( "%d\n" % (taxid) )
f_nodes.close()

for fout in fout_handles:
    fout.close()

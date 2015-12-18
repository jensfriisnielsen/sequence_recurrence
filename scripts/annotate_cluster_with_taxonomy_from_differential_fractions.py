#!/usr/bin/env python2.7

import sys
import getopt
import signal
import fileinput
import re
from collections import Counter
import glob
import os
import time

def usage():
    msg = """
Usage:
    {cmd} [options] <pars> [files]

Description:
    Goal: Annotate clusters with taxonomy, on a given level, based on difference between two maximal fractions.
    Taxonomy shows species and also what organism type the species fall into (Virus, Phage, Bacteria, etc.)
    Input format
    # Cluster Feature

    CDHIT
    # Cluster Feature C+F+ C+F- C-F+ C-F- P-value Odds-ratio Log-Odds-ratio 

    TAXONOMY tax-level
    # Cluster "tax-level" Taxid name rank count

    Output format
    # Cluster Feature C+F+ C+F- C-F+ C-F- P-value Odds-ratio Log-Odds-ratio Species Diff-fraction

Options:
    -r,--rank <string>
        default: species

    -h,--help
        print description

    -v,--verbose
"""
    print msg.format(cmd=sys.argv[0])

def parse_getopts(inargs, options={}):
    try:
        opts, args = getopt.getopt(inargs, ":hr:v",["help","rank=","verbose"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-r", "--rank"):
            options['rank'] = a
        elif o in ("-v", "--verbose"):
            options['verbose'] = True
        else:
            assert False, "unhandled option %s" % o

    return options, args

def difference_frac(counter,options=None):
    total = sum(counter[e] for e in counter)
    taxid,cnt = counter.most_common(1)[0]
    max_frac = cnt / float(total)
    if len(counter) > 1:
        cnt = counter.most_common(2)[1][1]
        next_frac = cnt / float(total)
        diff_frac = max_frac - next_frac
    else:
        diff_frac = max_frac
    return (taxid,max_frac,diff_frac)
    
        

def main():
    # PARSE OPTIONS
    options            = {}
    options['rank'] = 'species'
    options['verbose'] = False
    options, args      = parse_getopts(sys.argv[1:],options)

    sys.stderr.write( "# Reading taxonomy lists" )
    time1 = time.time()
    # read taxonomy lists
    taxdbtypes = {}
    roottaxdb = 'data/taxonomy'
    taxtypes_fn_list = glob.glob("%s/*.lst" % (roottaxdb))
    taxtypes = [ os.path.basename(fn).split('.')[0] for fn in taxtypes_fn_list ]
    for i,fn in enumerate(taxtypes_fn_list):
        f = open(fn)
        for line in f:
            taxid = int(line.strip())
            taxtype = taxtypes[i]
            taxdbtypes[taxid] = taxtype
        f.close()
    time2 = time.time()
    sys.stderr.write( " - %ds\n" % (time2-time1) )


    # PATHS (its a mess..)
    rootdir = '.'
    taxdir = rootdir + '/results/taxonomy' + pars
    assocdir = rootdir + '/results/feature_associations'
    basename_assoc = "pval.txt"
    basename_tax = "taxonomy.txt"

    fn_assoc = "%s/%s" % (assocdir,basename_assoc)
    fn_tax = "%s/%s" % (taxdir,basename_tax)
    
    clusters = set()
    associations = set()
    re_cluster = re.compile('Cluster\s+(\d+)')
    for line in fileinput.input(args):
        # cluster feature
        sline = line.strip("\n")
        fields = sline.split("\t")
        cstr = fields[0]
        feat = fields[1]
        if 'Cluster' in cstr:
            m = re_cluster.search(cstr)
            c = int(m.group(1))
        else:
            c = int(cstr)
        associations.add((c,feat))
        clusters.add(c)
    max_c = max(clusters)

    cluster_tax = {}
    taxdb = {}
    prev_c = None
    c_taxs = Counter()
    f_tax = open(fn_tax)
    for line in f_tax:
        # Cluster "tax-level" Taxid name rank count
        sline = line.strip("\n")
        fields = sline.split("\t")
        cstr = fields[0]
        if 'Cluster' in cstr:
            m = re_cluster.search(cstr)
            c = int(m.group(1))
        else:
            c = int(cstr)

        if c > max_c:
            break

        if c not in clusters:
            continue

        taxid = int(fields[2])
        name = fields[3]
        rank = fields[4]
        count = int(fields[5])

        if rank != options['rank']:
            continue

        if taxid not in taxdb:
            taxdb[taxid] = name

        if prev_c is not None and prev_c != c:
            # new c - calculate and store
            cluster_tax[prev_c] = difference_frac(c_taxs,options) # (taxid,frac,diff_frac)
            c_taxs = Counter()

        c_taxs[taxid] += count
        prev_c = c
    cluster_tax[prev_c] = difference_frac(c_taxs,options) # (taxid,frac,diff_frac)
    f_tax.close()


    f_assoc = open(fn_assoc)
    for line in f_assoc:
        sline = line.strip("\n")
        fields = sline.split("\t")
        cstr = fields[0]
        feat = fields[1]
        if 'Cluster' in cstr:
            m = re_cluster.search(cstr)
            c = int(m.group(1))
        else:
            c = int(cstr)

        if c > max_c:
            break

        if not c in clusters:
            continue 

        if not (c,feat) in associations:
            continue

        try:
            taxid,frac,diff_frac = cluster_tax[c]
            taxtype = taxdbtypes[taxid]
            name = taxdb[taxid]
            sys.stdout.write( "%s\t%s\t%d\t%.1f\t%.1f\t%s\n" % (sline,taxtype,taxid,frac,diff_frac,name) )
        except KeyError:
            # no taxonomy for this cluster
            sys.stdout.write( "%s\tNA\tNA\tNA\tNA\tNA\n" % (sline) )
    

    f_assoc.close()


if __name__ == "__main__":
    main()

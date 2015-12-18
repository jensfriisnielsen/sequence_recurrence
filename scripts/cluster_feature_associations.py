#!/usr/bin/env python2.7

import sys
import getopt
import signal
import fileinput
import os

def usage():
    msg = """
Usage:
    {cmd} [options] <files>

Description:
    print each cluster as a row of contigs

Options:
    --all-datasets <file>
        containing the 'baseline' of all datasets that is used for counting the different cells in the association matrix

    -h,--help
        print description

    -v,--verbose
"""
    print msg.format(cmd=sys.argv[0])

def parse_getopts(inargs, options={}):
    try:
        opts, args = getopt.getopt(inargs, ":hv",["all-datasets=","help","verbose"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("--all-datasets"):
            options['all-datasets'] = a
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-v", "--verbose"):
            options['verbose'] = True
        else:
            assert False, "unhandled option %s" % o

    return options, args

def main():
    # PARSE OPTIONS
    options            = {}
    options['all-datasets'] = None
    options['verbose'] = False
    options, args      = parse_getopts(sys.argv[1:],options)

    if options['all-datasets'] is None or not os.path.exists(options['all-datasets']):
        raise Exception("options --all-datasets must be a file")

    alldatasets = set()
    f = open(options['all-datasets'])
    for line in f:
        dataset = line.strip()
        alldatasets.add(dataset)
    f.close()

    featureset = set()
    featurepositive = {}
    featurenegative = {}
    for fn in args:
        fi = fn.split('/')[-1].split('.')[0] # features/f001.+
        sign = fn.split('/')[-1].split('.')[1]
        fset = set()
        f = open(fn)
        for line in f:
            dataset = line.strip()
            fset.add(dataset)
        f.close()
        if sign == '+':
            featurepositive[fi] = fset
        elif sign == '-':
            featurenegative[fi] = fset
        featureset.add(fi)

    features = sorted(list(featureset))
    for line in sys.stdin:
        sline = line.strip("\n")
        fields = sline.split("\t")

        cluster = fields[0]
        cp = set()
        for c in fields[1:]:
            cp.add(c)
        cn = alldatasets.difference(cp)

        for fi in features:
            fp = featurepositive[fi]
            fn = featurenegative[fi]
            cpfp = cp.intersection(fp)
            cpfn = cp.intersection(fn)
            cnfp = cn.intersection(fp)
            cnfn = cn.intersection(fn)
            sys.stdout.write( "%s\t%s\t%d\t%d\t%d\t%d\n" % (cluster,fi,len(cpfp),len(cpfn),len(cnfp),len(cnfn)) )


if __name__ == "__main__":
    main()

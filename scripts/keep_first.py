#!/usr/bin/env python2.7

import sys
import getopt
import signal
import fileinput

def usage():
    msg = """
Usage:
    {cmd} [options] <files>

Description:
    Iterates over <files>

Options:
    -d,--delimiter <char>

    -f,--fields <comma-separated-numbers>
        1-based

    -h,--help
        print description

    -v,--verbose
"""
    print msg.format(cmd=sys.argv[0])

def parse_getopts(inargs, options={}):
    try:
        opts, args = getopt.getopt(inargs, ":d:f:hv",["delimiter=","fields=","help","verbose"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for o, a in opts:
        if   o in ("-d", "--delimiter"):
            options['delimiter'] = a
        elif o in ("-f", "--fields"):
            options['fields'] = map(lambda x: int(x) - 1, a.split(','))
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
    options['delimiter'] = "\t"
    options['verbose'] = False
    options['fields']   = None
    options, args      = parse_getopts(sys.argv[1:],options)

    # ITERATE OVER FILES
    seen = set()
    for line in fileinput.input(args):
        if not options['fields'] is None:
            fields = line.split(options['delimiter'])
            s = fields[options['fields'][0]]
            if len(options['fields']) > 1:
                for i in options['fields'][1:]:
                    s += options['delimiter'] + "%s" % (fields[i])
        else:
            s = line
                
        
        if s not in seen:
            sys.stdout.write(line)
            seen.add(s)


if __name__ == "__main__":
    main()


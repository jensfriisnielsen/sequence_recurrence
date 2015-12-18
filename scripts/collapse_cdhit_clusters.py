#!/usr/bin/env python2.7

import sys
import getopt
import signal
import fileinput
#sys.path.append('/home/people/jef/local/lib/python2.7/jefutils')
sys.path.append('/home/projects/pr_46500/lib/python2.7/jefutils')
import jefutils

def usage():
    msg = """
Usage:
    {cmd} [options] <files>

Description:
    print each cluster as a row of contigs

Options:
    -h,--help
        print description

    --datasets
        print only the dataset part of each contig

    --circos
        NOT IMPLEMENTED YET

    -v,--verbose
"""
    print msg.format(cmd=sys.argv[0])

def parse_getopts(inargs, options={}):
    try:
        opts, args = getopt.getopt(inargs, ":hv",["help","datasets","verbose"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("--datasets"):
            options['datasets'] = True
        elif o in ("-v", "--verbose"):
            options['verbose'] = True
        else:
            assert False, "unhandled option %s" % o

    return options, args


class Alignment():
    """2x SeqPosition with IDs"""


    def __init__(self,qpos,spos,**kwargs):
        """assumes qpos and spos are SeqPosition()"""
        self.qpos = qpos
        self.spos = spos
        self.dict = kwargs


    def separation(self,other):
        return self.spos.intervals[0].separation(other.spos.intervals[0])


    def properly_paired(self,other,dist=1000):
        if self.sid == other.sid:
            return self.spos.intervals[0].properly_paired(other.spos.intervals[0])
        else:
            return False



def clusterelement(line):
    sline = line.strip()
    if '*' in sline:
        return None

    elif len(sline) > 0:
        fields = sline.split()
        i_a = int(fields[0])
        len_q = int(fields[1][:-3])
        qseqid = fields[2][1:-3]
        desc = fields[4]
        qstart = desc.split(':')[0]
        qend = desc.split(':')[1]
        sstart = desc.split(':')[2]
        desc2 = desc.split(':')[3]
        send = desc2.split('/')[0]
        qpos = seqposition.SeqPosition( (qstart,qend) )
        qpos.id = qseqid
        spos = seqposition.SeqPosition( (sstart,send) )
        pident = desc2[2][:-1]
        lenqp = len(qpos)
        lensp = len(spos)
        length = lenqp if lenqp < lensp else lensp
        a = Alignment(qpos,spos,pident=pident,length=length)
        return a


def clusterrepresentative(line):
    sline = line.strip()
    fields = sline.split()
    if '*' in sline:
        len_q = int(fields[1][:-3])
        sseqid = fields[2][1:-3]
        spos = seqposition.SeqPosition( (1,len_q) )
        spos.id = sseqid
        return spos

    elif len(sline) > 0:
        return None


class CdhitIterator():
    def __init__(self,handle):
        self.handle = handle
        self.i = 0

    def __iter__(self):
        return self

    def next(self):
        line = self.handle.readline()
        sline = line.strip()
        if line == '':
            raise StopIteration

        if self.i == 0:
            # first cluster
            line = self.handle.readline()
            sline = line.strip()

        c = Cluster()
        while sline[0] != '>':
            a = clusterelement(line)
            if a is not None:
                c.alignments.append(a)
            else:
                r = clusterrepresentative(line)
                c.representative = r
            line = self.handle.readline()
            sline = line.strip()
            if line == '':
                break

        for a in c.alignments:
            a.spos.id = c.representative

        self.i += 1
        return c
        

def parse_extended_format(handle):
    """http://weizhong-lab.ucsd.edu/cd-hit/wiki/doku.php?id=cd-hit_user_guide
    >Cluster 0
    0       726nt, >s0007_c000329... at 726:1:261597:262322/-/99.17%
    1       606nt, >s0007_c000437... at 606:1:333884:334489/-/99.01%
    2       515nt, >s0007_c000569... at 1:515:306103:306617/+/99.22%
    3       459nt, >s0007_c000666... at 459:1:178319:178777/-/99.78%
    4       452nt, >s0007_c000677... at 452:1:100805:101256/-/99.78%
    5       429nt, >s0007_c000725... at 429:1:94053:94481/-/99.30%
    """
    
    cdhit = CDHIT()
    for line in handle:
        sline = line.strip()
        fields = sline.split()

        if sline[0] == '>':
            c = Cluster()
            c.id = sline[1:]
            cdhit.clusters.append(c)
            if len(cdhit.clusters) > 1:
                c2 = cdhit.clusters[-2]
                for a in c2.alignments:
                    a.spos.id = c2.representative.id

        elif '*' in sline:
            # representative
            c = cdhit.clusters[-1]
            i_a = int(fields[0])
            len_q = int(fields[1][:-3])
            qseqid = fields[2][1:-3]
            spos = seqposition.SeqPosition( (1,len_q) )
            spos.id = qseqid
            c.representative = spos

        elif len(sline) > 0:
            i_a = int(fields[0])
            len_q = int(fields[1][:-3])
            qseqid = fields[2][1:-3]
            desc = fields[4]
            qstart = desc.split(':')[0]
            qend = desc.split(':')[1]
            sstart = desc.split(':')[2]
            desc2 = desc.split(':')[3]
            send = desc2.split('/')[0]
            qpos = seqposition.SeqPosition( (qstart,qend) )
            qpos.id = qseqid
            spos = seqposition.SeqPosition( (sstart,send) )
            pident = desc2[2][:-1]
            lenqp = len(qpos)
            lensp = len(spos)
            length = lenqp if lenqp < lensp else lensp
            a = alignment.Alignment(qpos,spos,pident=pident,length=length)
            c = cdhit.clusters[-1]
            c.alignments.append(a)

    return cdhit


class CDHIT():
    def __init__(self,clusters=[]):
        self.clusters = clusters
        

class Cluster():
    def __init__(self,alignments=None,representative=None,lengths=None,id=None):
        if alignments is None:
            alignments = []
        if lengths is None:
            lengths = []
        self.alignments = alignments
        self.lengths = lengths
        self.representative = representative
        self.id = id


def main():
    # PARSE OPTIONS
    options            = {}
    options['datasets'] = False
    options['verbose'] = False
    options, args      = parse_getopts(sys.argv[1:],options)

    f_cdhit = open(args[0])
#    cdhit = jefutils.cdhit.parse_extended_format(f_cdhit)

    i = -1
    for c in jefutils.cdhit.CdhitIterator(f_cdhit):
    #for i,c in enumerate(cdhit.clusters):
        i += 1
        id = c.representative.id
        if options['datasets']:
            id = id[:5]
            
        sys.stdout.write( "Cluster %d\t%s" % (i,id) )

        if not options['datasets']:
            for a in c.alignments:
                sys.stdout.write( "\t%s" % (a.qpos.id) )
        else:
            datasets = set()
            for a in c.alignments:
                qid = a.qpos.id[:5]
                datasets.add(qid)
            for ds in sorted(list(datasets)):
                sys.stdout.write( "\t%s" % (ds) )

        sys.stdout.write( "\n" )

    f_cdhit.close


if __name__ == "__main__":
    signal.signal(signal.SIGINT, jefutils.misc.quit_gracefully)
    main()

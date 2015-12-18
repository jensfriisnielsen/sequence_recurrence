#!/usr/bin/env python2.7

import sys
import getopt
import signal
import fileinput
import os
from string import maketrans
from collections import defaultdict
from collections import Counter
import time
import re

def usage():
    msg = """
Usage:
    {cmd} [options] <files>
     | {cmd} [options] 

Description:
    Iterates over stdin and or files.
    Downloads to default database if not specified.
    If already existing in database copies instead of downloads.

Options:
    --missing <outfilename>
        where to write the accessions that can't be found

    --skip <skipfilename>
        skip the entries in this file
"""
    print msg.format(cmd=sys.argv[0])

def parse_getopts(inargs, options={}):
    try:
        opts, args = getopt.getopt(inargs, "",["missing=","skip="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for o, a in opts:
        if   o in ("--missing"):
            options['missing'] = open(a,'w')
        elif o in ("--skip"):
            options['skip-file'] = a
            options['skip'] = set()
            f = open(a,'r')
            for line in f:
                sline = line.strip()
                acc = version(sline)
                options['skip'].add(acc)
            f.close()
        else:
            assert False, "unhandled option %s" % o

    return options, args


regex_identifier_gbfile = re.compile('GI:(\d+)')
regex_identifier = re.compile('gi[|](\d+)[|]')
def identifier(s):
    m = regex_identifier_gbfile.search(s)
    if m is None:
        m2 = regex_identifier.search(s)
        if m2 is None:
            raise Exception("'%s' doesn't contain genbank identifier" % (s))
        return m2.group(1)
    return m.group(1)
        

# http://www.ncbi.nlm.nih.gov/Sitemap/sequenceIDs.html
# some XM_ entries have as many as 9 digits
regex_version = re.compile('[A-Z]{1,3}_?\d{5,9}([.]\d+)?')
regex_pdb     = re.compile('pdb[|]([1-9][A-Z0-9]{3})[|]')
regex_sp1     = re.compile('sp[|][OPQ][0-9][A-Z0-9]{3}[0-9]([.]\d+)?[|]')
regex_sp2     = re.compile('sp[|][A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}([.]\d+)?[|]')
regex_prf     = re.compile('prf[|]{2}\w+(\b|$)')
regex_pir     = re.compile('pir[|]{2}\w+(\b|$)')
def version(s):
    m = regex_version.search(s)
    if m is None:
        m = regex_pdb.search(s)
        if m is None:
            m = regex_sp1.search(s)
            if m is None:
                m = regex_sp2.search(s)
                if m is None:
                    m = regex_prf.search(s)
                    if m is None:
                        m = regex_pir.search(s)
                        if m is None:
                            raise Exception("'%s' doesn't contain genbank version, pdb idcode, swiss-prot version, protein research foundation version, or PIR version" % (s))
                        return m.group(0)
                    return m.group(0)
                return m.group(0)
            return m.group(0)
        return m.group(0)
    return m.group(0)



def tax_string_write_lineage(lineage):
    return "\n".join([ "\t".join(map(str,tup)) for tup in lineage ])



def lineage_id(id,taxdb):
    l = []
    nid = id
    while nid != 1:
        l.append(nid)
        nid = taxdb[nid][0]
    return l[::-1]


def write_lineage(taxid,fn,taxdb):
    f_tax = open(fn,'w')
    lineage_ids = lineage_id(taxid,taxdb)
    lineage = [ (tid,taxdb[tid][3],taxdb[tid][1]) for tid in lineage_ids ] # taxid, name, rank
    tax_write_lineage(f_tax,lineage)
    f_tax.close()


def create_lineage(taxid,fn,taxdb):
    lineage_ids = lineage_id(taxid,taxdb)
    lineage = [ (tid,taxdb[tid][3],taxdb[tid][1]) for tid in lineage_ids ] # taxid, name, rank
    return lineage


def parse6(line):
    """-outfmt 6 (tab-delimited"""
    sline    = line.rstrip("\n")
    fields   = sline.split("\t")
    hit      = {}
    hit['qseqid']   = fields[0]
    hit['sseqid']   = fields[1]
    hit['pident']   = float(fields[2])
    hit['length']   = int(fields[3])
    hit['mismatch'] = int(fields[4])
    hit['gapopen']  = int(fields[5])
    hit['qstart']   = int(fields[6])
    hit['qend']     = int(fields[7])
    hit['sstart']   = int(fields[8])
    hit['send']     = int(fields[9])
    hit['evalue']   = float(fields[10])
    hit['bitscore'] = float(fields[11])
    return hit


def main():
    # PARSE OPTIONS
    options = {}
    options['missing'] = None
    options['skip'] = None
    options['skip-file'] = None
    options, args = parse_getopts(sys.argv[1:],options)

    time1 = time.time()
    # READ TAXDB
    taxnames = {}
    fn_tax = 'data/taxonomy/scientific_names.dmp'
    if os.path.exists(fn_tax):
        f_tax = open(fn_tax,'r')
        for line in f_tax:
            sline = line.strip("\t|\n")
            fields = sline.split("\t|\t")
            taxid = int(fields[0])
            name = fields[1]
            uname = fields[2]
            taxnames[taxid] = (name,uname)
        f_tax.close()
    time2 = time.time()
    sys.stderr.write( "%ds %s\n" % (time2-time1,fn_tax) )
    sys.stderr.flush()

    time1 = time.time()
    # READ TAXDB
    taxdb = {}
    fn_tax = 'data/taxonomy/nodes.dmp'
    if os.path.exists(fn_tax):
        f_tax = open(fn_tax,'r')
        for line in f_tax:
            sline = line.strip("\t|\n")
            fields = sline.split("\t|\t")
            taxid = int(fields[0])
            ptaxid = int(fields[1])
            rank = fields[2]
            division = fields[3]
            name = taxnames[taxid][0]
            uname = taxnames[taxid][1]
            taxdb[taxid] = (ptaxid,rank,division,name,uname)
        f_tax.close()
    time2 = time.time()
    sys.stderr.write( "%ds %s\n" % (time2-time1,fn_tax) )
    sys.stderr.flush()

    time1 = time.time()
    # READ TAXDB
    gi_taxid_db = {}
    fn_tax = 'data/taxonomy/gi_taxid_contigs.dmp'
    if os.path.exists(fn_tax):
        f_tax = open(fn_tax,'r')
        for line in f_tax:
            sline = line.strip("\n")
            fields = sline.split("\t")
            gi = int(fields[0])
            taxid = int(fields[1])
            gi_taxid_db[gi] = taxid
        f_tax.close()
    time2 = time.time()
    sys.stderr.write( "%ds %s\n" % (time2-time1,fn_tax) )
    sys.stderr.flush()

    time1 = time.time()
    # READ CDHIT
    clusters = {}
    re_contig = re.compile(' >(s\d{4}_c\d{6})[.]{3}')
    fn_cluster = args[0]
    f_cluster = open(fn_cluster) # *.bak.clstr
    for line in f_cluster:
        # 12162^I3574nt, >s0001_c000001... at 1:3571:407:3968/+/99.58%$
        # 11386^I3264nt, >s0001_c000002... at 1:3264:1:3264/+/100.00%$
        # 3294^I2985nt, >s0001_c000003... at 1:2985:7041:10025/+/100.00%$
        # 8457^I2897nt, >s0001_c000004... *
        # 6307^I2803nt, >s0001_c000005... at 2803:1:718:3520/-/99.64%$
        sline = line.strip("\n")
        fields = sline.split("\t")
        i_cluster = int(fields[0])
        if i_cluster not in clusters:
            clusters[i_cluster] = set()
        m = re_contig.search(fields[1])
        contig = m.group(1)
        clusters[i_cluster].add(contig)
    f_cluster.close()
    time2 = time.time()
    sys.stderr.write( "%ds %s\n" % (time2-time1,fn_cluster) )
    sys.stderr.flush()

    time1 = time.time()
    # READ BLASTN
    hits = defaultdict(Counter)
    fn_hits = args[1]
    f_hits = open(fn_hits)
    for line in f_hits:
        hit = parse6(line)
        qseqid = hit['qseqid']
        sseqid = hit['sseqid']
        hits[qseqid][sseqid] += 1
    f_hits.close()
    time2 = time.time()
    sys.stderr.write( "%ds %s\n" % (time2-time1,fn_hits) )
    sys.stderr.flush()

    time1 = time.time()
    # CALCULATE
    # following output is created for each cluster:
    # defline: the counts for each defline
    # tax-level: the counts for each level of taxonomy
    # tax-lineage: the counts for each complete lineage
    # aggregate: the number of hits, accessions, tax-levels and tax-lineages
    lineages = {}
    newlinetab = maketrans('\n','\t')
    ilen = len(clusters.keys())
    for i in sorted(clusters.keys()):
        if i % 1000 == 0:
            time2 = time.time()
            timediff = time2 - time1
            sys.stderr.write( "Cluster: %d / %d; Time: %ds\n" % (i,ilen,timediff) )
            sys.stderr.flush()
            time1 = time.time()

        cnt_hits = Counter() # counts the number of hits to all accessions of a given contig
        for c in clusters[i]:
            cnt_hits.update(hits[c])

        cnt_complete_lineages = Counter()
        cnt_level_lineages = Counter()
        for hit,cnt in cnt_hits.most_common():
            acc = version(hit)
            gi = int(identifier(hit))

            if acc in options['skip']:
                # don't try to process this it doesn't work and is a waste of time
                continue

            if acc not in lineages:
                lineage = create_lineage(taxid,fn_tax,taxdb)
                lineages[acc] = lineage
            else:
                lineage = lineages[acc]
            lineage_str = tax_string_write_lineage(lineage)

            cnt_level_l = Counter(lineage)
            for l,cnt_l in cnt_level_l.items():
                cnt_level_lineages[l] += cnt
            cnt_complete_lineages[lineage_str] += cnt


        for tup,cnt in cnt_level_lineages.most_common():
            taxid = tup[0]
            name = tup[1]
            rank = tup[2]
            sys.stdout.write( "%d\ttax-level\t%d\t%s\t%s\t%d\n" % (i,taxid,name,rank,cnt) )


if __name__ == "__main__":
    main()

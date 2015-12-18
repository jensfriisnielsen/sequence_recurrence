PAR = M0T28B0p1g0r1c090aS090G1
PARS = "-c 0.90 -aS 0.90 -G 1"
DATA = data
FEAT = $(DATA)/features
BONFERRONI = 3.01e-10

clustering: .FORCE
	scripts/cdhit.sh -M 0 -T 28 -B 0 -p 1 -g 0 -r 1 $(PARS)

feature_associations: .FORCE
	scripts/feature_association.sh results/clustering/contigs_all.cdhit.$(PAR).clstr

taxonomy: .FORCE
	scripts/cluster_taxonomy.sh results/clustering/contigs_all.cdhit.$(PAR).bak.clstr

bonferoni_significant.txt: results/feature_associations/pval.txt
	gawk -F $$'\t' '{if ($$7 < $(BONFERRONI)) print $$0}' results/feature_associations/pval.txt | grep -f $(FEAT)/featuresall > $@

bonferoni_significant.annotated.txt: bonferoni_significant.txt 
	scripts/annotate.sh bonferoni_significant.txt results/feature_associations/pval.txt > $@

bonferoni_significant.annotated.strongest.txt: bonferoni_significant.annotated.txt
	sort -t$$'\t' -k7,7g bonferoni_significant.annotated.txt | scripts/keep_first.py -f 1 > $@

bonferoni_significant.annotated.strongest.clstr.lst: bonferoni_significant.annotated.strongest.txt
	cut -f1 bonferoni_significant.annotated.strongest.txt | sort -u -k2,2n > $@

bonferoni_significant.annotated.strongest.clstr.bio.lst: bonferoni_significant.annotated.strongest.txt
	grep -f $(FEAT)/featuresbio bonferoni_significant.annotated.strongest.txt | cut -f1 | sort -u -k2,2n > $@

bonferoni_significant.annotated.strongest.clstr.met.lst: bonferoni_significant.annotated.strongest.txt
	grep -f $(FEAT)/featuresmet bonferoni_significant.annotated.strongest.txt | cut -f1 | sort -u -k2,2n > $@

bonferoni_significant.annotated.strongest.clstr.tec.lst: bonferoni_significant.annotated.strongest.txt
	grep -f $(FEAT)/featurestec bonferoni_significant.annotated.strongest.txt | cut -f1 | sort -u -k2,2n > $@

table.txt: bonferoni_significant.annotated.strongest.clstr.txt
	ln -s bonferoni_significant.annotated.strongest.clstr.txt table.txt 

.FORCE:

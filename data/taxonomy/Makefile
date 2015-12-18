all: download unzip cat divide

divide: unzip
	 ./divide_taxid.py division.dmp nodes.dmp

cat: unzip
	cat gi_taxid_{nucl,prot}.dmp > gi_taxid_contigs.dmp

unzip: download
	gunzip gi_taxid_prot.dmp.gz gi_taxid_nucl.dmp.gz scientific_names.dmp.gz nodes.dmp.gz

download:
	wget -m ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
	mv ftp.ncbi.nlm.nih.gov/pub/taxonomy/* .


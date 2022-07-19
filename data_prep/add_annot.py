#!/usr/bin/env python3
"""
Parse GTF annotation file and add to
WebSTR database

Usage:
./add_annot.py <build> <gtf> <dbfile>
"""

import gzip
import sqlite3
import sys

try:
	BUILD = sys.argv[1]
	GTF = sys.argv[2]
	DBFILE = sys.argv[3]
except:
	sys.stderr.write(__doc__)
	sys.exit(1)

def GetAnnotID():
	"""
	Get max ID of last annotation item
	so we can start counting from there
	"""
	conn = sqlite3.connect(DBFILE)
	curr = conn.cursor()
	nrow = curr.execute("select count(*) from GENEANNOTATIONS").fetchall()
	if nrow[0][0] == 0: return 0
	res = curr.execute("select max(feature_id) from GENEANNOTATIONS").fetchall()
	return int(res[0][0])

def ParseGTFLine(gtfline, build, annot_id):
	"""
	Parameters
	----------
	gtfline (str): Line from the GTF file

	Returns
	-------
	db_items (list): list of items to
	   insert into GENEANNOTATIONS

	db_items should have:
		feature_id int,
		genome_build varchar(10),
		feature_chrom varchar(10),
		feature_type varchar(20),
		feature_start int,
		feature_end int,
		feature_strand int,
		gene_name varchar(20),
		gene_id varchar(20),
		gene_type varchar(20)
	"""
	items = gtfline.strip().split("\t")
	# Parse feature coordinates
	feature_chrom = items[0]
	if build in ["hg19","hg38","rn7","mm10"] and "chr" not in feature_chrom:
		feature_chrom = "chr"+feature_chrom
	feature_type = items[2]
	feature_start = int(items[3])
	feature_end = int(items[4])
	feature_strand = int(items[6]=="+")

	# Parse additional metadata
	gene_name = None
	gene_id = None
	gene_type = None
	metadata = items[8].strip().split(";")
	for item in metadata:
		if item.strip() == "": continue
		key, val = item.strip().split()
		if key == "gene_name":
			gene_name = val.strip('"').strip()
		if key == "gene_id":
			gene_id = val.strip('"').strip()
		if key == "gene_biotype":
			gene_type = val.strip('"').strip()

	db_items = [annot_id, build, feature_chrom, feature_type, \
		feature_start, feature_end, feature_strand, \
		gene_name, gene_id, gene_type]

	return db_items

annot_id = GetAnnotID() + 1

with gzip.open(GTF, "r") as f:
	for line in f:
		line = line.decode("utf-8") 
		if line.startswith('#'): continue
		db_items = ParseGTFLine(line, BUILD, annot_id)
		#print(db_items)
		conn = sqlite3.connect(DBFILE)
		curr = conn.cursor()
		curr.execute('insert into GENEANNOTATIONS values (?,?,?,?,?,?,?,?,?,?)', db_items)		
		conn.commit()
		annot_id += 1

sys.exit(0)

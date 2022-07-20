"""
Utility functions for the WebSTR region page
"""

import utils

MAXREGIONSIZE = 1500000

def GetGenesFromRegion(db, genome, region):
	"""
	Get list of genes in a region

	Parameters
	----------
	db : Session
	   sqlalchemy session
	genome : str
	   Genome build id
	region : tuple of (str, int, int)
	    Contains chrom, start, end of region

	Returns
	-------
	genes : list of str
	    List of gene names to display
	"""
	# TODO not implemented
	return []

def GetRegionFromGenes(db, genome, genes):
	"""
	Get genomic region spanned by a set of genes

	Parameters
	----------
	db : Session
	   sqlalchemy session
	genome : str
	   Genome build id
	genes : list of str
	    List of gene names to display

	Returns
	-------
	region : tuple of (str, int, int)
	    Contains chrom, start, end of region
	"""
	chrom = "NA"
	start = 0
	end = 0
	# TODO not implemented
	return (chrom, start, end)

def ParseQuery(db, genome, query):
	"""
	Parse query for region view

	Parameters
	----------
	db : Session
	   sqlalchemy session
	genome : str
	   Genome build id
	query : str
	   Query string. Should be either:
	   (1) Gene name
	   (2) Region (chr:start-end)
	   (3) Region (chrom start end)

	Returns
	-------
	genes : list of str
	    List of gene names to display
	region : tuple of (str, int, int)
	    Contains chrom, start, end of region

	If invalid query, return [], None
	"""
	genes = []
	region = None
	if len(query.strip().split()) == 3:
		chrom, start, end = query.strip().split()
		try:
			start = int(start)
		except ValueError: return [], None
		try:
			end = int(end)
		except ValueError: return [], None
		region = (chrom, start, end)
		genes = GetGenesFromRegion(db, genome, region)
	elif ":" in query:
		try:
			chrom = "chr"+query.split(":")[0].replace("chr","").replace("CHR","")
			start = int(query.split(":")[1].split("-")[0])
			end = int(query.split(":")[1].split("-")[1])
		except: return [], None
		region = (chrom, start, end)
		genes = GetGenesFromRegion(db, genome, region)
	elif len(query.strip().split()) > 1:
		return [], None
	else:
		if query.startswith("ENS"):
			gene_attrib = "gene_id"
			gene = utils.GetGeneFromENS(db, genome, query)
			if gene is not None:
				genes = [gene]
			else:
				genes = []
			region = GetRegionFromGenes(db, genome, genes)
		else:
			gene_attrib = "gene_name"
			genes = [query]
			region = GetRegionFromGenes(db, genome, genes)
	return genes, region

def GetRegionItems(db, genome, query):
	"""
	Get items needed to render the region-level page

	Parameters
	----------
	db : Session
	   sqlalchemy session
	genome : str
	   Genome build id
	query : str
	   Query string. Should be either:
	   (1) Gene name
	   (2) Region (chr:start-end)
	   (3) Region (chrom start end)

	Returns
	-------
	template_items : dict
	   Dictionary passed to HTML rendering
	passed : bool
	   Return False if we could not parse the query
	   or failed to get region display info
	err : str
	   Contains an error message to display if passed=False
	"""
	# Initialize return values
	template_items = {
		"genome": genome,
		"query": query
	}
	passed = True
	msg = ""

	# Parse query
	genes, region = ParseQuery(db, genome, query)
	if region is None:
		return {}, False, "Failed to parse query"
	if (region[2]-region[1]) > MAXREGIONSIZE:
		return {}, False, "Region size bigger than " + str(MAXREGIONSIZE)
	template_items["genes"] = ",".join(genes)
	template_items["region"] = ",".join([str(item) for item in region])
	print(template_items) # TODO remove, for debugging
	return template_items, passed, msg
"""
Functions for interacting with UCSC API
for retrieving gene annotations

See: https://genome.ucsc.edu/goldenPath/help/api.html
"""

import requests

# Possible starting queries: 
# https://api.genome.ucsc.edu/search?search=brca1&genome=hg38
# https://api.genome.ucsc.edu/search?search=ENSG00000012048&genome=hg38
# https://api.genome.ucsc.edu/search?search=BRCA1&genome=rn7
def GetRegionFromGeneQuery(gene_query, genome_build):
	"""
	Given a gene query (e.g. BRCA1, or ENSG00000012048)
	and a genome build, return the chr:start-end of the region
	"""
	request_url = "" # TODO
	resp = requests.get(request_url)
	return None # TODO parse query results

# Possible starting queries:
# https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeq;chrom=chr1;start=11868;end=31109
def GetGeneFeaturesFromRegionQuery(region_query, genome_build):
	"""
	Given a region query (chr:start-end), get a list of genes, where
	each gene we have a list of exons with start, end
	e.g.
	genes = [
		{"start": 12345, "end": 4567, "exons": [
			{"start": 12345, "end": 12346}, {"start": 12350, "end": 12355"}...
		]"}
	]
	This should return something similar to what ExtractGeneFeaturesAPI
	currently returns
	"""
	return None # TODO
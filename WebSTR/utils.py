"""
Helper functions for WebSTR.py
"""

from sqlalchemy.orm import Session
import models, schemas

def get_genomes(db: Session, skip: int = 0, limit: int = 100):
    genomes = db.query(models.Genome).offset(skip).limit(limit).all()
    return [genome.genome_build for genome in genomes]

def GetGeneFromENS(db: Session, genome : str, ensid : str):
	"""
	Get the gene_name from a gene_id

	Parameters
	----------
	genome : str
	   Genome build id
	ensid : str
	   Ensemble gene_id

	Returns
	-------
	gene_name : str
	   Gene name
	   Return None if gene_id not found	
	"""
	annot = db.query(models.GeneAnnotations).filter(models.GeneAnnotations.gene_id == ensid and \
		models.GeneAnnotations.genome_build == genome).first()
	if annot is None: return None
	return annot.gene_name

from typing import List, Union

from pydantic import BaseModel

class GenomeBase(BaseModel):
	build: str

class Genome(GenomeBase):
    genome_id: int
    genome_build: str

    class Config:
    	orm_mode = True

class GeneAnnotationsBase(BaseModel):
	build: str

class GeneAnnotations(GeneAnnotationsBase):
    feature_id: int
    genome_build: str
    feature_chrom: str
    feature_type: str
    feature_start: int
    feature_end: int
    feature_strand: int
    gene_name: str
    gene_id: str
    gene_type: str

    class Config:
    	orm_mode = True
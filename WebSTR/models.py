from sqlalchemy import Boolean, Column, ForeignKey, Integer, String
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class Genome(Base):
	__tablename__ = "GENOMES"
	genome_id = Column(Integer, primary_key=True, index=True)
	genome_build = Column(String)

class GeneAnnotations(Base):
	__tablename__ = "GENEANNOTATIONS"
	feature_id = Column(Integer, primary_key=True, index=True)
	genome_buid = Column(String, ForeignKey("GENOMES.genome_build")),
	feature_chrom = Column(String)
	feature_type = Column(String)
	feature_start = Column(Integer)
	feature_end = Column(Integer)
	feature_strand = Column(Integer)
	gene_name = Column(String)
	gene_id = Column(String)
	gene_type = Column(String)
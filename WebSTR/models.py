from sqlalchemy import Boolean, Column, ForeignKey, Integer, String
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class Genome(Base):
	__tablename__ = "GENOMES"
	genome_id = Column(Integer, primary_key=True, index=True)
	genome_build = Column(String)
"""
Helper functions for WebSTR.py
"""

from sqlalchemy.orm import Session
import models, schemas

def get_genomes(db: Session, skip: int = 0, limit: int = 100):
    genomes = db.query(models.Genome).offset(skip).limit(limit).all()
    return [genome.genome_build for genome in genomes]
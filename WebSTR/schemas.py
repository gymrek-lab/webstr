from typing import List, Union

from pydantic import BaseModel

class GenomeBase(BaseModel):
	build: str

class Genome(GenomeBase):
    id: int
    build: str

    class Config:
    	orm_mode = True
#!/usr/bin/env python3
"""
WebSTR v2 database application
"""

import argparse
from fastapi import Depends, FastAPI, Request, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
import os
from sqlalchemy.orm import Session
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from starlette.responses import FileResponse 
import uvicorn

import utils
import region

#################### Paths ########################################
BASE_DIR = "/Users/melissagymrek/workspace/webstr/"
STATIC_DIR = os.path.join(BASE_DIR, "WebSTR", "static")
TEMPLATE_DIR = os.path.join(BASE_DIR, "WebSTR", "templates")

# NOTE: change this to modify database location
SQLALCHEMY_DATABASE_URL = "sqlite:////Users/melissagymrek/workspace/webstr/data_prep/webstr2.db"

#################### Set up the database ##########################
import utils, models, schemas

# NOTE: if not using sqlite, set "check_same_thread": True
engine = create_engine(
    SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False}
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

models.Base.metadata.create_all(bind=engine)

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

#################### Set up the app ##########################

app = FastAPI()
app.mount("/static", StaticFiles(directory=STATIC_DIR))
templates = Jinja2Templates(directory=TEMPLATE_DIR)

#################### Functions to render pages ##########################

@app.get("/", response_class=HTMLResponse)
def WebSTRHome(request: Request, db: Session = Depends(get_db)):
	template_items = {
		"request": request,
		"genomes": utils.get_genomes(db)
	}
	return templates.TemplateResponse("homepage.html", template_items)

@app.get("/region/", response_class=HTMLResponse)
@app.get("/region/{genome}/{query}", response_class=HTMLResponse)
async def WebSTRRegion(request: Request, genome: str, query: str, \
	db: Session = Depends(get_db)):
	if query.strip() == "":
		raise HTTPException(status_code=404, detail="No query found")
	if genome.strip() == "":
		raise HTTPException(status_code=404, detail="No genome selected")
	if genome.strip() not in utils.get_genomes(db):
		raise HTTPException(status_code=404, detail="Invalid genome")
	template_items, passed, err = region.GetRegionItems(db, genome, query)
	if not passed:
		raise HTTPException(status_code=404, detail=err)
	template_items["request"] = request
	return templates.TemplateResponse("region.html", template_items)

@app.get("/locus/", response_class=HTMLResponse)
@app.get("/locus/{genome}/{trsetid}/{strid}", response_class=HTMLResponse)
async def WebSTRLocus(request: Request, genome: str, trsetid: str, strid: str):
	template_items = {
		"request": request,
		"genome": genome,
		"trsetid": trsetid,
		"strid": strid
	}
	return templates.TemplateResponse("locus.html", template_items)

#################### Set up and run the server ###############
if __name__ == "__main__":
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--host", help="Host to run app", type=str, default="0.0.0.0")
	parser.add_argument("--port", help="Port to run app", type=int, default=5000)
	args = parser.parse_args()
	uvicorn.run(app, host=args.host, port=args.port)

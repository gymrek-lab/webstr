#!/usr/bin/env python3
"""
WebSTR v2 database application
"""

import argparse
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from starlette.responses import FileResponse 
import uvicorn

#################### Set up the app ##########################
app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")

#################### Functions to render pages ##########################

@app.get("/", response_class=HTMLResponse)
async def WebSTRHome(request: Request):
    return templates.TemplateResponse("homepage.html", {"request": request})

#################### Set up and run the server ###############
if __name__ == "__main__":
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--host", help="Host to run app", type=str, default="0.0.0.0")
	parser.add_argument("--port", help="Port to run app", type=int, default=5000)
	args = parser.parse_args()
	uvicorn.run(app, host=args.host, port=args.port)
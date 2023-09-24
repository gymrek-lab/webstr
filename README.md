# WebSTR web browser (front-end)

This repository contains code and instructions for [WebSTR](http://webstr.ucsd.edu/) - web browser of Human genome-wide variation in Short Tandem Repeats (STRs). Our goal is to make large STR genotype datasets used by the broader genomics community by facilitating open access to this data.

WebSTR is the result of collaboration between two scientific groups [Maria Anisimova’s Lab](https://github.com/acg-team) and [Melissa Gymrek’s Lab](https://github.com/gymrek-lab).

Source code for the WebSTR-API can be found here: https://github.com/acg-team/webSTR-API

## Building docker container for debug and production use

* You will need docker installed on your system.
* Docker container has several stages that allow to develop in container, run WebStr.py locally (debug stage) and production with gunicorn.
* Then run docker build n the project folder
```
docker build --target debug -t webstr:debug . 
docker build -t webstr . 

# example of running containers
docker run -it --rm -t webstr:debug
docker run -it --rm -t webstr
```

## Instructions on how to develop WebSTR web app locally with docker (for development)

#### 1. For local development you can use [DevContainer](https://code.visualstudio.com/docs/devcontainers/containers) extention for Visual Studio Code, .devcontainer file is already in repository.
#### 2. Or you can rename .vscode/lanuch.json.example to .vscode/lanuch.json and .vscode/tasks.json.example to .vscode/tasks.json and use Visual Studio Code Debug in docker feature.

## Instructions on how to set-up WebSTR web app locally without docker (for development)

### Set up python3 and virtualenv on your machine:
[For Mac, follow instructions here.](https://gist.github.com/pandafulmanda/730a9355e088a9970b18275cb9eadef3)

### Create new virtual env with python3 and install all the requirements with the following command:
`pip install -r requirements.txt`

### Copy data files to data directory

We provide neccesssary data files upon request and currently working on incorporating them fully into the backend. 

You will need:
* dbSTR.db  (contains older panels mapped to hg19 assembly)
* hg19.fa and hg19.fai
* hg38.fa and hg38.fai

###  Start web server

To run for testing and development:
```
export DATAPATH=*full data directory path*
python ./WebSTR/WebSTR.py --host 0.0.0.0 --port <port>
```

When working with docker container you will need to mount files to the container and set the DATAPATH variable to the mounted directory path inside the container.
```
docker build -t  webstr_frontend .
docker run --mount type=bind,src=PATHTOYOURDATA,dst=/data --env DATABATH=/data webstr_frontend
```

For production mode, we use gunicorn + nginx. 

### WebSTR Backend - database and API

By default WebSTR will be using WebSTR-API hosted on our server. If you would like to set up the database and WebSTR backend as well (for example if you would like to make a new endpoint or add your own data to the database), please follow instructions to set up WebSTR database and API locally provided here [https://github.com/acg-team/webSTR-API](https://github.com/acg-team/webSTR-API). Proceed then to modify API_URL variable in WebSTR browser so the app communicates with your local instance of WebSTR-API. 

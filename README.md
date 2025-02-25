# WebSTR web browser (front-end)

This repository contains code and instructions for [WebSTR](http://webstr.ucsd.edu/) - web browser of Human genome-wide variation in Short Tandem Repeats (STRs). Our goal is to make large STR genotype datasets used by the broader genomics community by facilitating open access to this data.

WebSTR is the result of collaboration between two scientific groups [Maria Anisimova’s Lab](https://github.com/acg-team) and [Melissa Gymrek’s Lab](https://github.com/gymrek-lab).

Source code for the WebSTR-API can be found here: https://github.com/acg-team/webSTR-API

## Contributing

If you would like to make changes to WebSTR, you must:
1. Create a new branch off of main and make your edits
2. Submit a pull request to merge your new branch to main
3. The pull request must be reviewed and approved by a WebSTR developer prior to merging with main.

## Instructions for setting up WebSTR for local development (without docker)

1. Set up python3 and virtualenv on your machine:
[For Mac, follow instructions here.](https://gist.github.com/pandafulmanda/730a9355e088a9970b18275cb9eadef3)

2. Create a new conda environment with python3 and install all the requirements with the following command:

`conda env create -f environment.yml`

3. Copy data files to data directory

WebSTR looks for certain files at `BASEPATH`. You will need the following:
* `$BASEPATH/hg19/hg19.fa` (hg19 reference genome)
* `$BASEPATH/hg38/hg38.fa` (hg38 reference genome)
* `$BASEPATH/dbSTR.db` (legacy hg19 version of the database, can be obtained [here](https://drive.google.com/file/d/1Lm-nx-G2V726Re67EnOHTWTYhgDo-W38/view?usp=sharing))

The hg38 database is managed by the backend and API. This is described below, along with instructions on how to test WebSTR with a non-production version of the backend database.

4. To run for testing and development:

```
git clone https://github.com/gymrek-lab/webstr
cd webstr
# optionally, checkout a specific branch to test
export BASEPATH=*full data directory path*
export FLASK_DEBUG=1 # run in debug mode
python ./WebSTR/WebSTR.py --host 0.0.0.0 --port <port>
```

You can then access the application at `localhost:$port` in your web browser.

## Instructions for setting up WebSTR for local development (with docker)

1. Clone the WebSTR repository

```
git clone https://github.com/gymrek-lab/webstr
cd webstr
# optionally, checkout a specific branch to test
```

2. Build the debug version of the docker (requires docker to be installed)

```
docker build --target debug -t webstr:debug .
```

3. Run the docker

You will need to mount files to the container and set the `$BASEPATH` variable:

```
docker run --mount type=bind,src=${BASEPATH},dst=/data --env BASEPATH=/data  -it --rm -t webstr:debug
```

You can then access the application at `localhost:5000` in your web browser.

## Running WebSTR in production mode with docker

For production mode, we use gunicorn + nginx. 

1. Build the docker in production mode

```
docker build -t webstr .
```

2. Run the docker

```
docker run --mount type=bind,src=${BASEPATH},dst=/data --env BASEPATH=/data  -it --rm -t webstr
```

## WebSTR Backend - database and API

WebSTR access its database through an API. By default, it uses https://webstr-api.lsfm.zhaw.ch. To set a custom location for the API, for example if you are testing a non-production version of the database, you can set a different location using the `WEBSTR_API_URL` environment variable.

The code for the WebSTR-API backend is maintained [here](https://github.com/acg-team/webSTR-API). If you have your own version of the database and want to test the backend locally you will need to run something like the following from inside the `webSTR-API` repo:

```
# Set the path to your local database
export DATABASE_URL="postgres://webstr:webstr@localhost:5432/strdb"

# Launch the API
uvicorn strAPI.main:app --host=0.0.0.0 --port=${PORT:-5000} --reload
```

Now the API should be available at localhost:5000. Before running WebSTR you can set:

```
export WEBSTR_API_URL=http://0.0.0.0:5000
```
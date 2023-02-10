# webstr
WebSTR database (hg19) and web application code

## Running the webserver

### Set up python3 and virtualenv on your machine:
[For Mac, follow instructions here.](https://gist.github.com/pandafulmanda/730a9355e088a9970b18275cb9eadef3)

### Create new virtual env with python3 and install all the requirements with the following command:
`pip install -r requirements.txt`

### Copy data files to data directory

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

For production mode, use a WSGI server like waitress or gunicorn. 

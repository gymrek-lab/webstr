## Instructions for testing WebSTR locally

1. Connect to snorlax (or whatever server you are using) allowing port forwarding on 5001:

```
ssh -L 5001:localhost:5001 $USER@snorlax.ucsd.edu
```

2. If you are using the production version of the database, set `WEBSTR_API_URL`:

```
export WEBSTR_API_URL=https://webstr-api.lsfm.zhaw.ch
```

Otherwise, follow step 3.
   
3. If you are testing a local version of the database, you must first set up the API. If you have already cloned the webSTR-API repo, `cd` to that directory and do `git pull` to get the latest changes. If not do:

```
git clone https://github.com/acg-team/webSTR-API
cd webSTR-API
export DATABASE_URL="postgres://webstr:webstr@localhost:5432/strdb"
uvicorn strAPI.main:app --host=0.0.0.0 --port=${PORT:-5000} --reload
export WEBSTR_API_URL=http://0.0.0.0:5000
```

3. If you have already cloned the WebSTR repo, `cd` to that directory and do `git pull` to get the latest changes. If not do:

```
git clone https://github.com/gymrek-lab/webstr/
cd webstr
```

Optionally, checkout the branch you want to test.

4.   Start the webstr site

```
python ./WebSTR/WebSTR.py --host 0.0.0.0 --port 5001
```

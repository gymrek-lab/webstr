## Instructions for Running webSTR on snorlax

**Please note that port --5000 must be used for the API and port --5001 should be used for snorlax**

1. Connect to snorlax allowing port forwarding on 5001:

```
ssh -L 5001:localhost:5001 mgymrek@snorlax.ucsd.edu
```

2. If you have alread cloned the webSTR-API repo, `cd` to that directory and do `git pull` to get the latest changes. If not do:

```
git clone https://github.com/acg-team/webSTR-API
cd webSTR-API
```

3. Start the API

```
export DATABASE_URL="postgres://webstr:webstr@localhost:5432/strdb"
uvicorn strAPI.main:app --host=0.0.0.0 --port=${PORT:-5000} --reload
```

4. If you have already cloned the WebSTR repo, `cd` to that directory and do `git pull` to get the latest changes. If not do:

```
git clone https://github.com/gymrek-lab/webstr/
cd webstr
```

5.   Start the webstr site

```
python ./WebSTR/WebSTR.py --host 0.0.0.0 --port 5001
```

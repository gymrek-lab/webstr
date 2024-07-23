## Instructions for Running webSTR on snorlax

**Please note that port --5000 must be used for the API and port --5001 should be used for snorlax**

- Install any needed dependencies 

- Start API
  ```
  export DATABASE_URL="postgres://webstr:webstr@localhost:5432/strdb"
  ```
  ```
  uvicorn strAPI.main:app --host=0.0.0.0 --port=${PORT:-5000} --reload
  ```

- Run on snorlax
  ```
  python ./WebSTR/WebSTR.py --host 0.0.0.0 --port 5001
  ```

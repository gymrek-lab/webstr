# ---------------------------------------------
# Image with packages installed
# ---------------------------------------------
FROM python:3.11-slim as base

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

COPY requirements.txt .

RUN apt-get update && \
    apt-get install -y gcc python3-dev build-essential  && \
    pip install -r requirements.txt --no-cache-dir  && \
    apt-get remove -y gcc python3-dev build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# ---------------------------------------------
# Image with source code that runs python entrypoint
# ---------------------------------------------
FROM base as local
COPY . /app
WORKDIR /app/WebSTR

ARG PORT=5000
ENV FLASK_PORT=${PORT}
# Overrides ling to WebSTR api that python modules use
ENV WEBSTR_API_URL=https://webstr-api.lsfm.zhaw.ch
# Tells flask command where to look for the app
ENV FLASK_APP=WebSTR:server
# Turns on Flask debug
ENV FLASK_DEBUG=1

EXPOSE ${FLASK_PORT}
ENTRYPOINT ["python", "WebSTR.py"]


# -------------------------------------------------------------------
# Image with source code that runs in debug without gunicorn
# -------------------------------------------------------------------
FROM local as debug
# Turns on Flask debug
ENV FLASK_DEBUG=1


# -------------------------------------------------------------------
# Image with source code that runs prod setup with gunicorn
# -------------------------------------------------------------------
FROM local as prod
# Turns off Flask debug
ENV FLASK_DEBUG=0
# Creates a non-root user with an explicit UID and adds permission to access the /app folder
# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
RUN adduser -u 1000 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser
# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
ENTRYPOINT ["bash", "-c", "gunicorn --bind 0.0.0.0:${FLASK_PORT} ${FLASK_APP}"]

```
# Set up dependencies
conda create -n webstr python=3.11
conda activate webstr
conda install flask numpy pandas dash plotly dash gunicorn
pip3 install install dash_bio pyfaidx

# Run webstr
cd ~/webstr
git pull

# Test
python3 WebSTR/WebSTR.py --port=5000 

# Run
nohup gunicorn --workers 3 --bind unix:WebSTR.sock -m 007 WebSTR:server &
```
import os
from dbutils import *
from utils import *
import requests
import json
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
import sys

API_URL = os.environ.get("WEBSTR_API_URL",'http://webstr-api.ucsd.edu')

seqbuf = 120
seqbreakline = 100

################ General functions ################
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def GetSTRInfo(str_query, genome_query, dbSTRPath, reffa):
    logging.debug(f"Getting STR info for str_query: {str_query}, genome_query: {genome_query}")

    if genome_query is None or genome_query == "hg19":
        strinfo = GetSTRMetadataHg19(str_query, dbSTRPath)
    else:
        strinfo = GetSTRMetadataAPI(str_query)

    if strinfo is None:
        logging.warning(f"No STR info found for {str_query}")
        return None

    try:
        chrom = strinfo["chrom"]
        start = strinfo["start"]
        end = strinfo["end"]


        lflank = str(reffa[chrom][start-seqbuf:start]).upper()
        strseq = str(reffa[chrom][start:end]).upper()
        rflank = str(reffa[chrom][end:end+seqbuf]).upper()

        strinfo["chrom"] = strinfo["chrom"].replace("chr", "")
        strinfo["seq"] = GetSTRSeqHTML(lflank, strseq, rflank)
        strinfo["motif"] = GetMotifAndComplement(strinfo["motif"])


    except AttributeError as e:
        logging.error(f"Error while processing STR info for {str_query}: {e}")
        return None

    except Exception as e:
        logging.error(f"Error while processing STR info for {str_query}: {e}")
        return None

    return strinfo




def GetSTRSeqHTML(lflank, strseq, rflank, charbreak=50):
    ret = '<font size="3" color="black">...'
    numchar = 0
    for i in range(len(lflank)):
        ret += lflank[i]
        numchar += 1
        if numchar%seqbreakline == 0: ret += "<br>"
    ret += "</font>"
    ret += '<font size="4" color="red"><b>'
    for i in range(len(strseq)):
        ret += strseq[i]
        numchar += 1
        if numchar%seqbreakline == 0: ret += "<br>"
    ret += "</font></b>"
    ret += '<font size="3" color="black">'
    for i in range(len(rflank)):
        ret += rflank[i]
        numchar += 1
        if numchar%seqbreakline == 0: ret += "<br>"
    ret += "</font>..."
    return ret

def GetFreqPlotlyJSON(genome_query, freq_dist):
    if genome_query is None or genome_query == "hg19":
        return GetFreqPlotHg19(freq_dist)
    else:
        return GetFreqPlotHg38(freq_dist)

################ Fetching STR metadata  ################
def GetSTRMetadataHg19(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    squery = ("select str.chrom, str.start, str.end,str.motif, (str.end-str.start+1)/str.period copies from strlocmotif str where str.strid = '{}'").format(strid)
    df = ct.execute(squery).fetchall()
    if len(df) == 0: return None
    return {"chrom": df[0][0],
            "start": df[0][1],
            "end": df[0][2],
            "motif": df[0][3],
            "copies": df[0][4],
            "gene_name": "",
            "gene_desc": "",
            "crc_data": None,
            "gtex_data": GetGTExInfo(strid, DbSTRPath),
            "mut_data": GetMutInfo(strid, DbSTRPath),
            "imp_data": GetImputationInfo(strid, DbSTRPath),
            "imp_allele_data": GetImputationAlleleInfo(strid, DbSTRPath),
            "seq_data": None,
            "freq_dist": GetFreqSTRInfo(strid, DbSTRPath)
    }

def GetSTRMetadataAPI(repeat_id):
    repeat_url = API_URL + '/repeatinfo/?repeat_id=' + repeat_id
    print(f"Fetching repeat info from: {repeat_url}")  # Debug statement
    resp = requests.get(repeat_url)
    
    # Check if the API call fails or the response is empty
    try:
        repeat = json.loads(resp.text)
    except json.JSONDecodeError:
        sys.stderr.write(f"ERROR: Could not decode JSON for repeat_id {repeat_id}\n")
        return None

    # Check if `repeat` is a dictionary before calling .get()
    if isinstance(repeat, dict):
        strinfo = {
            "chrom": repeat.get("chr", "Unknown"),
            "start": repeat.get("start", "Unknown"),
            "end": repeat.get("end", "Unknown"),
            "motif": repeat.get("motif", "Unknown"),
            "copies": repeat.get("copies", "Unknown"),
            "gene_name": repeat.get("gene_name", "Unknown"),
            "gene_desc": repeat.get("gene_desc", "No description available"),
            "crc_data": None,
            "gtex_data": None,
            "mut_data": None,
            "imp_data": None,
            "imp_allele_data": None,
            "seq_data": GetSeqDataAPI(repeat_id) or "No sequence information currently available",
            "freq_dist": GetFreqSTRInfoAPI(repeat_id) or []
        }
    else:
        print(f"Expected dictionary, got {type(repeat)} instead.")
        return None

    return strinfo


################# Fetch data from hg19 specific tables #################
def GetGTExInfo(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select estr.gene,estr.genename,ti.tissue,estr.beta,estr.beta_se,estr.pval,estr.caviar "
              "from estr_gtex2 estr, "
              "tissues ti, "
              "strlocmotif str "
              "where estr.chrom = str.chrom "
              "and estr.strstart = str.start "
              "and estr.strend = str.end "
              "and estr.signif = 'True' "
              "and estr.tissue_cd = ti.tissue_cd "
              "and strid = '{}' order by estr.caviar desc").format(strid)
    df = ct.execute(gquery).fetchall()
    if len(df) == 0: return None
    return df
 
def GetMutInfo(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select mut.est_logmu_ml, mut.est_beta_ml, mut.est_pgeom_ml, mut.up, mut.down, mut.p, mut.zscore_1, mut.zscore_2"
              " from mutrates mut where mut.str_id = '{}'").format(strid)
    df = ct.execute(gquery).fetchall()
    if len(df) != 1:
        return None
    else:
        df = list(df[0])
        df[0] = 10**df[0]
        return df

def GetImputationInfo(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select imp.loo_concordance,imp.loo_r,imp.wgs_eur_concordance,imp.wgs_eur_r,imp.wgs_afr_concordance,imp.wgs_afr_r,"
              " imp.wgs_eas_concordance,imp.wgs_eas_r"
              " from locstat imp where imp.str_id = '{}'").format(strid)
    df = ct.execute(gquery).fetchall()
    if len(df) != 1:
        return None
    else:
        return list(df[0])

def GetImputationAlleleInfo(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select al.allele, al.r2, al.pval"
              " from allelstat al where al.str_id = '{}'").format(strid)
    df = ct.execute(gquery).fetchall()
    return df




############## Functions  for plotting freq distributions #############
def GetFreqSTRInfo(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select cohort_id, (end-start+1+af.length)/period copies,sum(nvals) nvals from"
              " allelefreq af,"
              " strlocmotif strm"
              " where af.str_id = strm.strid"
              " and str_id = '{}' "
              " group by cohort_id, copies").format(strid)
    df = ct.execute(gquery).fetchall()
    if len(df) == 0: return None
    else: return df


def GetFreqSTRInfoAPI(repeat_id):
    freq_url = API_URL + '/allfreqs/?repeat_id=' + repeat_id
    print(f"Fetching frequency data from: {freq_url}")  # Debug statement
    resp = requests.get(freq_url)

    
    # Attempt to load response as JSON
    try:
        freq_data = json.loads(resp.text)
       
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON for frequency data: {e}, raw response: {resp.text}")
        return None

    if not freq_data:
        print(f"No frequency data found for repeat_id: {repeat_id}")  # Debug statement
        return None

    # Process and return grouped frequency data
    df = pd.DataFrame.from_records(freq_data)
    print(f"Frequency DataFrame for repeat_id {repeat_id}: {df.head()}")  # Debug statement
    print(f"Type of DataFrame: {type(df)}")  # Check DataFrame type

    if df.empty:
        print(f"Empty frequency data for repeat_id: {repeat_id}")  # Debug statement
        return None

    df["percentage"] = df["frequency"] * 100
    df["copies"] = df["n_effective"]
    grouped_df = df[["population", "copies", "percentage"]].groupby(by="population")
 
    return grouped_df






def GetFreqPlotHg19(freq_dist):
    data1 = pd.DataFrame(np.array(freq_dist).reshape(-1,3), columns = list("abc"))
    minx = min(data1['b'])-1
    maxx = max(data1['b'])+1
    x1=data1.loc[data1['a'] == 1]
    x2=data1.loc[data1['a'] == 2]
    x3=data1.loc[data1['a'] == 3]
    x4=data1.loc[data1['a'] == 4]
    #for Cohort, X in data1.groupby('a'):
    trace1 = go.Bar(
        x=x1['b'],
        y=x1['c'],
        name = "Gtex"
    )

    trace2 = go.Bar(
        x=x2['b'],
        y=x2['c'],
        xaxis='x2',
        yaxis='y2',
        name = "1000 Genomes Africa"
    )

    trace3 = go.Bar(
        x=x3['b'],
        y=x3['c'],
        xaxis='x3',
        yaxis='y3',
        name = "1000 Genomes East Asia"
    )

    trace4 = go.Bar(
        x=x4['b'],
        y=x4['c'],
        xaxis='x4',
        yaxis='y4',
        name = "1000 Genomes Europe"
    )

    data = [trace1, trace2, trace3, trace4]

    layout = go.Layout(
        showlegend=True,
        legend_title="Populations",
        margin=dict(l=20, r=20, t=20, b=20),
        width=1200,
         
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        xaxis=dict(
            automargin = True,
            domain=[0, 0.24],
            titlefont=dict(size=20),
            title="Number of motif copies",
            range=[minx, maxx]
        ),
        yaxis=dict(
            automargin = True,
            title_text="Count in a population",
            titlefont=dict(size=20),
            showline=True
        ),
        xaxis2=dict(
            domain=[0.25, 0.49],
            anchor='y2',
            
            range=[minx, maxx]
        ),
        yaxis2=dict(
            anchor='x2',
        ),
        xaxis3=dict(
            domain=[0.50, 0.74],
            anchor='y3',
           
            range=[minx, maxx]
        ),
        yaxis3=dict(
            anchor='x3'
        ),
        xaxis4=dict(
            domain=[0.75, 1.0],
            anchor='y4',
             
            range=[minx, maxx]
        ),
        yaxis4=dict(
            anchor='x4'
        ))

    plotly_plot_json_datab = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)
    plotly_plot_json_layoutb = json.dumps(layout, cls=plotly.utils.PlotlyJSONEncoder)
    return plotly_plot_json_datab, plotly_plot_json_layoutb

def GetFreqPlotHg38(freq_dist):
    if freq_dist is None: return dict(), dict()
    data = []

    for cohort in freq_dist.groups.keys():
        items = freq_dist.get_group(cohort)
        trace = go.Bar(
            x=items['copies'],
            y=items['percentage'],
            name = cohort
        )
        data.append(trace)
    
    layout = go.Layout(
        width = 1200,
        barmode = 'group',     
              
        xaxis=dict(
            automargin = True,
            titlefont=dict(size=20),
            title="Number of motif copies",
        ),
        yaxis=dict(
            automargin = True,
            title_text="Fraction in in a population (%)",
            titlefont=dict(size=20),
            showline=True
        )
    )
    plotly_plot_json_datab = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)
    plotly_plot_json_layoutb = json.dumps(layout, cls=plotly.utils.PlotlyJSONEncoder)
    return plotly_plot_json_datab, plotly_plot_json_layoutb

################ Get sequence info (specific to hg38) #######

def GetSeqDataAPI(repeat_id):
    seq_url = API_URL + '/allseq/?repeat_id=' + repeat_id
    print(f"Fetching sequence data from: {seq_url}")  # Debug statement
    resp = requests.get(seq_url)

    print(f"Raw sequence response: {resp.text}")  # Debug raw response
    print(f"Response status code: {resp.status_code}")  # Check HTTP status code

    # Check if the API call fails with a 404 or another error
    if resp.status_code == 404:
        print("No sequence data available, API returned 404.")  # Debug statement
        return "No sequence information currently available"

    # Handle empty response gracefully
    json_data = resp.json()
    if not json_data:
        print(f"No sequence data found for repeat_id: {repeat_id}")  # Debug statement
        return "No sequence information currently available"

    # Process the data if it's not empty
    df = pd.DataFrame.from_records(json_data)

    if df.empty:
        print(f"No sequence data found after fetching for repeat_id: {repeat_id}")  # Debug statement
        return "No sequence information currently available"
    
    # Pivot and return the sequence data
    agg_df = df.pivot_table(index='sequence', columns='population', values='frequency', aggfunc='first').reset_index()
    agg_df = agg_df.fillna(0)

    return agg_df.to_dict(orient='records')

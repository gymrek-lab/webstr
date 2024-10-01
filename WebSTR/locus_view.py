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

API_URL = os.environ.get("WEBSTR_API_URL", 'http://webstr-api.ucsd.edu')

seqbuf = 120
seqbreakline = 100

################ General functions ################
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def GetSTRInfo(str_query, genome_query, dbSTRPath, reffa):
    logging.debug(f"Getting STR info for str_query: {str_query}, genome_query: {genome_query}")
    print(f"Getting STR info for {str_query} and genome {genome_query}")  # Debug statement

    # Fetch STR metadata based on genome
    if genome_query is None or genome_query == "hg19":
        strinfo = GetSTRMetadataHg19(str_query, dbSTRPath)
    else:
        strinfo = GetSTRMetadataAPI(str_query)

    # Ensure valid STR info
    if strinfo is None:
        logging.warning(f"No STR info found for {str_query}")
        return None

    try:
        chrom = strinfo["chrom"]
        start = strinfo["start"]
        end = strinfo["end"]
        print(f"STR info: chrom={chrom}, start={start}, end={end}")  # Debug statement

        # Process STR sequence
        lflank = str(reffa[chrom][start - seqbuf:start]).upper()
        strseq = str(reffa[chrom][start:end]).upper()
        rflank = str(reffa[chrom][end:end + seqbuf]).upper()

        strinfo["seq"] = GetSTRSeqHTML(lflank, strseq, rflank)
        strinfo["motif"] = GetMotifAndComplement(strinfo["motif"])

        # Check if frequency data exists and process accordingly
        freq_dist = strinfo.get("freq_dist", pd.DataFrame())  # Default to empty DataFrame
        if isinstance(freq_dist, pd.DataFrame) and freq_dist.empty:
            print("Frequency data is empty.")  # Debug statement
            strinfo["freq_dist"] = pd.DataFrame()  # Ensure it’s an empty DataFrame

        print(f"STR Info after processing: {strinfo}")  # Debug statement

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
        if numchar % seqbreakline == 0: ret += "<br>"
    ret += "</font>"
    ret += '<font size="4" color="red"><b>'
    for i in range(len(strseq)):
        ret += strseq[i]
        numchar += 1
        if numchar % seqbreakline == 0: ret += "<br>"
    ret += "</font></b>"
    ret += '<font size="3" color="black">'
    for i in range(len(rflank)):
        ret += rflank[i]
        numchar += 1
        if numchar % seqbreakline == 0: ret += "<br>"
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

    # Check if the response is valid JSON
    try:
        repeat = json.loads(resp.text)
        print(f"Data received from API for repeat info: {repeat}")  # Debug statement
    except json.JSONDecodeError:
        print(f"ERROR: Could not decode JSON for repeat_id {repeat_id}")  # Debug statement
        return None

    if not isinstance(repeat, dict):
        print(f"Invalid data structure for repeat_id {repeat_id}: {type(repeat)}")  # Debug statement
        return None

    # Handle missing or None fields safely
    strinfo = {
        "chrom": repeat.get("chr", "Unknown"),
        "start": repeat.get("start", "Unknown"),
        "end": repeat.get("end", "Unknown"),
        "motif": repeat.get("motif", "Unknown"),
        "copies": repeat.get("copies", "Unknown"),
        "gene_name": repeat.get("gene_name", "Unknown"),
        "gene_desc": repeat.get("gene_desc", "No description available"),
        "crc_data": None,
        "seq_data": GetSeqDataAPI(repeat_id) or "No sequence information currently available",
        "freq_dist": GetFreqSTRInfoAPI(repeat_id) or pd.DataFrame()  # Return empty DataFrame if no frequency data
    }

    if repeat.get("total_calls") is not None:
        strinfo["crc_data"] = [
            repeat.get('total_calls', 0),
            repeat.get('frac_variable', 0),
            repeat.get('avg_size_diff', 0)
        ]

    print(f"Processed STR metadata for repeat_id {repeat_id}: {strinfo}")  # Debug statement
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

    # Parse the response
    try:
        freq_data = json.loads(resp.text)
        print(f"Data received from API for frequency: {freq_data}")  # Debug statement
    except json.JSONDecodeError:
        print(f"Error decoding JSON for frequency data for repeat_id {repeat_id}")  # Debug statement
        return pd.DataFrame()  # Return an empty DataFrame if JSON decoding fails

    if len(freq_data) == 0:  # Check if the frequency data list is empty
        print(f"No frequency data found for repeat_id: {repeat_id}")  # Debug statement
        return pd.DataFrame()  # Return an empty DataFrame if no frequency data

    # Process the frequency data into a DataFrame
    try:
        df = pd.DataFrame.from_records(freq_data)
        if df.empty:  # Use .empty to check if the DataFrame is empty
            print(f"No frequency data found after fetching for repeat_id: {repeat_id}")  # Debug statement
            return pd.DataFrame()  # Return an empty DataFrame if no data exists

        df["percentage"] = df["frequency"] * 100
        df["copies"] = df["n_effective"]
        grouped_df = df.groupby(by="population")
        print(f"Grouped frequency data: {grouped_df}")  # Debug statement
        return grouped_df

    except Exception as e:
        print(f"Error processing frequency data for repeat_id {repeat_id}: {e}")  # Debug statement
        return pd.DataFrame()  # Return an empty DataFrame on error


def GetFreqPlotHg38(freq_dist):
    # Ensure freq_dist is a valid DataFrame with groups
    if freq_dist is None or (isinstance(freq_dist, pd.DataFrame) and freq_dist.empty):
        print("No frequency data available or freq_dist is an empty DataFrame.")  # Debug statement
        return dict(), dict()  # Return empty plot data to avoid breaking the page

    # Debug print statements for understanding freq_dist content
    print(f"Type of freq_dist: {type(freq_dist)}")  # Debug statement
    if isinstance(freq_dist, pd.DataFrame):
        print(f"Content of freq_dist (head if DataFrame): {freq_dist.head()}")  # Debug statement
    else:
        print(f"Contents of freq_dist: {freq_dist}")  # Debug statement
    
    # Checking if freq_dist has 'groups'
    if not hasattr(freq_dist, 'groups'):
        print(f"freq_dist does not have 'groups'. Contents: {freq_dist}")  # Debug statement
        return dict(), dict()

    data = []
    for cohort in freq_dist.groups.keys():
        items = freq_dist.get_group(cohort)
        print(f"Processing cohort: {cohort}, items: {items}")  # Debug statement
        trace = go.Bar(
            x=items['copies'],
            y=items['percentage'],
            name=cohort
        )
        data.append(trace)

    layout = go.Layout(
        width=1200,
        barmode='group',
        xaxis=dict(
            automargin=True,
            titlefont=dict(size=20),
            title="Number of motif copies",
        ),
        yaxis=dict(
            automargin=True,
            title_text="Fraction in a population (%)",
            titlefont=dict(size=20),
            showline=True
        )
    )
    
    try:
        plotly_plot_json_datab = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)
        plotly_plot_json_layoutb = json.dumps(layout, cls=plotly.utils.PlotlyJSONEncoder)
    except Exception as e:
        print(f"Error generating plotly JSON: {e}")  # Debug statement
        return dict(), dict()

    return plotly_plot_json_datab, plotly_plot_json_layoutb




################ Get sequence info (specific to hg38) #######

def GetSeqDataAPI(repeat_id):
    seq_url = API_URL + '/allseq/?repeat_id=' + repeat_id
    print(f"Fetching sequence data from: {seq_url}")  # Debug statement
    resp = requests.get(seq_url)

    if resp.status_code == 404:
        print("No sequence data available, API returned 404.")  # Debug statement
        return "No sequence information currently available"

    json_data = resp.json()
    print(f"Data received from API for sequence: {json_data}")  # Debug statement

    if not json_data:
        print(f"No sequence data found for repeat_id: {repeat_id}")  # Debug statement
        return "No sequence information currently available"

    df = pd.DataFrame.from_records(json_data)
    if df.empty:
        print(f"No sequence data found after fetching for repeat_id: {repeat_id}")  # Debug statement
        return "No sequence information currently available"

    agg_df = df.pivot_table(index='sequence', columns='population', values='frequency', aggfunc='first').reset_index()
    agg_df = agg_df.fillna(0)
    return agg_df.to_dict(orient='records')

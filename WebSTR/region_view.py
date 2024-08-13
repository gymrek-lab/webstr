import os
import requests
import pandas as pd
import numpy as np
import re
from dbutils import *
from utils import *

API_URL = os.environ.get("WEBSTR_API_URL",'http://webstr-api.ucsd.edu')

def GetRegionData(region_query, region_genome, DbSTRPath=None):
    if region_genome == "hg19":
        return GetRegionDataHg19(region_query, DbSTRPath)
    else:
        return GetRegionDataAPI(region_query)

############## Functions for hg19 version ##############

def GetRegionDataHg19(region_query, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    region_query = CleanRegionQuery(region_query)
   
    ######## Get coordinates to search #########
    genebuf = 0.1 # increase region width by this much
    if region_query.find(":") < 0: # search is by gene
        the_attrib = "gene_name"
        if region_query[:3] == "ENS":
            the_attrib = "gene_id"
            region_query = region_query.split(".")[0]
        else: 
            the_attrib = "gene_name"
        gene_query_sql = ("select fe.seqid,min(fe.start),max(fe.end)"
                      " from features fe, newattrib at "
                      " where at.value='{}' and at.attrib='{}' and fe.id=at.id").format(region_query, the_attrib)
        gene_df = ct.execute(gene_query_sql).fetchall()
        if len(gene_df) == 0 or None in gene_df[0]:
            chrom = None
            start = None
            end = None
        else:
            chrom = "chr"+gene_df[0][0].replace("chr","")
            start = int(gene_df[0][1])
            end = int(gene_df[0][2])
            buf = int((end-start)*(genebuf))
            start = start-buf
            end = end+buf
    else:
        try:
            chrom = "chr"+region_query.split(":")[0]
            start = int(region_query.split(":")[1].split("-")[0])
            end = int(region_query.split(":")[1].split("-")[1])
        except:
            chrom, start, end = None, None, None
    if chrom is None: return pd.DataFrame({}) # No data found

    ########## Load STRs in this region ###########
    region_query_sql = ("select str.chrom,str.strid,str.motif,str.start,str.end,str.period,str.length"
                    " from"
                    " strlocmotif str"
                    " where str.chrom = '{}' and str.end >= {} and str.start <= {}").format(chrom, start, end)
    df = ct.execute(region_query_sql).fetchall()
    if len(df) == 0: return pd.DataFrame({}) # No data found 

    df_hg19 = pd.DataFrame.from_records(df)
    df_hg19.columns = ["chr","strid", "motif", "start","end","period","length"] 
    df_hg19["featuretype"] = "NA"
    df_hg19["chr"] = df_hg19["chr"].apply(lambda x: x.replace("chr",""))
    df_hg19["length"] = df_hg19["length"].round(2)
    df_hg19 = df_hg19[["chr","start","end","motif","period","length","strid","featuretype"]].sort_values("start")
    df_hg19.drop_duplicates(inplace=True)

    ########## Load heterozygosity data ###########
    het_data = GetHCalc(df_hg19.strid.unique(), DbSTRPath)
    df_hg19 = pd.merge(df_hg19, het_data, left_on='strid', right_on='str_id', how='left')

    ########## Load eSTR data ###########
    estr_data = GetestrCalc(df_hg19.strid.unique(), DbSTRPath)
    df_hg19 = pd.merge(df_hg19, estr_data, left_on='strid', right_on='str_id', how='left')

    ########## Clean up the df and return ######
    df_hg19 = df_hg19.replace(np.nan, '', regex=True)
    df_hg19.rename(columns={'chrom':'chr', 'str.start':'start', 'str.end': 'end'}, inplace=True)

    return df_hg19


# set up for hg19 eSTR badges
def createret(thecolor, betav, tissue, gene):
    ret = f'<span class="badge popover-badge" data-toggle="popover" data-content="{tissue} ({gene})" style="background-color:{thecolor}; cursor: pointer;">{str(betav)}</span>'
    return ret


def GetestrHTML(df):
    df2 = pd.DataFrame(np.array(df).reshape(-1, 2), columns=list("TR"))
    t1 = df2["R"].str.split(pat=":", n=-1, expand=True)
    df2['thtml'] = None
    delim = [', ', ';']
    nrows = t1.shape[0]
    for i in range(nrows):
        ret = '<h5>'
        df2t = df2.iloc[i]
        t2 = df2t["R"].split(":")
        for j in range(len(t2)):
            t2t = t1.iloc[i, j]
            t2 = re.split(r'(?:' + '|'.join(delim) + r')', t2t)
            t2num = float(t2[1])
            t2tissue = t2[0]
            t2gene = t2[2]
            color = "black"
            if 'Adipose' in t2tissue:
                color = "darkorange"
            elif 'Artery-Aorta' in t2tissue:
                color = "salmon"
            elif 'Artery-Tibial' in t2tissue:
                color = "red"
            elif 'Brain_Caud' in t2tissue:
                color = "lemonchiffon"
            elif 'Brain_Cere' in t2tissue:
                color = "yellow"
            elif 'Cells' in t2tissue:
                color = "skyblue"
            elif 'Esophagus-Mucosa' in t2tissue:
                color = "sienna"
            elif 'Esophagus-Muscularis' in t2tissue:
                color = "burlywood"
            elif 'Heart' in t2tissue:
                color = "darkviolet"
            elif 'Lung' in t2tissue:
                color = "greenyellow"
            elif 'Muscle' in t2tissue:
                color = "mediumslateblue"
            elif 'Nerve' in t2tissue:
                color = "gold"
            elif 'Skin-Not' in t2tissue:
                color = "blue"
            elif 'Skin-Sun' in t2tissue:
                color = "cornflowerblue"
            elif 'Thyroid' in t2tissue:
                color = "green"
            elif 'WholeBlood' in t2tissue:
                color = "purple"

            ret += f'<span class="badge popover-badge" data-toggle="popover" data-content="{t2tissue} ({t2gene})" style="background-color:{color}; cursor: pointer;">{t2num}</span>'
        
        ret += '&nbsp;</h5>'
        df2.iloc[i, 2] = ret
    df2.columns = ["str_id", "tissues", "thtml"]
    return df2



def GetHCalc(strid,DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    stridt = tuple(strid)
    x1 = stridt
    if len(stridt) == 1: x1 = "('" + ''.join(stridt) + "')"
    stridt= x1
    gquery = (" select str_id, group_concat(name || ';' || H,':') Hvals from"
              " (select name,str_id,round(1-sum(fi*fi),1) H"
              " from"
              " (select cohort_id,str_id,copies,sum(cast(nvals as FLOAT)/cast(totvals as FLOAT)) as fi"
              " from"
              " (select af.cohort_id, af.str_id, (end-start+1+af.length)/period copies,sum(nvals) nvals, tvals.totvals from"
              " allelefreq af,"
              " strlocmotif strm,"
              " (select cohort_id,str_id,sum(nvals) totvals"
              " from allelefreq"
              " where str_id in  {} "
              " group by str_id,cohort_id"
              " order by str_id,cohort_id) tvals"
              " where af.str_id = strm.strid"
              " and af.str_id in  {} "
              " and af.str_id = tvals.str_id"
              " and af.cohort_id = tvals.cohort_id"
              " group by af.cohort_id, af.str_id, copies)"
              " group by cohort_id,str_id,copies) sum1,"
              " COHORTS co"
              " where sum1.cohort_id = co.cohort_id"
              " group by str_id , name ) group by str_id").format(stridt,stridt)
    df = ct.execute(gquery).fetchall()
    df2 = GetHvalSeqHTML2(df)
    return df2

def GetHvalSeqHTML2(df):
    df2 = pd.DataFrame(np.array(df).reshape(-1,2), columns= list("TR"))
    t1 = df2["R"].str.split(pat=":",n= -1, expand = True)
    df2['ret'] = None
    nrows =t1.shape[0]
    delim = [', ',';']
    for i in range(nrows):
        ret = '<h5>'
        t2ta = df2.iloc[i]
        t2 = t2ta["R"].split(":")
        for j in range(len(t2)):
            t2t = t1.iloc[i,j]
            t2=re.split(r'(?:' + '|'.join(delim) + r')', t2t )
            t2num = float(t2[1])
            if t2num == 0.0: 
                #thecolor = "gray"
                thecolor = "label-default"
            elif t2num <= 0.1:
                #thecolor = "blue"
                thecolor = "label-info"
            elif t2num <= 0.5:
                #thecolor = "navy"
                thecolor = "label-primary"
            else:
                thecolor = "label-danger"
            ret += '<span class="label ' + thecolor + '">' + t2[0] + ':' + t2[1] + '</span>'
            ret += '&nbsp;'
        ret += '</h5>'
        df2.iloc[i,2] = ret
    df2.columns = ["str_id","Hvals","Hhtml"]
    return df2            

def GetestrCalc(strid,DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    stridt = tuple(strid)
    x1 = stridt
    if len(stridt) == 1: x1 = "('" + ''.join(stridt) + "')"
    stridt= x1
    gquery = ("select strid str_id,  group_concat(tissue || ';' || round(beta,1) || ';' || genename, ':') tissues "
              "from estr_gtex2 estr, "
              "tissues ti, "
              "strlocmotif str "
              "where estr.chrom = str.chrom "
              "and estr.strstart = str.start "
              "and estr.strend = str.end "
              "and estr.signif = 'True' "
              "and estr.tissue_cd = ti.tissue_cd "
              "and strid in {} group by strid ").format(stridt)
    df = ct.execute(gquery).fetchall()
    if len(df) > 0:
        df2 = GetestrHTML(df)
    else:
        df2 = pd.DataFrame(columns = ["str_id","tissues","thtml"])
    return df2

################## Functions for hg38 version ####################
def GetRegionDataAPI(region_query):
    region_query = CleanRegionQuery(region_query)
    if region_query.find(":") < 0: # Search is by gene
        if region_query.find("ENSG") == 0:
            strexp_url = API_URL + '/repeats/?ensembl_ids=' + region_query
        else:
            strexp_url = API_URL + '/repeats/?gene_names=' + region_query
    else:
        strexp_url = API_URL + '/repeats/?region_query=' + region_query        
    resp = requests.get(strexp_url)
    df_hg38 = pd.DataFrame.from_records(resp.json())
    df_hg38.rename(columns={'repeat_id': 'strid'}, inplace=True)
    return df_hg38


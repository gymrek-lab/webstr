import os
from dbutils import *
import pyfaidx
import requests
import json
import numpy as np

seqbuf = 120
seqbreakline = 100
API_URL = os.environ.get("WEBSTR_API_URL",'http://webstr-api.ucsd.edu')
#API_URL = 'http://0.0.0.0:5000'

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

def GetSTRInfo(strid, DbSTRPath, reffa):
    ct = connect_db(DbSTRPath).cursor()
    squery = ("select str.chrom, str.start, str.end,str.motif, (str.end-str.start+1)/str.period copies from strlocmotif str where str.strid = '{}'").format(strid)
    df = ct.execute(squery).fetchall()
    if len(df) == 0: return None, None, None, None
    chrom = df[0][0]
    start = int(df[0][1])
    end = int(df[0][2])
    motif = df[0][3]
    copies = df[0][4]
    lflank = str(reffa[chrom][start-seqbuf:start]).upper()
    strseq = str(reffa[chrom][start:end]).upper()
    rflank = str(reffa[chrom][end:end+seqbuf]).upper()
    seq = GetSTRSeqHTML(lflank,strseq,rflank)
    return chrom, start, end, motif, copies, seq

def GetSTRInfoAPI(repeat_id, reffa):
    
    repeat_url = API_URL + '/repeatinfo/?repeat_id=' + repeat_id 
    
    resp = requests.get(repeat_url)
    repeat = json.loads(resp.text)


    #if len(df) == 0: return None, None, None, None
    chrom = repeat['chr']
    start =  repeat['start']
    end = repeat['end']
    motif = repeat['motif']
    copies = repeat['copies']
    gene_name = repeat['gene_name']
    gene_desc = repeat['gene_desc']
    crc_data = []
    if repeat['total_calls'] is not None:
        crc_data = [repeat['total_calls'], repeat['frac_variable'], repeat['avg_size_diff']]


    lflank = str(reffa[chrom][start-seqbuf:start]).upper()
    strseq = str(reffa[chrom][start-1:end]).upper()
    rflank = str(reffa[chrom][end:end+seqbuf]).upper()
    seq = GetSTRSeqHTML(lflank,strseq,rflank)
    return chrom, start, end, seq, gene_name, gene_desc, motif, copies, crc_data

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
    return df

def GetMutInfo(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select mut.est_logmu_ml, mut.est_beta_ml, mut.est_pgeom_ml, mut.up, mut.down, mut.p, mut.zscore_1, mut.zscore_2"
              " from mutrates mut where mut.str_id = '{}'").format(strid)
    df = ct.execute(gquery).fetchall()
    return df

def GetImputationInfo(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select imp.loo_concordance,imp.loo_r,imp.wgs_eur_concordance,imp.wgs_eur_r,imp.wgs_afr_concordance,imp.wgs_afr_r,"
              " imp.wgs_eas_concordance,imp.wgs_eas_r"
              " from locstat imp where imp.str_id = '{}'").format(strid)
    df = ct.execute(gquery).fetchall()
    return df

def GetFreqSTRInfo(strid,DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select cohort_id, (end-start+1+af.length)/period copies,sum(nvals) nvals from"
              " allelefreq af,"
              " strlocmotif strm"
              " where af.str_id = strm.strid"
              " and str_id = '{}' "
              " group by cohort_id, copies").format(strid)
    df = ct.execute(gquery).fetchall()
    return df

def GetFreqSTRInfoAPI(repeat_id):
    repeat_url = API_URL + '/allfreqs/?repeat_id=' + repeat_id 
    
    resp = requests.get(repeat_url)
    df = pd.DataFrame.from_records(json.loads(resp.text))

    if not df.empty:
        df["percentage"] = df["frequency"] * 100
        df["copies"] = df["n_effective"]
        grouped_df = df[["population", "copies", "percentage"]].groupby(by="population")
    

        return grouped_df
    else:
        return None
    #print(pd.DataFrame(np.array(grouped_df).reshape(-1,5), columns = list("abcdf")))

    #return grouped_df

def GetImputationAlleleInfo(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select al.allele, al.r2, al.pval"
              " from allelstat al where al.str_id = '{}'").format(strid)
    df = ct.execute(gquery).fetchall()
    return df


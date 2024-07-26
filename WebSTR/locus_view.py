import os
from dbutils import *
import pyfaidx
import requests
import json
import numpy as np
import pandas as pd


seqbuf = 120
seqbreakline = 100

API_URL = os.environ.get("WEBSTR_API_URL",'http://webstr-api.ucsd.edu')

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





def GetSeqDataAPI(repeat_id):
   
    seq_url = API_URL + '/allseq/?repeat_id=' + repeat_id
    print(seq_url)
    resp = requests.get(seq_url)
    
    # Create a DataFrame from the API response
    df = pd.DataFrame.from_records(resp.json())
    
    # Aggregating frequencies by sequence
    # This pivots the dataframe to have sequences as rows and populations as columns with frequencies as values
    agg_df = df.pivot_table(index='sequence', columns='population', values='frequency', aggfunc='first').reset_index()
    
    # Fill NaN values with 0 (or any other appropriate value indicating no frequency)
    agg_df = agg_df.fillna(0)
    
    # Convert the aggregated DataFrame back to a dictionary for rendering
    seq_data = agg_df.to_dict(orient='records')
    
    # Optional: If you still want to sort by sequence length and then by a population attribute, you'll need a different approach
    # since populations are now columns. One way to handle this might be to sort by sequence length directly:
    seq_data = sorted(seq_data, key=lambda x: len(x['sequence']))

    return seq_data



   

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


def GetImputationAlleleInfo(strid, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    gquery = ("select al.allele, al.r2, al.pval"
              " from allelstat al where al.str_id = '{}'").format(strid)
    df = ct.execute(gquery).fetchall()
    return df

def GetFreqPlotlyJSON2(freq_dist):
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

def GetFreqPlot(freq_dist):
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


"""
Functions for making gene plots on the region page
"""

import json
import os
import plotly
import plotly.graph_objs as go
import requests
from dbutils import *
from utils import *

# Global plotting variables
GENEBUFFER = 0.1
EXON_WIDTH = 0.3
GENE_WIDTH = 0.03
GENE_COLOR = "black"

# Set up API
API_URL = os.environ.get("WEBSTR_API_URL",'http://webstr-api.ucsd.edu')

def GetGeneShapes(region_query, region_genome, DbSTRPath=None):
    region_query = CleanRegionQuery(region_query)
    # Extract gene features
    if region_genome == "hg19":
        gene_features = ExtractGeneFeaturesHg19(region_query, DbSTRPath)
    else:
        gene_features = ExtractGeneFeaturesAPI(region_query)

    # Get shapes
    min_start = 99999999999
    max_end = 0
    shapes = []
    num_genes = len(gene_features)
    for i in range(num_genes):
        gene = gene_features[i]
        if (gene['start'] < min_start): min_start = gene['start']
        if (gene['end'] > max_end): max_end = gene['end'] 
        buf = int((max_end-min_start)*(GENEBUFFER))
        max_end = max_end + buf
        min_start = min_start - buf

        shape = {
            "type": "rect",
            "x0": gene['start'],
            "x1": gene['end'],
            "y0": (i+1)-GENE_WIDTH/2,
            "y1": (i+1)+GENE_WIDTH/2,
            "fillcolor": GENE_COLOR,
            "line": {"width": 0}
            }
        shapes.append(shape)
        # Then put each feature (exons only)
        for exon in gene['exons']:
            shape = {
                "type": "rect",
                "x0": exon['start'],
                "x1": exon['end'],
                "y0": (i+1)-EXON_WIDTH/2,
                "y1": (i+1)+EXON_WIDTH/2,
                "fillcolor": GENE_COLOR,
                "line": {"width": 0}
                }
            shapes.append(shape)
    if num_genes != 0:
        trace = go.Scatter(
            x = [gene['start'] for gene in gene_features],
            y = [(i+1+EXON_WIDTH) for i in range(num_genes)],
            mode = "text",
            hoverinfo="none",
            textposition='middle right',
            text = [GetGeneText(gene['name'], gene['strand']) for gene in gene_features],
            textfont=dict(
                family='sans serif',
                size=20,
                color='black')
        )
    else:
        trace = []
    return trace, shapes, num_genes

def GetSTRColor(period):
    colors = ["gray","red","gold","blue","purple","green","magenta", "pink", "yellow"]
    color_index = int(period)-1
    if color_index <= len(colors):
        return colors[int(period)-1]
    else:
        return "black"

def GetGenePlotlyJSON(region_data, gene_trace, gene_shapes, numgenes):
    # Draw gene info
    region_data2 = region_data
    
    #chrom = region_data2["chrom"].values[0].replace("chr","")
    chr = region_data2["chr"].values[0].replace("chr","")
    
    # Get points for each STR
    trace1 = go.Scatter(
        x = (region_data2["start"]+region_data2["end"])/2,
        y = [0]*region_data2.shape[0],
        mode="markers",
        marker=dict(size=10, color=region_data2["period"].apply(lambda x: GetSTRColor(x)), line=dict(width=2)),
        text=region_data2.apply(lambda x: x["chr"]+":"+str(x["start"]) + " ("+x["motif"]+")", 1),
        hoverinfo='text'
    )
    plotly_data = [trace1, gene_trace]
    plotly_layout= go.Layout(
        height=300+50*numgenes,
        hovermode= 'closest',
        showlegend= False,
        #legend=dict(orientation="h"),
        shapes=gene_shapes,
        xaxis=dict(
            title="Position (chr%s)"%chr ,
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=True,
            showticklabels=True,
            tickformat = '.0f'
        ),
        yaxis=dict(
            fixedrange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            showticklabels=False
        )
    )
    plotly_plot_json = json.dumps(plotly_data, cls=plotly.utils.PlotlyJSONEncoder) 
    plotly_layout_json = json.dumps(plotly_layout, cls=plotly.utils.PlotlyJSONEncoder)
    return plotly_plot_json, plotly_layout_json

def GetGeneText(gene_name, strand):
    txt = "<i>"
    if strand == "+":
        txt += gene_name + " " + "&#8594;"
    else:
        txt += "&#8592;" + " " + gene_name
    txt += " "*3
    txt += "</i>"
    return txt

######### Functions for hg19 versions ##########
def ExtractGeneFeaturesHg19(region_query, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    genes = []
    if region_query.find(":") < 0: # Search is by gene
        if region_query[:3] == "ENS":
            gene_query = ("select at.value from newattrib at where at.id='{}' and at.attrib='gene_name'").format(region_query.split(".")[0])
            gene_df = ct.execute(gene_query).fetchall()
            genes = [item[0] for item in gene_df]
        else:
            genes.append(region_query)
    else: # Search is by region
        chrom = "chr"+region_query.split(":")[0].replace("chr","")
        start = int(region_query.split(":")[1].split("-")[0])
        end = int(region_query.split(":")[1].split("-")[1])
        gene_query = ("select value from newattrib where attrib = 'gene_name' and id in (select id from features where seqid='{}' and start>={} and end<={} group by id) ").format(chrom, start, end)
        gene_df = ct.execute(gene_query).fetchall()
        genes = list(set([item[0] for item in gene_df]))

    gene_features = []
    for g in genes:
        # Get exons
        feature_query = ("select fe.id,fe.start,fe.end,fe.strand from features fe, newattrib at where at.attrib='gene_name' and at.value='{}' and fe.id=at.id").format(g)
        feature_df = ct.execute(feature_query).fetchall()
        # Get gene info
        gene_info = {"name": g, "exons": [], \
                     "start": min([int(item[1]) for item in feature_df]), \
                     "end": max([int(item[2]) for item in feature_df]), \
                     "strand": feature_df[0][3]}
        # Fill in exon info
        for f in feature_df:
            if "CDS" in f[0]: continue
            gene_info["exons"].append({"start": f[1], "end": f[2]})
        gene_features.append(gene_info)
    return gene_features

######### Functions for hg38 versions ##########
def ExtractGeneFeaturesAPI(region_query):
    if region_query.find(":") < 0:
        if region_query.find("ENSG") == 0:
            gene_url = API_URL + '/genefeatures/?ensembl_ids=' + region_query
        else:
            gene_url = API_URL + '/genefeatures/?gene_names=' + region_query
    else:
        gene_url = API_URL + '/genefeatures/?region_query=' + region_query
    resp = requests.get(gene_url)
    gene_features = json.loads(resp.text)
    if gene_features is None: return []
    return gene_features


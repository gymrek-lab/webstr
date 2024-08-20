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
    chrom = region_data["chr"].values[0].replace("chr","")
    
    # Get points for each STR
    trace1 = go.Scatter(
        x = (region_data["start"]+region_data["end"])/2,
        y = [0]*region_data.shape[0],
        mode="markers",
        marker=dict(size=10, color=region_data["period"].apply(lambda x: GetSTRColor(x)), line=dict(width=2)),
        text=region_data.apply(lambda x: x["chr"]+":"+str(x["start"]) + " ("+x["motif"]+")", 1),
        hoverinfo='text'
    )
    plotly_data = [trace1, gene_trace]
    plotly_layout= go.Layout(
        height=300+50*numgenes,
        hovermode= 'closest',
        showlegend= False,
        shapes=gene_shapes,
        xaxis=dict(
            title="Position (chr%s)"%chrom ,
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
   """
   Extracts gene features from the API based on the given region query or gene names/Ensembl IDs.
   Args:
       region_query (str): The region query or gene names/Ensembl IDs.
   Returns:
       list: A list of gene features extracted from the API.
   """
   # Threshold distance to merge nearby exons
   merge_distance_threshold = 50


   # Check if the query is for a region
   if ":" in region_query:
       # Use the /repeats endpoint to get all STRs in the region
       gene_url = f"{API_URL}/repeats?region_query={region_query}"
      
       # Fetch data from the API
       resp = requests.get(gene_url)
      
       # Check if the API call was successful
       if resp.status_code != 200:
           print(f"Error: API call failed with status code {resp.status_code}")
           return []


       try:
           str_data = resp.json()
       except json.JSONDecodeError:
           print("Failed to parse JSON response")
           return []
       # Check if the response is empty
       if not str_data: 
           return []


       gene_dict = {}
       for entry in str_data:
           gene_id = entry.get("ensembl_id")
           gene_name = entry.get("gene_name")
           strand = entry.get("strand")
           start = entry.get("start")
           end = entry.get("end")
          
           # Check if gene_id and gene_name are present
           if gene_id and gene_name: 
               if gene_id not in gene_dict:
                   gene_dict[gene_id] = {
                       "name": gene_name,
                       "strand": strand,
                       "start": start,
                       "end": end,
                       "exons": [{"start": start, "end": end}]
                   }
               else:
                   # Update start and end boundaries
                   gene_dict[gene_id]["start"] = min(gene_dict[gene_id]["start"], start)
                   gene_dict[gene_id]["end"] = max(gene_dict[gene_id]["end"], end)
                   gene_dict[gene_id]["exons"].append({"start": start, "end": end})


       # Merge nearby exons within the merge_distance_threshold
       for gene_id, gene_info in gene_dict.items():
           merged_exons = []
           sorted_exons = sorted(gene_info["exons"], key=lambda x: x["start"])
          
           current_exon = sorted_exons[0]
           for next_exon in sorted_exons[1:]:
               if next_exon["start"] - current_exon["end"] <= merge_distance_threshold:
                   # Extend the current exon
                   current_exon["end"] = max(current_exon["end"], next_exon["end"])
               else:
                   # Append the current exon and move to the next
                   merged_exons.append(current_exon)
                   current_exon = next_exon
           # Append the last exon
           merged_exons.append(current_exon)


           # Update the gene info with merged exons
           gene_info["exons"] = merged_exons


       # Convert gene_dict to a list of gene features
       gene_features = list(gene_dict.values())
      
   else:
       # If it's not a region query, use the original logic for gene names and Ensembl IDs
       if region_query.startswith("ENSG"):
           gene_url = f"{API_URL}/genefeatures/?ensembl_ids={region_query}"
       else:
           gene_url = f"{API_URL}/genefeatures/?gene_names={region_query}"
      
       resp = requests.get(gene_url)
       gene_features = json.loads(resp.text)
       if gene_features is None:
           return []
  
   # Print the merged gene features
   for gene_info in gene_features:
       print(f"Gene: {gene_info['name']}, Start: {gene_info['start']}, End: {gene_info['end']}, Exons: {gene_info['exons']}")


   return gene_features


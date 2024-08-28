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
GENEBUFFER = 0.06
EXON_WIDTH = 0.3
GENE_WIDTH = 0.03
GENE_COLOR = "black"


API_URL = os.environ.get("WEBSTR_API_URL",'http://localhost:5000')

def GetGeneShapes(region_query, region_genome, DbSTRPath=None):
    region_query = CleanRegionQuery(region_query)
    
    # Extract gene features
    if region_genome == "hg19":
        gene_features = ExtractGeneFeaturesHg19(region_query, DbSTRPath)
    else:
        gene_features = ExtractGeneFeaturesAPI(region_query)

    if len(gene_features) == 0:
        print(f"No gene features found for query: {region_query}")
        return [], [], [], 0  # Return empty lists and 0 genes

    # Determine if the query is a region or a gene due to added functionality of partial drawings allowing redirect to gene search
    if ":" in region_query:
        region_split = region_query.split(':')
        min_start = int(region_split[1].split('-')[0])  
        max_end = int(region_split[1].split('-')[1])   
    else:
        # Handle gene queries or other cases
        min_start = min(gene['start'] for gene in gene_features)
        max_end = max(gene['end'] for gene in gene_features)

    buf = int((max_end - min_start) * GENEBUFFER)
    buffered_max_end = max_end + buf
    buffered_min_start = min_start - buf

    shapes = []
    annotations = []
    num_genes = len(gene_features)

    for i in range(num_genes):
        gene = gene_features[i]

        # Check if gene is partial
        is_partial_left = gene['start'] < min_start
        is_partial_right = gene['end'] > max_end

        # Adjust gene start and end to be within the buffered region
        gene_start = max(gene['start'], buffered_min_start)
        gene_end = min(gene['end'], buffered_max_end)
        shape = {
            "type": "rect",
            "x0": gene_start,
            "x1": gene_end,
            "y0": (i + 1) - GENE_WIDTH / 2,
            "y1": (i + 1) + GENE_WIDTH / 2,
            "fillcolor": GENE_COLOR,
            "line": {"width": 0}
        }
        shapes.append(shape)

        # Draw each exon that falls within the buffered region
        for exon in gene['exons']:
            exon_start = max(exon['start'], buffered_min_start)
            exon_end = min(exon['end'], buffered_max_end)
            if exon_start <= exon_end:  # Ensure valid exon range
                shape = {
                    "type": "rect",
                    "x0": exon_start,
                    "x1": exon_end,
                    "y0": (i + 1) - EXON_WIDTH / 2,
                    "y1": (i + 1) + EXON_WIDTH / 2,
                    "fillcolor": GENE_COLOR,
                    "line": {"width": 0}
                }
                shapes.append(shape)

        # Add left arrow if the gene extends left
        if is_partial_left:
            annotations.append(dict(
                x=gene_start,
                y=(i + 1 + GENE_WIDTH),
                xref="x",
                yref="y",
                text="&#8606;", 
                showarrow=False,
                font=dict(size=40, color="#4874e4"),
                hovertext=f"{gene['name']} start coordinate: {gene['start']}",
                hoverlabel=dict(bgcolor="white", font_size=16),
                xshift=-15  # Shift the arrow left
            ))

        # Add right arrow if the gene extends right
        if is_partial_right:
            annotations.append(dict(
                x=gene_end,
                y=(i + 1 + GENE_WIDTH),
                xref="x",
                yref="y",
                text="&#8608;",  # Right arrow symbol
                showarrow=False,
                font=dict(size=40, color="#4874e4"),
                hovertext=f"{gene['name']} end coordinate: {gene['end']}",
                hoverlabel=dict(bgcolor="white", font_size=16),
                xshift=15  # Shift the arrow right
            ))

    # Create the gene text trace
    if num_genes != 0:
        trace = go.Scatter(
            x=[max(buffered_min_start, gene['start']) for gene in gene_features],
            y=[(i + 1 + EXON_WIDTH) for i in range(num_genes)],
            mode="text",
            hoverinfo="none",
            textposition='middle right',
            
            text=[GetGeneText(gene['name'], gene['strand']) for gene in gene_features],
            textfont=dict(
                family='sans serif',
                size=20,
                color='black'
            )
        )
    else:
        trace = []

    return trace, shapes, annotations, num_genes


def GetSTRColor(period):
    colors = ["gray","red","gold","blue","purple","green","magenta", "pink", "yellow"]
    color_index = int(period)-1
    if color_index <= len(colors):
        return colors[int(period)-1]
    else:
        return "black"

def GetGenePlotlyJSON(region_data, gene_trace, gene_shapes, annotations, numgenes):
    chrom = region_data["chr"].values[0].replace("chr","")

    # Get points for each STR
    trace1 = go.Scatter(
        x=(region_data["start"] + region_data["end"]) / 2,
        y=[0] * region_data.shape[0],
        mode="markers",
        marker=dict(size=10, color=region_data["period"].apply(lambda x: GetSTRColor(x)), line=dict(width=2)),
        text=region_data.apply(lambda x: x["chr"] + ":" + str(x["start"]) + " (" + x["motif"] + ")", 1),
        hoverinfo='text'
    )

    plotly_data = [trace1, gene_trace]

    plotly_layout = go.Layout(
        height=300 + 50 * numgenes,
        hovermode='closest',
        showlegend=False,
        shapes=gene_shapes,
        annotations=annotations,
        xaxis=dict(
            title="Position (chr%s)" % chrom,
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=True,
            showticklabels=True,
            tickformat='.0f'
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
    # Default to empty string if gene name is None
    if gene_name is None:
        gene_name = ""

    # Initialize the text with italic tags
    txt = "<i>"

    
    # Add the gene name with directional arrow
    if strand == "+":
        txt += gene_name + " " + "&#8594;"
    else:
        txt += "&#8592;" + " " + gene_name

    
    # Close italic tags
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

    try:
        gene_features = json.loads(resp.text)
    except json.JSONDecodeError as e:
        print(f"Failed to decode")
        gene_features = None
    
    if gene_features is None:
        return []

    
    return gene_features

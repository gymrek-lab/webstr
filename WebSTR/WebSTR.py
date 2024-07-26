#!/usr/bin/env python3
"""
WebSTR database application
"""

import sys
import os
import argparse
from flask import Flask, redirect, render_template, request, session, url_for, jsonify
import pandas as pd
import numpy as np

from locus_view import *
from region_view import *
from gene_plots import *

#from dash_graphs import add_dash_graphs_to_flask_server

API_URL = os.environ.get("WEBSTR_API_URL",'http://webstr-api.ucsd.edu')

#################### Database paths ###############
PLATFORM = "snorlax" # or AWS
BASEPATH =  "/storage/resources/dbase/human/"
if PLATFORM == "snorlax":
    DbSTRPath = BASEPATH
    RefFaPath_hg19 = BASEPATH + "hg19/hg19.fa"
    RefFaPath_hg38 = BASEPATH + "hg38/hg38.fa"
else:
    sys.stderr.write("Could not locate database files\n")
    sys.exit(1)

#################### Set up flask server ###############
server = Flask(__name__)
server.secret_key = 'dbSTR' 
#add_dash_graphs_to_flask_server(server)
 
#################### Render region page ###############
@server.route('/search')
def search():
    region_genome = request.args.get('genome')
    region_query = request.args.get('query').upper()
    region_data = GetRegionData(region_query, region_genome, DbSTRPath)
    if region_data.shape[0] == 0: return render_template('view2_nolocus.html')
    gene_trace, gene_shapes, numgenes = GetGeneShapes(region_query, region_genome, DbSTRPath)
    plotly_plot_json, plotly_layout_json = GetGenePlotlyJSON(region_data, gene_trace, gene_shapes, numgenes)
    return render_template('view2.html',
                           table = region_data.to_records(index=False),
                           graphJSON = plotly_plot_json, layoutJSON = plotly_layout_json,
                           chrom = region_data["chr"].values[0].replace("chr",""),
                           strids = list(region_data["strid"]),
                           genome = region_genome) 

#################### Render locus page ###############
@server.route('/locus')
def locusview():
    str_query = request.args.get('repeat_id')
    genome_query = request.args.get('genome')

    # Extract STR info
    if genome_query == "hg38":
        reffa = pyfaidx.Fasta(RefFaPath_hg38)
    else:
        reffa = pyfaidx.Fasta(RefFaPath_hg19)
    strinfo = GetSTRInfo(str_query, genome_query, DbSTRPath, reffa)
    if strinfo is None: return render_template('view2_nolocus.html')
    
    # Plotting info
    plotly_plot_json_datab, plotly_plot_json_layoutb = GetFreqPlotlyJSON(genome_query, strinfo["freq_dist"])

    # Render
    return render_template('locus.html', strid=str_query,
                           graphJSONx=plotly_plot_json_datab, graphlayoutx=plotly_plot_json_layoutb, 
                           chrom=strinfo["chrom"], start=strinfo["start"], end=strinfo["end"], 
                           strseq=strinfo["seq"], gene_name=strinfo["gene_name"], gene_desc=strinfo["gene_desc"],
                           motif=strinfo["motif"], copies=strinfo["copies"], crc_data=strinfo["crc_data"],
                           estr=strinfo["gtex_data"], mut_data=strinfo["mut_data"],
                           imp_data=strinfo["imp_data"], imp_allele_data=strinfo["imp_allele_data"],
                           seq_data=strinfo["seq_data"])

#################### Render other HTML pages ###############

#### Static pages ####
@server.route('/')
def dbSTRHome():
    return render_template('homepage.html')

@server.route('/faq')
def dbSTRFAQ():
    return render_template("faq.html")

@server.route('/contact')
def dbSTRContact():
    return render_template("contact.html")

@server.route('/about')
def dbSTRAbout():
    return render_template("about.html")

@server.route('/downloads')
def dbSTRDownloads():
    return render_template("downloads.html")

@server.route('/terms')
def dbSTRTerms():
    return render_template("terms.html")

#### Predefined locus set pages #####
@server.route('/pathogenic')
def dbSTRpathogenic():
    return render_template("pathogenic.html")

@server.route('/GWAS')
def dbSTRGWAS():
    return render_template("GWAS.html")

#### CRC research ####
@server.route('/crc_research')
def graphs():    
    return render_template("crc_research.html", api_url = API_URL)

#### Error pages #####
@server.route('/url')
def my_method():
    try:
        call_method_that_raises_exception()
    except Exception as e:
        render_template("500.html", error= str(e))

@server.errorhandler(404)
def internal_server_error(error):
    server.logger.error('Server Error: %s', (error))
    return render_template('500.htm', emsg = error), 404

@server.errorhandler(Exception)
def unhandled_exception(e):
    server.logger.error('Unhandled Exception: %s', (e))
    return render_template('500.html', emsg = e), 500
 
#################### Set up and run the server ###############
def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--host", help="Host to run app", type=str, default="0.0.0.0")
    parser.add_argument("--port", help="Port to run app", type=int, default=int(os.environ.get("FLASK_PORT", "5000"))) 

    args = parser.parse_args()

    # FLASK_DEBUG is not mandatory but can be handy to configure debugging in a container
    FLASK_DEBUG = os.environ.get("FLASK_DEBUG", 0) == "1"
    if FLASK_DEBUG:
        server.config.update(
            TEMPLATES_AUTO_RELOAD=True
        )
    server.run(debug=FLASK_DEBUG, host=args.host, port=args.port)

if __name__ == '__main__':
    main()
 

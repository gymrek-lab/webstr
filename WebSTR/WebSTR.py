#!/usr/bin/env python3
"""
WebSTR database application
"""

import argparse
#import dash
from flask import Flask, redirect, render_template, request, session, url_for
#from dash.dependencies import Output, Input, State
from collections import deque
import pandas as pd
import numpy as npa
import json
from textwrap import dedent as d
import sys
import os

#from locus_view_dash import *
from locus_view import *
from region_view import *

#################### Database paths ###############
PLATFORM = "snorlax" # or AWS
BASEPATH =  "/storage/resources/dbase/human/"
#os.environ['DATAPATH']
if PLATFORM == "snorlax":
    #BasePath = "/storage/resources/dbase/dbSTR/SS1/" # TODO this is allele freq. not used now
    #DbSTRPath = "/storage/resources/dbase/dbSTR/"
    #RefFaPath_hg19 = "/storage/resources/dbase/human/hg19/hg19.fa"
    DbSTRPath = BASEPATH
    RefFaPath_hg19 = BASEPATH + "hg19/hg19.fa"
    RefFaPath_hg38 = BASEPATH + "hg38/hg38.fa"
elif PLATFORM == "AWS":
    #BasePath = ""
    DbSTRPath = ""
    RefFaPath_hg19 = "" # TODO
else:
    sys.stderr.write("Could not locate database files\n")
    sys.exit(1)

#################### Set up flask server ###############
server = Flask(__name__)
server.secret_key = 'dbSTR' 

#################### Render locus page ###############
#app = dash.Dash(__name__, server=server, url_base_pathname='/dashapp')
#app.config['suppress_callback_exceptions']=True
#SetupDashApp(app)
#
#@app.callback(dash.dependencies.Output('field-dropdown','value'),
#              [dash.dependencies.Input('url', 'href')])
#def main_display_page(href): return display_page(href)
#
#@app.callback(Output('table2', 'rows'), [Input('field-dropdown', 'value')])
#def main_update_table(user_selection): return update_table(user_selection, BasePath)
#
#@app.callback(Output('STRtable', 'rows'), [Input('field-dropdown', 'value')])
#def main_getdata(user_selection): return getdata(user_selection, BasePath)
#
#@app.callback(Output('Main-graphic','figure'),
#              [Input('table2','rows')])
#def main_update_figure(rows): return update_figure(rows)

#################### Render region page ###############

@server.route('/search')
def search():
    region_queryGenome = request.args.get('genome')
    region_queryOrg = request.args.get('query')
    region_query = region_queryOrg.upper()

    if (region_queryGenome == 'hg19'):
        region_data = GetRegionData(region_query, DbSTRPath)
        
        if region_data.shape[0] > 0:
            strs_id = region_data.strid.unique()

            H_data = GetHCalc(strs_id,DbSTRPath)
            estr_data = GetestrCalc(strs_id,DbSTRPath)
            Regions_data = pd.merge(region_data, H_data, left_on='strid', right_on = 'str_id')
            Regions_data = pd.merge(Regions_data, estr_data, left_on='strid', right_on = 'str_id', how='left')
            Regions_data = Regions_data.replace(np.nan, '', regex=True)
            # Get the STRs on the plotly graph
            Regions_data.rename(columns = {'chrom':'chr', 'str.start':'start', 'str.end': 'end'}, inplace = True)
            #chrom = Regions_data["chr"].values[0].replace("chr","")
            gene_trace, gene_shapes, numgenes, min_gene_start, max_gene_end = GetGeneShapes(region_query, DbSTRPath)
            region_data2 = Regions_data
            #if (max_gene_end) > 0:
            #    region_data2 =  GetRegionData(region_data["chr"].values[0] + ":" + str(min_gene_start) + "-" + str(max_gene_end), DbSTRPath)

            plotly_plot_json, plotly_layout_json = GetGenePlotlyJSON(region_data2, gene_trace, gene_shapes, numgenes)

            return render_template('view2.html',
                                table = Regions_data.to_records(index=False),
                                graphJSON = plotly_plot_json, layoutJSON = plotly_layout_json,
                                chrom = region_data["chr"].values[0].replace("chr",""),
                                strids = list(Regions_data["strid"]),
                                genome = region_queryGenome) 
        else:
            return render_template('view2_nolocus.html')
    else:
        # Use the API, because hg38 is requested
        region_data_hg38 = GetRegionDataAPI(region_query)
        if region_data_hg38.shape[0] > 0:
            gene_trace_hg38, gene_shapes_hg38, numgenes_hg38, min_gene_start_hg38, max_gene_end_hg38 = GetGeneGraph(region_query)
            plotly_plot_json_hg38, plotly_layout_json_hg38 = GetGenePlotlyJSON(region_data_hg38, gene_trace_hg38, gene_shapes_hg38, numgenes_hg38)
            return render_template('view2.html',
                                    table = region_data_hg38.to_records(index=False),
                                    graphJSON = plotly_plot_json_hg38, layoutJSON = plotly_layout_json_hg38,
                                    chrom = region_data_hg38["chr"].values[0].replace("chr",""),
                                    strids = list(region_data_hg38["repeat_id"]),
                                    genome = region_queryGenome)
        else:
            return render_template('view2_nolocus.html')



@server.route('/locus')
def locusview():
    str_query = request.args.get('repeat_id')
    genome_query = request.args.get('genome')
    mut_data = []
    imp_data = []
    gtex_data = []
    imp_allele_data = []
    freq_dist = []
    crc_data = []
    gene_name = ""
    gene_desc = ""
    motif = ""
    copies = ""
    plotly_plot_json_datab = dict()
    plotly_plot_json_layoutb = dict()

    if ((genome_query is None) or (genome_query == 'hg19')):
        reffa = pyfaidx.Fasta(RefFaPath_hg19)

        chrom, start, end, motif, copies, seq = GetSTRInfo(str_query, DbSTRPath, reffa)
        gtex_data = GetGTExInfo(str_query, DbSTRPath)
        mut_data = GetMutInfo(str_query, DbSTRPath)
        imp_data = GetImputationInfo(str_query, DbSTRPath)
        imp_allele_data = GetImputationAlleleInfo(str_query, DbSTRPath)
        freq_dist = GetFreqSTRInfo(str_query, DbSTRPath)
        if len(freq_dist) > 0:
            plotly_plot_json_datab, plotly_plot_json_layoutb = GetFreqPlotlyJSON2(freq_dist)
        
    elif (genome_query == 'hg38'):
        reffa = pyfaidx.Fasta(RefFaPath_hg38)
        chrom, start, end, seq, gene_name, gene_desc, motif, copies, crc_data = GetSTRInfoAPI(str_query, reffa)
        freq_dist = GetFreqSTRInfoAPI(str_query)
        if freq_dist:
            plotly_plot_json_datab, plotly_plot_json_layoutb = GetFreqPlot(freq_dist)
        
    if len(mut_data) != 1: mut_data = None
    else:
        mut_data = list(mut_data[0])
        mut_data[0] = 10**mut_data[0]
    if len(imp_data) != 1: imp_data = None
    else:
        imp_data = list(imp_data[0])

    if len(gtex_data) == 0: gtex_data = None
    if len(crc_data) == 0: crc_data = None
    if len(imp_allele_data) == 0: imp_allele_data = None

    
    return render_template('locus.html', strid=str_query,
                           graphJSONx=plotly_plot_json_datab,graphlayoutx=plotly_plot_json_layoutb, 
                           chrom=chrom.replace("chr",""), start=start, end=end, strseq=seq,
                           gene_name=gene_name, gene_desc=gene_desc,
                           estr=gtex_data, mut_data=mut_data, motif=motif, copies=copies, crc_data = crc_data,
                           imp_data=imp_data, imp_allele_data=imp_allele_data,freq_dist=freq_dist)

#################### Render HTML pages ###############
@server.route('/')
@server.route('/dbSTR')
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

@server.route('/url')
def my_method():
    try:
        call_method_that_raises_exception()
    except Exception as e:
            render_template("500.html", error= str(e))

@server.route('/pathogenic')
def dbSTRpathogenic():
    return render_template("pathogenic.html")


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
    parser.add_argument("--port", help="Port to run app", type=int, default=5000)
    parser.add_argument("--ref-hg38", help="Address to hg38 ref genome", type=str, default="/storage/resources/dbase/human/hg38/hg38.fa")
    parser.add_argument("--ref-hg19", help="Address to hg19 ref genome", type=str, default="/storage/resources/dbase/human/hg19/hg19.fa")
    args = parser.parse_args()
    server.run(debug=False, host=args.host, port=args.port)

if __name__ == '__main__':
    main()

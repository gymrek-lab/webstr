import json
import logging
import os
import dash
import dash_bio as dashbio
from dash import dcc, html, Input, Output, callback
from dash.dependencies import Input, Output, State
from flask import render_template
import pandas as pd
import requests

API_URL = os.environ.get("WEBSTR_API_URL",'http://webstr-api.ucsd.edu')


def GetCrcExprRepeatLenCorrInfoAPI():
    """
    retrieves data from an API endpoint crc_expr_repeatlen_corr and returns it as a Pandas DataFrame.
    """

    repeat_url = API_URL + '/crc_expr_repeatlen_corr/'

    logging.warn(f"GetCrcExprRepeatLenCorrInfoAPI is loading data from {repeat_url}")
    
    resp = requests.get(repeat_url)
    data_frame = pd.read_json(resp.text)

    logging.warn(f"GetCrcExprRepeatLenCorrInfoAPI recieved {len(data_frame)} entities")

    return data_frame

DATASET_CRC_EXPR_LEN_CORR = {
            'label': 'STR',
            'dataframe': None,
            'dataage': None,
            'datasource': '',
            'dataprops': dict(
                #effect_size='log2_(L3i/L1)_signal_ratio',
                effect_size='coefficient',
                p='p_value_corrected',
                snp=None,
                gene='annotation',
                ylabel='FDR corrected p-value (-log10)',
                xlabel='Gene expression coefficient',
                height=500,
                title='Correlation between gene expression and STR length in CRC patients',
                genomewideline_value = 2,
                #annotation='annotation'
            ),
            'annotations': {}
    }

def load_crc_corr():
    data =  GetCrcExprRepeatLenCorrInfoAPI()
    annotations = {}

    for i, row in data.iterrows():
        data.at[i, 'annotation'] = f"{row['name']} {row['ensembl_id']} <br>REPEAT: {row['chr']}_{row['start']}" 
        annotations[data.at[i, 'annotation']] = data.iloc[i]

    logging.warn(f"padded annotations to {i} rows")
    DATASET_CRC_EXPR_LEN_CORR['dataframe'] = data
    DATASET_CRC_EXPR_LEN_CORR['annotations'] = annotations


@callback(
        Output('dashbio-default-volcanoplot', 'figure'),
        Input('default-volcanoplot-input', 'value'),
    )
def render_volcano_plot(effects):
    """Update rendering of data points upon changing x-value of vertical dashed lines."""

    if DATASET_CRC_EXPR_LEN_CORR['dataframe'] is None:
        load_crc_corr()

    return dashbio.VolcanoPlot(
        DATASET_CRC_EXPR_LEN_CORR['dataframe'],     
        effect_size_line = effects,
        **DATASET_CRC_EXPR_LEN_CORR['dataprops'],
    )

@callback(
        Output('volcanoplot-current-selection', 'children'),
        Input('dashbio-default-volcanoplot', 'clickData'),
    )
def render_volcano_selection(clickData):
    """Update rendering of data points upon changing x-value of vertical dashed lines."""

    if clickData is not None:
        annotation = clickData["points"][0]["text"]
        """
        Store rows in a dictionary by display value, as this is the only way I found to get the whole data from the clickData later on
        """
        if annotation and annotation.startswith("<br>GENE: "):
            annotation = annotation[10:]
            data_object = DATASET_CRC_EXPR_LEN_CORR['annotations'][annotation]  

            return html.Span([
                f"Selected GENE: ",  
                html.A(target='_blank', href=f"/search?query={data_object['name']}&genome=hg38",          children=f"{data_object['name']} ({data_object['ensembl_id']})"), 
                " and STR: ", 
                html.A(target='_blank', href=f"/locus?repeat_id={data_object['repeat_id']}&genome=hg38",  children=f"{data_object['chr']}_{data_object['start']}")
            ])
    return html.Span("Nothing selected")


def render_crc_expr_len_corr_volcano_plot():
    """
    This function returns a Div element that contains a range slider and a volcano plot.
    """
    return html.Div([
        'Effect sizes',
        dcc.RangeSlider(
            id='default-volcanoplot-input',
            min=-3,
            max=3,
            step=0.05,
            marks={i: {'label': str(i)} for i in range(-3, 3)},
            value=[-0.5, 0.5]
        ),
        html.Div(id='volcanoplot-current-selection'),
        html.Div(
            dcc.Graph(
                id='dashbio-default-volcanoplot'
            )
        ),
        html.Br(),
    ])

#################### Dash section used to render Volcano plot ###############
def add_dash_graphs_to_flask_server(server):    
    """
    Adds Dash endpoints to a Flask server.

    Parameters:
    server (Flask): The Flask server to add the Dash graphs to.

    Returns:
    None
    """

    try:
        load_crc_corr()
    except Exception as ex:
        logging.error(f"Failed to load CRC expression length correlation data, next attempt will be made on first request to /crc_research page. {ex}") 

    dash_app = dash.Dash(
        server=server, 
        routes_pathname_prefix="/dash/", 
        serve_locally=False
    )
    dash_app.layout = render_crc_expr_len_corr_volcano_plot()

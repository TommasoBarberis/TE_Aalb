"""
Author: Tommaso Barberis
Date: 27/04/2022
Description: script to inspect results from extension tests
Usage:
    On remote:
        conda activate TE_Aalb && python3 app.py -d test_dir/ -p 8000
    On local:
        ssh -L 8050:localhost:8000 user@hostname

    Then you can browse to localhost:8050 in your local machine to see the app
"""

from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import glob, re, shutil, os, flask, dash, argparse
from natsort import natsorted

server = flask.Flask('app')
server.secret_key = os.environ.get('secret_key', 'secret')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", metavar='', help=('''
    directory having results from test_extender.sh
    '''))
    parser.add_argument("-p", metavar='', default='8000', help=('''
    port where the server will be run
    '''))

    args = parser.parse_args()
    
    workdir = args.d        
    workdir = os.path.abspath(workdir)
    workdir = workdir + "/"
    print(workdir)


    port = int(args.p)

    currentdir = os.getcwd()
    seq_dirs = natsorted(glob.glob(workdir + "seq_*"))
    seq_dirs = [i.split("/")[-1] for i in seq_dirs]

    # Tabs
    tabs = []
    c = 0
    for seq in seq_dirs:
        val = "tab-" + str(c)
        tabs.append(dcc.Tab(label=seq, value=val, className='tab'))
        c += 1

    app = dash.Dash(__name__)
    app.title = '“Test for extension and edge polishing of consensi”'
    server = app.server

    app.layout = html.Div([
        dcc.Tabs(id='tabs',
            value='tab-0',
            vertical=True,
            children=tabs,
            persistence=True
        ),
        html.Div(id='tabs-content-example-graph')
    ])

    @app.callback(Output('tabs-content-example-graph', 'children'),
                Input('tabs', 'value'))
    def render_content(tab):
        tab_id = str(tab).split('-')[-1]
        prefix = workdir + "seq_" + tab_id + "/"
        
        title = "SEQUENCE " + tab_id
        
        # starting sequence
        with open(prefix + "seq.fasta", "r") as f:
            seq = f.readlines()[1]
        
            # copy number
        copy_number = 0
        with open(prefix + "seq.fasta.blast.bed.fa", "r") as f:
            copy_number = len(f.readlines())/2
            
            # divergence
        with open(prefix + "log.txt") as f:
            div_data = f.readlines()
        div_first = div_data[0].split("\t")[-1]

            # first TE-Aid
        origin_teaid = prefix + "seq.fasta.c2g.pdf"
        shutil.copyfile(origin_teaid, currentdir + "/assets/img/seq.fasta_" + tab_id + ".c2g.pdf")
        
        # second sequence
        with open(prefix + "no_clean_consensus.fasta", "r") as f:
            seq_no_clean = f.readlines()[1]

            # divergence
        div_no_clean = div_data[1].split("\t")[-1]

            # alignement
        no_clean_svg = prefix + "no_clean_output.svg"
        shutil.copyfile(no_clean_svg, currentdir + "/assets/img/no_clean_output_" + tab_id + ".svg")
            # TE-Aid
        no_clean_teaid = prefix + "no_clean_consensus.fasta.c2g.pdf"
        shutil.copyfile(no_clean_teaid, currentdir + "/assets/img/no_clean_consensus.fasta_" + tab_id + ".c2g.pdf")
        
        # third sequence
        with open(prefix + "no_ins_consensus.fasta", "r") as f:
            seq_no_ins = f.readlines()[1]
            
            #divergence
        div_no_ins = div_data[2].split("\t")[-1]

            # alignement
        no_ins_svg = prefix + "no_ins_output.svg"
        shutil.copyfile(no_ins_svg, currentdir + "/assets/img/no_ins_output_" + tab_id + ".svg")
            # TE-Aid
        no_ins_teaid = prefix + "no_ins_consensus.fasta.c2g.pdf"
        shutil.copyfile(no_ins_teaid, currentdir + "/assets/img/no_ins_consensus.fasta_" + tab_id + ".c2g.pdf")
        
        return html.Div([
            html.H1(title, id='tab_title'),
            html.H3("Starting consensus", className='old_cons'),
            html.Div([
                html.P("Sequence:", className='info'),
                html.P(seq, id='seq'),
                html.P("Length:", className='info'),
                html.P(len(seq)),
                html.P("Number of copies in the genome:", className='info'),
                html.P(copy_number),
                html.P("Divergence:", className='info'),
                html.P(div_first),
                html.P("TE-Aid results", className='info'),
                html.Iframe(src="/assets/img/seq.fasta_" + tab_id + ".c2g.pdf", className='pdf_viewer'),
            ]),
            html.H3("Alignement with flanking regions (1500pb)", className='old_cons'),
            html.Div([
                html.P("Sequence:", className='info'),
                html.P(seq_no_clean, id='seq'),
                html.P("Length:", className='info'),
                html.P(len(seq_no_clean)),
                html.P("Divergence:", className='info'),
                html.P(div_no_clean),
                html.P("Alignement:", className='info'),
                html.Img(src='assets/img/no_clean_output_' + tab_id + '.svg', className='ali'),
                html.P("TE-Aid results", className='info'),
                html.Iframe(src="/assets/img/no_clean_consensus.fasta_" + tab_id + ".c2g.pdf", className='pdf_viewer'),
            ]),
            html.H3("Alignement removing insertions", className='old_cons'),
            html.Div([
                html.P("Sequence:", className='info'),
                html.P(seq_no_ins, id='seq'),
                html.P("Length:", className='info'),
                html.P(len(seq_no_ins)),
                html.P("Divergence:", className='info'),
                html.P(div_no_ins),
                html.P("Alignement:", className='info'),
                html.Img(src='assets/img/no_ins_output_' + tab_id + '.svg', className='ali'),
                html.P("TE-Aid results", className='info'),
                html.Iframe(src="/assets/img/no_ins_consensus.fasta_" + tab_id + ".c2g.pdf", className='pdf_viewer'),
            ]),
        ])

    app.run_server(debug=True, port=port)
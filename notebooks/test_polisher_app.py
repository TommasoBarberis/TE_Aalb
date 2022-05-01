"""
Author: Tommaso Barberis
Date: #)/04/2022
Description: script to inspect results from polisher tests
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

server = flask.Flask('app')
server.secret_key = os.environ.get('secret_key', 'secret')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", metavar='', default='.', help=('''
    directory having results from the test of extender
    '''))
    parser.add_argument("-p", metavar='', default='8000', help=('''
    port where the server will be run on remote
    '''))

    args = parser.parse_args()
    
    workdir = args.d        
    workdir = os.path.abspath(workdir)
    workdir = workdir + "/"

    port = int(args.p)

    app = dash.Dash(__name__)
    app.title = '“Test for extension and edge polishing of consensi”'
    server = app.server

    drop_list = []
    for i in range(1,101):
        drop_list.append('seq ' + str(i))

    
    app.layout = html.Div(id='main-div', children=[
        html.Div(id='intro', children=[
            html.H1('Test for extension and edge polishing of consensi')
        ]), 
        dcc.Dropdown(drop_list, drop_list[0], id='seq-selecter'),
        html.Div(id='output_div', children=[])
    ])

    @app.callback(Output('output_div', 'children'), Input('seq-selecter', 'value'))
    def render_content(selected_seq):
        seq_id = selected_seq.split()[1]

        currentdir = os.getcwd()

        # it is necessary to copy the files outside the app
        #shutil.copyfile(currentdir + "/seq_" + seq_id + "/test.fasta.c2g.pdf", currentdir + "/assets/img/test.fasta_" + seq_id + ".c2g.pdf")
        #shutil.copyfile(currentdir + "/seq_" + seq_id + "/new_consensus.fasta.c2g.pdf", currentdir + "/assets/img/new_consensus.fasta_" + seq_id + ".c2g.pdf")

        return html.Div(id = 'result_div', children = [
            html.Div(id='pdf_viewer_left', children=[
                html.H2("Avant", className="h2-pdf-div"),
                html.Iframe(src="/assets/img/test.fasta_" + seq_id + ".c2g.pdf"),
            ]),
            html.Div(id='pdf_viewer_right', children=[
                html.H2("Après", className="h2-pdf-div"),
                html.Iframe(src="/assets/img/new_consensus.fasta_" + seq_id + ".c2g.pdf"),
            ])
        ])

    app.run_server(debug=True, port=port)
# import
# from jupyter_dash import JupyterDash
# import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import glob, re, shutil, os, flask, dash
from natsort import natsorted

server = flask.Flask('app')
server.secret_key = os.environ.get('secret_key', 'secret')

workdir="/home/lerat/TE_Aalb/test/extension/"
notebookdir = os.getcwd()
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
        value='tab-1',
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
    
    with open(prefix + "seq.fasta", "r") as f:
        seq = f.readlines()[1]
    
    copy_number = 0
    with open(prefix + "seq.fasta.blast.bed.fa", "r") as f:
        copy_number = len(f.readlines())/2
        
    origin_teaid = prefix + "seq.fasta.c2g.pdf"
    shutil.copyfile(origin_teaid, notebookdir + "/assets/img/seq.fasta_" + tab_id + ".c2g.pdf")
    
    with open(prefix + "no_clean_consensus.fasta", "r") as f:
        seq_no_clean = f.readlines()[1]
        
    no_clean_svg = prefix + "no_clean_output.svg"
    shutil.copyfile(no_clean_svg, notebookdir + "/assets/img/no_clean_output_" + tab_id + ".svg")
    
    no_clean_teaid = prefix + "no_clean_consensus.fasta.c2g.pdf"
    shutil.copyfile(no_clean_teaid, notebookdir + "/assets/img/no_clean_consensus.fasta_" + tab_id + ".c2g.pdf")
    
    with open(prefix + "no_ins_consensus.fasta", "r") as f:
        seq_no_ins = f.readlines()[1]
        
    no_ins_svg = prefix + "no_ins_output.svg"
    shutil.copyfile(no_ins_svg, notebookdir + "/assets/img/no_ins_output_" + tab_id + ".svg")
    
    no_ins_teaid = prefix + "no_ins_consensus.fasta.c2g.pdf"
    shutil.copyfile(no_ins_teaid, notebookdir + "/assets/img/no_ins_consensus.fasta_" + tab_id + ".c2g.pdf")
    
    return html.Div([
        html.H1(title, id='tab_title'),
        html.H3("Startinggg consensus", className='old_cons'),
        html.Div([
            html.P("Sequence:", className='info'),
            html.P(seq, id='seq'),
            html.P("Length:", className='info'),
            html.P(len(seq)),
            html.P("Number of copies in the genome:", className='info'),
            html.P(copy_number),
            html.P("TE-Aid results", className='info'),
            html.Iframe(src="/assets/img/seq.fasta_" + tab_id + ".c2g.pdf", className='pdf_viewer'),
        ]),
        html.H3("Alignement with flanking regions (1500pb)", className='old_cons'),
        html.Div([
            html.P("Sequence:", className='info'),
            html.P(seq_no_clean, id='seq'),
            html.P("Length:", className='info'),
            html.P(len(seq_no_clean)),
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
            html.P("Alignement:", className='info'),
            html.Img(src='assets/img/no_ins_output_' + tab_id + '.svg', className='ali'),
            html.P("TE-Aid results", className='info'),
            html.Iframe(src="/assets/img/no_ins_consensus.fasta_" + tab_id + ".c2g.pdf", className='pdf_viewer'),
        ]),
    ])

if __name__ == '__main__':
    app.run_server()
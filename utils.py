from Bio import SeqIO
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter
import pysam
import pandas as pd
import statistics

# parse data
def cons_parser(cons_file):
    seq_dict = {}
    for seq_record in SeqIO.parse(cons_file, "fasta"):
        name = seq_record.name
        try:
            name = name.split("#")[1].split("/")[0]
        except:
            name = "MITE"
        seq_len = len(seq_record.seq)

        if name not in seq_dict.keys():
            seq_dict[name] = [seq_len]
        else:
            tmp_list = seq_dict[name]
            tmp_list.append(seq_len)
            seq_dict[name] = tmp_list 
    
    df = pd.DataFrame(columns=["Annotation", "Length"])
    for k in seq_dict.keys():
        for val in seq_dict[k]:
            df2 = pd.DataFrame({'Annotation': [k], 'Length': [val]})

            df = pd.concat([df, df2], ignore_index = True, axis = 0)
            #df = df.append({"Annotation": k, "Length": val}, ignore_index=True)
    return df

# Generate histogram for length distribution
def length_hist(df, title):            
    fig = px.histogram(df, x="Length", color="Annotation", marginal="box", title=title)
    fig.show()

# Generate stats for length distribution
class Stats():
    def __init__(self, seq_info):
        self.seq_info = seq_info
        try:
            self.nb_consensi = seq_info.shape[0]
            self.longest = seq_info['Length'].max()
            self.shortest = seq_info['Length'].min()
            self.average = round(seq_info['Length'].mean(), 1)
            self.median = seq_info['Length'].median()
        except:
            self.nb_consensi = len(seq_info)
            self.longest = max(seq_info)
            self.shortest = min(seq_info)
            self.average = round(sum(seq_info)/len(seq_info), 1)
            self.median = statistics.median(seq_info)
                
    def _repr_html_(self):
        htmlString = f"""
        <h2>Summary</h2>
        <div style="width:50%; float:left">
            <table>
                <thead><tr><td></td><td></td></tr></thead>
                <tbody>
                    <tr><td><b>Number of consensi:</b></td><td>{self.nb_consensi}</td></tr>
                    <tr><td><b>Longest sequence (pb):</b></td><td>{self.longest}</td></tr>
                    <tr><td><b>Shortest sequence (pb):</b></td><td>{self.shortest}</td></tr>
                    <tr><td><b>Average length (pb):</b></td><td>{self.average}</td></tr>
                    <tr><td><b>Median length (pb):</b></td><td>{self.median}</td></tr>
                </tbody>
                <thead><tr><td></td><td></td></tr></thead>
            </table>
        </div>
        """

        htmlString += """
        <div>
            <table>
                <thead><tr><td></td><td></td></tr></thead>
                <tbody>
        """

        count = self.seq_info['Annotation'].value_counts()
        for te_type in count.keys():
            htmlString += f"""
                <tr><td><b>{te_type}</b></td><td>{count[te_type]}</td></tr>
            """
        htmlString += f"""
                <thead><tr><td></td><td></td></tr></thead>
                </tbody>
            </table>
        </div>
        """
        return htmlString

# Generate histogram for reads distribution
def reads_hist(samfile, consensi, title):
    sam = pysam.AlignmentFile(samfile, "r")
    df = pd.DataFrame(columns=['mapped reads', 'length'])
    
    mapped_cons = []
    for read in sam:
        if read.reference_name != None:
            mapped_cons.append(read.reference_name)
    
    mapped_cons = dict(Counter(mapped_cons))
    
    seq_dict = {}
    for seq_record in SeqIO.parse(consensi, "fasta"):
        name = seq_record.name
        seq_len = len(seq_record.seq)
        seq_dict[name] = seq_len
    
    for cons in mapped_cons.keys():
        mapped_cons[cons] = (mapped_cons[cons]*150)/seq_dict[cons]
    
    mapped_cons = list(mapped_cons.values())
    layout=go.Layout(
        title=title, 
        xaxis=dict(title="N.° mapped reads / consensus length"), 
        yaxis=dict(title="Number of consensi")
    )
    #fig = go.Figure(layout=layout)
    #fig.add_trace(go.Histogram(x=mapped_cons))
    df = pd.DataFrame(mapped_cons, columns =['Coverage of the consensus'])
    fig = px.histogram(df, x="Coverage of the consensus", marginal="box", title=title)
    #fig = px.scatter(df, x="mapped reads", y="length")
    fig.show()
    
# Parse data from RepBase
def parse_repbase(filename):
    LTR = ['Copia', 'Gypsy', 'BEL', 'LTR']
    LINE = ['R1', 'Jockey', 'I', 'CR1', 'L1', 'L2B', 'RTE', 'R4', 'Loa', 'R2', 'Hero', 'NeSL', 'L2', 'CRE', 'Tx1', \
            'Proto2', 'RandI', 'Proto1', 'L2A']
    SINE = ['SINE', 'SINE2/tRNA', 'SINE3/5S']
    non_LTR = ['Outcast', 'MINIME_DN', 'Non-LTR', 'PEN1', 'PEN4', 'PEN2', 'Daphne', 'RTEX', 'Penelope', 'Ingi', \
               'Kiri', 'Vingi', 'Crack', 'DIRS', 'Nimb', 'Rex1']
    DNA = ['Mariner/Tc1', 'P', 'DNA', 'hAT', 'Transib', 'piggyBac', 'EnSpm/CACTA', 'Helitron', 'Harbinger', \
         'IKIRARA1', 'Zator', 'VEGE_DW', 'P-element', 'MuDR', 'ISL2EU', 'Sola1', 'Ginger2/TDD', 'Academ', 'Merlin', \
            'PAT', 'Sola3', 'Ginger1', 'Sola2', 'Polinton', 'Kolobok', 'CryptonV',  'Crypton', 'CryptonI', 'IS3EU', \
         'Dada', 'CryptonA', 'MITE', 'STREPE_PF', 'Chapaev', 'Zisupton']
    UNKNOWN = ['Transposable', 'TELREP_AG', 'ALAD', 'DMFTZ', 'DMHMR2', 'DEC1_DS', 'DMHETRP', '5S_DM' ,'AY1', \
               'DMHMR1', 'DMFUSHI', 'ISFUN1', 'Nonautonomous', 'SCAR_MA', 'Repetitive', 'SNAPBACK_TC', 'LIRP1', \
               'OFU85403', 'RSAI', 'Interspersed', 'MBOI', 'LGRP1', 'ECORI_Hm', 'CCRP1', 'LDRP1', 'LARP1', 'OVRP1', \
              'LMRP1', 'Minicircle', 'GQRP1', 'GPRP1', 'R1A_SS', 'FR1', 'LARRP1', 'R1B_SS', 'STTREP_Mp', 'AT-rich', \
              'RP1_GL', 'EcoR1', 'APO1_AP', 'BR6_CP', 'HTE1', 'LDRP2', 'HHA1_BT', 'Origin', 'SAL_CL', 'PFRP3', \
               'STREPB_FA', 'APO2_AP', 'SCAI_EH', 'Multicopy', 'R1B_DS', 'REP-1_HM', 'AlKe1_AL', 'HHAI', 'PFRP5', \
               'REP-1_HMM', 'LARP2', 'RS3', 'BMRP1', 'ALBAMH1', 'AlKe6_AL', 'OARP1', 'EcoRI', 'transposon', 'tandem', \
              'SIRE', 'SCAR_MI', 'MINEX_Le', 'pSOS', 'TANDREP_TG', 'Ed_ERE1', 'Ribosomal', 'AFRP1', 'PENT_PU', \
               'EHINV1', 'DDTDD', 'REP-1_PXu', 'PFRP1', 'PFH76', 'SAU3A_TR' ]
    MSAT = ['MSAT']
    SAT = ['Satellite', 'ARS406', 'SAT', 'SZ23_TC']
    simple = ['Simple']
    pseudogene = ['rRNA', 'EHINV2', 'tRNA']
    seq_dict = {}
    for seq_record in SeqIO.parse(filename, "fasta"):
        try:
            name = seq_record.description.split()[1]
            if name in LTR:
                name = "LTR"
            if name in LINE:
                name = "LINE"
            if name in SINE:
                name = "SINE"
            if name in non_LTR:
                name = "non_LTR"
            if name in DNA:
                name = "DNA"
            if name in UNKNOWN:
                name = "UNKNOWN"
            if name in MSAT:
                name = "MSAT"
            if name in SAT:
                name = "SAT"        
            if name in simple:
                name = "simple"
            if name in pseudogene:
                name = "pseudogene"

            seq_len = len(seq_record.seq)

            if name not in seq_dict.keys():
                seq_dict[name] = [seq_len]
            else:
                tmp_list = seq_dict[name]
                tmp_list.append(seq_len)
                seq_dict[name] = tmp_list
        except:
            pass
    df = pd.DataFrame(columns=["Annotation", "Length"])
    for k in seq_dict.keys():
        for val in seq_dict[k]:
            df2 = pd.DataFrame({'Annotation': [k], 'Length': [val]})

            df = pd.concat([df, df2], ignore_index = True, axis = 0)
    return df
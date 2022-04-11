from re import S
from Bio import SeqIO
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter
import pysam
import pandas as pd
import statistics
from glob import glob
import re
import Levenshtein as lv

# parse data
def cons_parser(cons_file):
    df = pd.DataFrame(columns=["Seq name", "Annotation", "Length"])

    for seq_record in SeqIO.parse(cons_file, "fasta"):
        name = seq_record.name
        try:
            name = name.split("#")[1].split("/")[0]
        except:
            name = "MITE"
        seq_len = len(seq_record.seq)

        df2 = pd.DataFrame({'Seq name': [seq_record.name],'Annotation': [name], 'Length': [seq_len]})
        df = pd.concat([df, df2], ignore_index = True, axis = 0)            
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
    df = pd.DataFrame(columns=['Name', 'Coverage', 'length'])
    
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
        try:
            annot = name.split("#")[1].split("/")[0]
        except:
            annot = "MITE"
        seq_dict[name] = (annot, seq_len)
    
    for cons in mapped_cons.keys():
        nb_bases = mapped_cons[cons]*300/seq_dict[cons][1]
        
        df2 = pd.DataFrame({'Name': [cons], 'Annotation': [seq_dict[cons][0]],'Coverage': [nb_bases], 'length': [seq_dict[cons][1]]})
        df = pd.concat([df, df2], ignore_index = True, axis = 0)
    
    fig = px.scatter(df, x="length", y="Coverage", color="Annotation", title=title, hover_data=['Name'], marginal_x="histogram")
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

    df = pd.DataFrame(columns=['Name', 'Annotation', 'Length'])

    for seq_record in SeqIO.parse(filename, "fasta"):
        name = seq_record.id
        try:
            annot = seq_record.description.split()[1]
            if annot in LTR:
                annot = "LTR"
            if annot in LINE:
                annot = "LINE"
            if annot in SINE:
                annot = "SINE"
            if annot in non_LTR:
                annot = "non_LTR"
            if annot in DNA:
                annot = "DNA"
            if annot in UNKNOWN:
                annot = "UNKNOWN"
            if annot in MSAT:
                annot = "MSAT"
            if annot in SAT:
                annot = "SAT"        
            if annot in simple:
                annot = "simple"
            if annot in pseudogene:
                annot = "pseudogene"

            seq_len = len(seq_record.seq)
            
            df2 = pd.DataFrame({'Seq name': [name], 'Annotation': [annot], 'Length': [seq_len]})
            df = pd.concat([df, df2], ignore_index = True, axis = 0)
            
        except:
            pass
    return df

def reads_hist_repBase(samfile, consensi, title):
    sam = pysam.AlignmentFile(samfile, "r")
    df = pd.DataFrame(columns=['Name', 'Annotation', 'Coverage', 'length'])
    mapped_cons = []

    for read in sam:
        if read.reference_name != None:
            mapped_cons.append(read.reference_name)
    
    mapped_cons = dict(Counter(mapped_cons))

    seq_df = parse_repbase(consensi)

    for cons in mapped_cons.keys():
        nb_bases = mapped_cons[cons]*300/seq_df["Name"][cons]
        
        df2 = pd.DataFrame({'Name': [cons], 'Annotation': [seq_df['Name']['Annotation']],'Coverage': [nb_bases], 'length': [seq_dict[cons][1]]})
        df = pd.concat([df, df2], ignore_index = True, axis = 0)

    fig = px.scatter(df, x="length", y="Coverage", color="Annotation", title=title, hover_data=['Name'], marginal_x="histogram")
    fig.show()

    
def refiner_table(work_dir):
    """
    Parse directory with Refiner benchmark.
    """
    
    sampled_clst_dir = sorted(list(glob(work_dir+"/cluster_*", recursive = True)))
    clst_type = ["100", "150", "200", "250", "300", "350", "400", "450", "500"]
    info = ["nb_seq", "cons_length", "div"]

    columns = ["cluster"]
    for t in clst_type:
        for i in info:
            columns.append(i + "-" + t)

    refiner = pd.DataFrame(columns=columns)
    
    sizes = [100, 150, 200, 250, 300, 350, 400, 450, 500]

    for clst in sampled_clst_dir:
        row = []

        cluster_name = clst.split("/")[-1]
        row.append(cluster_name)

        cons_500 = clst + "/cluster-500.clst.fasta.refiner_cons"
        for seq_record in SeqIO.parse(cons_500, "fasta"):
            seq_500 = seq_record.seq

        for s in sizes:
            sample_name = clst + "/cluster-" + str(s) + ".clst.fasta"

            nb_seq = 0
            with open(sample_name, "r") as f:
                data = f.read()
                nb_seq = len(re.findall(r'>', data))
                row.append(nb_seq)

            cons_file = clst + "/cluster-" + str(s) + ".clst.fasta.refiner_cons"

            for seq_record in SeqIO.parse(cons_file, "fasta"):
                cons_len = len(seq_record.seq)
                row.append(cons_len)
                dist = lv.distance(str(seq_500), str(seq_record.seq))
                row.append(dist)
        
        refiner.loc[-1] = row
        refiner.index = refiner.index + 1
        refiner = refiner.sort_index()

    refiner = refiner.set_index('cluster')
    
    return refiner    
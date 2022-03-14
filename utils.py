from Bio import SeqIO
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter
import pysam

# parse data
def cons_parser(cons_file):
    seq_dict = {}
    for seq_record in SeqIO.parse(cons_file, "fasta"):
        name = seq_record.name
        name = name.split("#")[1].split("/")[0]
        seq_len = len(seq_record.seq)

        if name not in seq_dict.keys():
            seq_dict[name] = [seq_len]
        else:
            tmp_list = seq_dict[name]
            tmp_list.append(seq_len)
            seq_dict[name] = tmp_list 
    return seq_dict

# Generate histogram for length distribution
def length_hist(seq_dict):
    layout = go.Layout(
        title="Consensi length distribution",
        barmode='overlay',
        xaxis=dict(
            title='sequence length'
        ),
        yaxis=dict(
            title='count'
        )
    )
    
    fig = go.Figure(layout=layout)
    for k in seq_dict.keys():
        fig.add_trace(go.Histogram(x=seq_dict[k], name=k, nbinsx = 200))

    fig.update_layout(barmode='overlay')
    fig.show()

# Generate stats for length distribution
class Stats():
    def __init__(self, seq_info, classification):
        self.seq_info = seq_info
        self.classification = classification
        if classification:
            self.nb_consensi = sum(len(lst) for lst in seq_info.values())
            self.longest = max(max(lst) for lst in seq_info.values())
            self.shortest = min(min(lst) for lst in seq_info.values())
            self.average = round(sum(sum(lst) for lst in seq_info.values())/self.nb_consensi, 0)
        else:
            self.nb_consensi = len(seq_info)
            self.longest = max(seq_info)
            self.shortest = min(seq_info)
            self.average = round(sum(seq_info)/self.nb_consensi, 0)
                
    def _repr_html_(self):
        htmlString = f"""
        <h2>Summary</h2>
        <div style="width:50%; float:left">
            <table>
                <thead><tr><td></td><td></td></tr></thead>
                <tbody>
                    <tr><td><b>Number of consensi:</b></td><td>{self.nb_consensi}</td></tr>
                    <tr><td><b>Longest sequence:</b></td><td>{self.longest}</td></tr>
                    <tr><td><b>Shortest sequence:</b></td><td>{self.shortest}</td></tr>
                    <tr><td><b>Average length:</b></td><td>{self.average}</td></tr>                    
                </tbody>
                <thead><tr><td></td><td></td></tr></thead>
            </table>
        </div>
        """
        if self.classification:
            htmlString += """
            <div>
                <table>
                    <thead><tr><td></td><td></td></tr></thead>
                    <tbody>
            """
            for te_type in self.seq_info.keys():
                htmlString += f"""
                    <tr><td><b>{te_type}</b></td><td>{len(self.seq_info[te_type])}</td></tr>
                """
            htmlString += f"""
                    <thead><tr><td></td><td></td></tr></thead>
                    </tbody>
                </table>
            </div>
            """
        return htmlString

# Generate histogram for reads distribution
def reads_hist(samfile, consensi):
    sam = pysam.AlignmentFile(samfile, "r")
    
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
        mapped_cons[cons] = mapped_cons[cons]/seq_dict[cons]
    
    mapped_cons = list(mapped_cons.values())
    layout=go.Layout(title="Mapped reads distribution", xaxis=dict(title="N.° mapped reads / consensus length"), yaxis=dict(title="Number of consensi"))
    fig = go.Figure(layout=layout)
    fig.add_trace(go.Histogram(x=mapped_cons))
    fig.show()
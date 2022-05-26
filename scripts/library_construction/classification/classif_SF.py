import argparse


def classif(tsv_file):
    LTR = ['Copia', 'Gypsy', 'BEL', 'LTR', 'ATCOPIA32B_I', 'COPIA1', 'ERV', 'RLTR17D', 'GYPSY', 'AACOPIA1', 'Dp_DGLT', \
        'MER101E', 'ATCOPIA', 'GGERVL']
    LINE = ['R1', 'Jockey', 'I', 'CR1', 'L1', 'L2B', 'RTE', 'R4', 'Loa', 'R2', 'Hero', 'NeSL', 'L2', 'CRE', 'Tx1', \
            'Proto2', 'RandI', 'Proto1', 'L2A', 'JAM1', 'LOA', 'LINE', 'RT1', 'Odin', 'CER8', 'LIN4', 'TART', 'HAL']
    SINE = ['SINE', 'SINE2/tRNA', 'SINE3/5S', 'FEILAI', 'FEILAI_B', '5S', 'Gecko']
    non_LTR = ['Outcast', 'MINIME_DN', 'Non-LTR', 'PEN1', 'PEN4', 'PEN2', 'Daphne', 'RTEX', 'Penelope', 'Ingi', \
            'Kiri', 'Vingi', 'Crack', 'DIRS', 'Nimb', 'Rex1', 'Loner', 'Waldo', 'Lian', 'ORTE', 'HATN16_DR', 'Tad1', \
            'WUKONG', 'MosquI_Aa2', 'HACR1', 'Athena', 'Dong', 'SARTTc', 'HOPE', 'R6Ag', 'Neptune']
    DNA = ['Mariner/Tc1', 'P', 'DNA', 'hAT', 'Transib', 'piggyBac', 'EnSpm/CACTA', 'Helitron', 'Harbinger', \
        'IKIRARA1', 'Zator', 'VEGE_DW', 'P-element', 'MuDR', 'ISL2EU', 'Sola1', 'Ginger2/TDD', 'Academ', 'Merlin', \
            'PAT', 'Sola3', 'Ginger1', 'Sola2', 'Polinton', 'Kolobok', 'CryptonV',  'Crypton', 'CryptonI', 'IS3EU', \
        'Dada', 'CryptonA', 'STREPE_PF', 'Chapaev', 'Zisupton', 'Shinagawa', 'ENSPM', 'HATx1', 'EnSpm1', 'MuDrx', \
        'Tx_mos', 'Tc1', 'HELITRON', 'HAT', 'Mariner', 'MARINER', 'CHARLIE', 'Novosib', 'nPIF']
    MITE = ['MITE', 'otherMITEs', 'MIMO', 'WUJIN']
    UNKNOWN = ['Transposable', 'TELREP_AG', 'ALAD', 'DMFTZ', 'DMHMR2', 'DEC1_DS', 'DMHETRP', '5S_DM' ,'AY1', \
            'DMHMR1', 'DMFUSHI', 'ISFUN1', 'Nonautonomous', 'SCAR_MA', 'Repetitive', 'SNAPBACK_TC', 'LIRP1', \
            'OFU85403', 'RSAI', 'Interspersed', 'MBOI', 'LGRP1', 'ECORI_Hm', 'CCRP1', 'LDRP1', 'LARP1', 'OVRP1', \
            'LMRP1', 'Minicircle', 'GQRP1', 'GPRP1', 'R1A_SS', 'FR1', 'LARRP1', 'R1B_SS', 'STTREP_Mp', 'AT-rich', \
            'RP1_GL', 'EcoR1', 'APO1_AP', 'BR6_CP', 'HTE1', 'LDRP2', 'HHA1_BT', 'Origin', 'SAL_CL', 'PFRP3', \
            'STREPB_FA', 'APO2_AP', 'SCAI_EH', 'Multicopy', 'R1B_DS', 'REP-1_HM', 'AlKe1_AL', 'HHAI', 'PFRP5', \
            'REP-1_HMM', 'LARP2', 'RS3', 'BMRP1', 'ALBAMH1', 'AlKe6_AL', 'OARP1', 'EcoRI', 'transposon', 'tandem', \
            'SIRE', 'SCAR_MI', 'MINEX_Le', 'pSOS', 'TANDREP_TG', 'Ed_ERE1', 'Ribosomal', 'AFRP1', 'PENT_PU', \
            'EHINV1', 'DDTDD', 'REP-1_PXu', 'PFRP1', 'PFH76', 'SAU3A_TR', 'AMPLICON_AA', 'mTA_Ele42', 'ZEBEDEE', \
            'AeHerves', 'SLF14']
    SAT = ['Satellite', 'ARS406', 'SAT', 'SZ23_TC', 'TguSat1', 'MOSAT', 'HSATII']
    MSAT = ['MSAT', 'TREP60', 'MINISAT']
    SIMPLE = ['Simple', 'CACTA', 'MonoRep163']
    # PSEUDOGENE = ['rRNA', 'EHINV2', 'tRNA']
    rRNA = ['LSU', 'SSU']
    snRNA = ['U5B1']

    new_data = {}
    with open(tsv_file, 'r') as f:
        f.readline()
        data = f.readlines()
        for line in data:
            fields = line.split()
            q_name = fields[0]
            score = fields[1]
            div = fields[2]
            match_file = fields[3]
            classif = fields[4]
            
            try:
                match = match_file.split('-')
                match = match[0]
            except:
                match = match_file
                    
        
            for ltr in LTR:                    
                if match in ltr or ltr.startswith(match) or match.startswith(ltr):
                    classif = "LTR"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for line in LINE:
                if match in line or line.startswith(match) or match.startswith(line):
                    classif = "LINE"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for sine in SINE:
                if match in sine or sine.startswith(match) or match.startswith(sine):
                    classif = "SINE"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for non_ltr in non_LTR:
                if match in non_ltr or non_ltr.startswith(match) or match.startswith(non_ltr):
                    classif = "non_LTR"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for dna in DNA:
                if match == dna or dna.startswith(match) or match.startswith(dna):
                    classif = "DNA"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for mite in MITE:
                if match == mite or mite.startswith(match) or match.startswith(mite):
                    classif = "MITE"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for unknown in UNKNOWN:
                if match in unknown or unknown.startswith(match) or match.startswith(unknown):
                    classif = "UNKNOWN"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for msat in MSAT:
                if match in msat or msat.startswith(match) or match.startswith(msat):
                    classif = "MSAT"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for sat in SAT:
                if match in sat or sat.startswith(match) or match.startswith(sat):
                    classif = "SAT"        
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for simple in SIMPLE:
                if match in simple or simple.startswith(match) or match.startswith(simple):
                    classif = "simple"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            # for pseudogene in PSEUDOGENE:
            #     if match in pseudogene or pseudogene.startswith(match) or match.startswith(pseudogene):
            #         classif = "pseudogene"
            #         new_data[q_name] = (score, div, match_file, classif)
            #         break
            for rrna in rRNA:
                if match in rrna or rrna.startswith(match) or match.startswith(rrna):
                    classif = "rRNA"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            for snrna in snRNA:
                if match in snrna or snrna.startswith(match) or match.startswith(snrna):
                    classif = "snRNA"
                    new_data[q_name] = (score, div, match_file, classif)
                    break
            
    with open(tsv_file + '.reclassified', 'w') as f:
        for seq in new_data.keys():
            f.write(seq + '\t' + new_data[seq][0] + '\t' + new_data[seq][1] + '\t' + new_data[seq][2] + '\t' + new_data[seq][3] + '\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", metavar='', help=('''
    Filtered .out file at the .tsv format
    '''))

    args = parser.parse_args()

    # assign parameters
    tsv_file = args.f

    if tsv_file:
        classif(tsv_file)
    pass
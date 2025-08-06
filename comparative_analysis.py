#!/usr/bin/env python3
"""
Comparison Pipeline

Summarizes observed amino acid variant counts per sample and generates a count heatmap.

Features:
1. Aggregate observed amino acid variants into a single vector per sample.
2. Visualize count distributions across samples as a heatmap.

Usage:
    python comparison_and_modeling.py \
        --observed_dir <directory of observed AA matrix CSVs> \
        --comparison_dir <directory of comparison AA matrix CSVs> \
        --manifest <manifest file> \
        --out_dir <output directory> \
        [--vector_file <observed summary CSV>] \
        [--step summarize|heatmap|all]

Arguments:
    --observed_dir      Directory of observed AA CSVs named {sample-id}.csv
    --comparison_dir    Directory of comparison AA matrix CSVs named {sample-id}.csv
    --manifest          Tab-delimited manifest file listing sample IDs (first column)
    --out_dir           Directory to save outputs: summary and heatmap
    --vector_file       Optional path to save or load single summary CSV (default: out_dir/observed_summary.csv)
    --step              Step to run: summarize, compare, heatmap, single_file or all (default: all)

Dependencies:
    pandas
    numpy
    matplotlib
    seaborn
"""

import argparse
import logging
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

# Configure logging
logging.basicConfig(level=logging.INFO)

# Columns for summarizing observed variants
SUB_VECTOR = ['ID','SUM','AE','AG','AP','AS','AT','AV','CF','CG','CR','CS','CW','CX','CY',
              'DA','DE','DG','DH','DN','DV','DY','EA','ED','EG','EK','EQ','EV','EX','FC',
              'FI','FL','FS','FV','FY','GA','GC','GD','GE','GR','GS','GV','GW','GX','HD',
              'HL','HN','HP','HQ','HR','HY','IF','IK','IL','IM','IN','IR','IS','IT','IV',
              'IW','KE','KI','KM','KN','KQ','KR','KT','KX','LF','LH','LI','LM','LP','LQ',
              'LR','LS','LV','LW','LX','MI','MK','ML','MR','MT','MV','ND','NH','NI','NK',
              'NS','NT','NY','PA','PH','PL','PQ','PR','PS','PT','QE','QH','QK','QL','QP',
              'QR','QX','RC','RG','RH','RI','RK','RL','RM','RP','RQ','RS','RT','RW','RX',
              'SA','SC','SF','SG','SI','SL','SN','SP','SR','ST','SW','SX','SY','TA','TI',
              'TK','TM','TN','TP','TR','VA','VD','VE','VF','VG','VI','VL','VM','WC','WG',
              'WL','WR','WS','WX','XC','XE','XG','XK','XL','XQ','XR','XS','XW','XY','YC',
              'YD','YF','YH','YN','YS','YX']

# Columns for summarizing variants
HEADING = ['AE','AG','AP','AS','AT','AV','CF','CG','CR','CS','CW','CX','CY',
              'DA','DE','DG','DH','DN','DV','DY','EA','ED','EG','EK','EQ','EV','EX','FC',
              'FI','FL','FS','FV','FY','GA','GC','GD','GE','GR','GS','GV','GW','GX','HD',
              'HL','HN','HP','HQ','HR','HY','IF','IK','IL','IM','IN','IR','IS','IT','IV',
              'IW','KE','KI','KM','KN','KQ','KR','KT','KX','LF','LH','LI','LM','LP','LQ',
              'LR','LS','LV','LW','LX','MI','MK','ML','MR','MT','MV','ND','NH','NI','NK',
              'NS','NT','NY','PA','PH','PL','PQ','PR','PS','PT','QE','QH','QK','QL','QP',
              'QR','QX','RC','RG','RH','RI','RK','RL','RM','RP','RQ','RS','RT','RW','RX',
              'SA','SC','SF','SG','SI','SL','SN','SP','SR','ST','SW','SX','SY','TA','TI',
              'TK','TM','TN','TP','TR','VA','VD','VE','VF','VG','VI','VL','VM','WC','WG',
              'WL','WR','WS','WX','XC','XE','XG','XK','XL','XQ','XR','XS','XW','XY','YC',
              'YD','YF','YH','YN','YS','YX']

# Substitution matrix for expected AA vectors
SUB_MATRIX = [
    ['A',0,0,'AD','AE',0,'AG',0,0,0,0,0,0,'AP',0,0,'AS','AT','AV',0,0,0],
    ['C',0,0,0,0,'CF','CG',0,0,0,0,0,0,0,0,'CR','CS',0,0,'CW','CY','CX'],
    ['D','DA',0,0,'DE','DG','DH',0,0,0,0,'DN',0,0,0,0,0,'DV',0,'DY',0],
    ['E','EA',0,'ED',0,0,'EG',0,0,'EK',0,0,0,0,'EQ',0,0,0,'EV',0,0,'EX'],
    ['F',0,'FC',0,0,0,0,0,'FI',0,'FL',0,0,0,0,0,'FS',0,'FV',0,'FY',0],
    ['G','GA','GC','GD','GE',0,0,0,0,0,0,0,0,0,0,'GR','GS',0,'GV','GW',0,'GX'],
    ['H',0,0,'HD',0,0,0,0,0,'HL',0,'HN','HP','HQ','HR',0,0,0,0,'HY',0],
    ['I',0,0,0,0,'IF',0,0,0,'IK','IL','IM','IN',0,0,'IR','IS','IT','IV','IW',0,0],
    ['K',0,0,'KE',0,0,0,'KI',0,0,'KM','KN',0,'KQ','KR',0,'KT',0,0,0,'KX'],
    ['L',0,0,0,0,'LF',0,'LH','LI',0,0,'LM',0,'LP','LQ','LR','LS',0,'LV','LW',0,'LX'],
    ['M',0,0,0,0,0,0,0,'MI','MK','ML',0,0,0,0,'MR',0,'MT','MV',0,0,0],
    ['N',0,0,'ND',0,0,'NH','NI','NK',0,0,0,0,0,0,'NS','NT',0,0,'NY',0],
    ['P','PA',0,0,0,0,0,'PH',0,0,'PL',0,0,0,'PQ','PR','PS','PT',0,0,0,0],
    ['Q',0,0,0,'QE',0,0,'QH','QK','QL',0,0,'QP',0,'QR',0,0,0,0,0,'QX'],
    ['R',0,'RC',0,0,0,'RG','RH','RI','RK','RL','RM',0,'RP','RQ',0,'RS','RT',0,'RW',0,'RX'],
    ['S','SA','SC',0,0,'SF','SG',0,'SI',0,'SL',0,'SN','SP',0,'SR',0,'ST',0,'SW','SY','SX'],
    ['T','TA',0,0,0,0,0,0,'TI','TK',0,'TM','TN','TP',0,'TR',0,0,0,0,0,0],
    ['V','VA',0,'VD','VE','VF','VG',0,'VI',0,'VL','VM',0,0,0,0,0,0,0,0,0,0],
    ['W',0,'WC',0,0,0,'WG',0,0,0,'WL',0,0,0,0,'WR','WS',0,0,0,0,'WX'],
    ['Y',0,'YC','YD',0,'YF',0,'YH',0,0,0,0,'YN',0,0,0,'YS',0,0,0,0,'YX'],
    ['STOP',0,'XC',0,'XE',0,'XG',0,0,'XK','XL',0,0,0,'XQ','XR','XS',0,0,'XW','XY',0]
]


def parse_args():
    parser = argparse.ArgumentParser(description="Comparison Pipeline")
    parser.add_argument('--observed_dir', required=True,
                        help='Directory of observed AA matrix CSVs')
    parser.add_argument('--comparison_dir', required=True,
                        help='Directory of comparison AA matrix CSVs')
    parser.add_argument('--manifest', required=True,
                        help='Manifest file listing sample IDs')
    parser.add_argument('--out_dir', required=True,
                        help='Directory for outputs')
    parser.add_argument('--vector_file', required=False,
                        help='Path to summary CSV to skip compare step')
    parser.add_argument('--step', choices=['summarize','compare','heatmap','single-file','all'], default='all',
                        help='Step(s) to run')
    return parser.parse_args()




def load_observed(observed_dir, manifest):
    ids = pd.read_table(manifest, header=None).iloc[:,0].astype(str)
    obs = {}
    for sid in ids:
        fp = Path(observed_dir)/f"{sid}.csv"
        if fp.exists(): obs[sid] = pd.read_csv(fp, index_col=0)
        else: logging.warning(f"Missing observed for {sid}")
    return obs


def summarize_observed(observed, out_dir, vector_file=None):
    rows = []
    for sid, df in observed.items():
        row = {'ID': sid, 'SUM': len(df)}
        for code in SUB_VECTOR[2:]:
            r, c = code[0], code[1:]
            row[code] = ((df['ST']==r)&(df['END']==c)).sum()
        rows.append(row)
    summary = pd.DataFrame(rows).set_index('ID')
    summary = summary.reindex(columns=SUB_VECTOR[1:])
    out_path = Path(vector_file) if vector_file else Path(out_dir)/'observed_summary.csv'
    summary.to_csv(out_path)
    logging.info(f"Saved observed summary to {out_path}")
    return summary


def compare_vectors(observed_summary, out_dir, vector_file):
    sims = pd.DataFrame(index=observed_summary.index, columns=expected.index)
    for sid in observed_summary.index:
        obs_vec = observed_summary.loc[sid].values.astype(float)
        for sig in expected.index:
            exp_vec = expected.loc[sig].astype(float).values
            sims.loc[sid, sig] = np.dot(exp_vec, obs_vec)/(np.linalg.norm(exp_vec)*np.linalg.norm(obs_vec))
    out_path = Path(comparison_file) if comparison_file else Path(out_dir)/"similarity_matrix.csv"
    sims.to_csv(out_path)
    logging.info(f"Similarity matrix saved to {out_path}")
    return sims


def plot_heatmap(similarity, out_dir, vector_file=None):
    fig, ax = plt.subplots(figsize=(10,10))
    sb.heatmap(similarity.astype(float), cmap='viridis', ax=ax)
    hm_path = Path(vector_file) if vector_file else Path(out_dir)/"heatmap.png"
    plt.tight_layout()
    plt.savefig(hm_path)
    plt.close(fig)
    logging.info(f"Heatmap saved to {hm_path}")


def plot_cluster(observed_dir, out_dir, vector_file=None):
    ids = pd.read_table(manifest, header=None).iloc[:,0].astype(str)
    file1 = observed_dir + id + '.csv'
    M = pd.read_csv(file1)
    ID = M.iloc[:,0]
    M = M.iloc[:,1:]
    M = np.matrix(M, dtype=float)
    M = M/M.sum
    MOD = MOD.T
    MOD = MOD + 0.0000000000000000000001
    cbar_kws= {'label':'Percentage of Substitutions'}
    heat_map4 = sb.clustermap(MOD, vmin=0, vmax=10, metric='cosine', figsize=(30,40), cbar_kws=cbar_kws)
    plt.savefig(outpath)
    logging.info(f"Heatmap saved to {hm_path}")
   
  
def compare_multiple(observed_dir, comparison_dir, out_dir, vector_file=None):
    ids = pd.read_table(manifest, header=None).iloc[:,0].astype(str)
    file1 = observed_dir + id + '.csv'
    file2 = comparison_dir + id + '.csv'
    outpath = out_dir + id + '-compared.csv'
  
    M = pd.read_csv(file1)
    ID = M.iloc[:,0]
    M = M.iloc[:,1:]
    M = np.matrix(M, dtype=float)
    M = M/M.sum
    
    N = pd.read_csv(file2)
    N = N.iloc[:,1:]
    N = np.matrix(N, dtype=float)
    N = N/N.sum

    delt = M-N
    delt = delt + 0.0000000000000000000001
    delt = pd.DataFrame(delt)
    delt.to_csv(outpath, index=False, header=False)
    
    text1 = HEADING
    text2 = ID.T
    outpath = f"{Out}{id}-compared.png"
    delt = delt.T
    sb.set(rc={'figure.figsize':(30,40)})
    cbar_kws= {'shrink':0.25,'ticks':[-5,-4,-3,-2,-1,0,1,2,3,4,5],'label':'Percent of Substitutions','orientation':'vertical',"use_gridspec":False}
    heat_map = sb.heatmap(delt, vmin=-5, vmax=5,annot_kws = {'size':15}, cbar_kws=cbar_kws)
    heat_map.set_yticklabels(text1)
    heat_map.set_xticklabels(text2, rotation=90)
    plt.savefig(outpath)

  
def single_file_count(vector_file):
    df = pd.read_csv(vector_file)
    ID = df.iloc[:,0]
    text1 = HEADING
    text2 = ID.T
    out_path = out_dir + 'aa-count.png'
    MOD = df.iloc[:,2:]
    MOD = MOD.T
    MOD = MOD + 0.000000000000000000001
    sb.set(rc={'figure.figsize':(30,40)})
    cbar_kws= {'shrink':0.25,'ticks':[0,1,2,3,4,5,6,7,8,9,10],'label':'Number of Substitutions','orientation':'vertical',"use_gridspec":False}
    heat_map = sb.heatmap(MOD, vmin=0, vmax=10,annot_kws = {'size':15}, cbar_kws=cbar_kws)
    heat_map.set_yticklabels(text1)
    heat_map.set_xticklabels(text2, rotation=90)
    plt.savefig(outpath)


def single_file_proportion(vector_file):
    df = pd.read_csv(vector_file)
    ID = df.iloc[:,0]
    text1 = HEADING
    text2 = ID.T
    out_path = out_dir + 'aa-proportion.png'
    NSC = df.iloc[:,1]
    MOD = df.iloc[:,2:]
    MOD = MOD/NSC
    MOD = MOD.T
    MOD = MOD + 0.000000000000000000001
    sb.set(rc={'figure.figsize':(30,40)})
    cbar_kws= {'shrink':0.25,'ticks':[0,1,2,3,4,5],'label':'Percent of Substitutions','orientation':'vertical',"use_gridspec":False}
    heat_map = sb.heatmap(MOD, vmin=0, vmax=5,annot_kws = {'size':15}, cbar_kws=cbar_kws)
    heat_map.set_yticklabels(text1)
    heat_map.set_xticklabels(text2, rotation=90)
    plt.savefig(outpath)


def main():
    args = parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    signatures = load_signatures(args.signature_file)
    observed = load_observed(args.observed_dir,args.manifest)

    if args.step in ['summarize','all']:
        observed_summary = summarize_observed(observed,args.out_dir,args.vector_file)
    else:
        observed_summary = pd.read_csv(Path(args.vector_file) if args.vector_file else Path(args.out_dir)/'observed_summary.csv', index_col=0)

    if args.step in ['compare','all']:
        sims = compare_vectors(expected,observed_summary,args.out_dir,args.comparison_file)
    else:
        sims = pd.read_csv(Path(args.comparison_file) if args.comparison_file else Path(args.out_dir)/"similarity_matrix.csv", index_col=0)

    if args.step in ['heatmap','all']:
        plot_heatmap(sims,args.out_dir,args.comparison_file)

if __name__=='__main__':
    main()

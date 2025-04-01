import os
import sys
import gzip
import pandas as pd
import numpy as np
from Bio import SeqIO

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def ComputeRanges(distances, max_dist = 300000):
    ranges = []
    start_idx = -1
    end_idx = -1
    for idx, dist in enumerate(distances):
        if dist <= max_dist:
            if start_idx == -1:
                start_idx = idx
                end_idx = idx
            end_idx += 1
        else:
            if start_idx != -1:
                ranges.append((start_idx, end_idx))
            start_idx = -1
    if start_idx != -1:
        ranges.append((start_idx, end_idx))
    return ranges

def ComputeSummaryDF(dfs, contig_len_dict, shift):
    locus_df = {'LocusID' : [], 'Locus' : [], 'Contig' : [], 'StartPos' : [], 'EndPos' : [], 'Length' : [], 'GeneTypes' : [], 'NumV' : [], 'NumProdV' : [], 'FracProdV' : [], 'RelStart' : [], 'RelEnd' : []}
    locus_id = 1
    for df in dfs:
        if len(df) == 0:
            continue
        locus = df['Locus'][0]
        contigs = set(df['Contig'])
        for contig in contigs:
            contig_df = df.loc[df['Contig'] == contig]
            if len(contig_df) <= 2:
                continue
            positions = sorted(contig_df['Pos'])
            distances = []
            for i in range(1, len(positions)):
                distances.append(positions[i] - positions[i - 1])
            ranges = ComputeRanges(distances)
            for start_idx, end_idx in ranges:
                start_pos = positions[start_idx]
                end_pos = positions[end_idx]
                range_df = contig_df.loc[(contig_df['Pos'] >= start_pos) & (contig_df['Pos'] <= end_pos)]
                range_v_df = range_df.loc[range_df['GeneType'] == 'V']
                range_prod_v_df = range_v_df.loc[range_v_df['Productive']]
                locus_df['LocusID'].append(locus_id)
                locus_df['Locus'].append(locus)
                locus_df['Contig'].append(contig)
                refined_start = max(0, start_pos - shift)
                locus_df['StartPos'].append(refined_start)
                refined_end = min(end_pos + shift, contig_len_dict[contig])
                locus_df['EndPos'].append(refined_end)
                locus_df['Length'].append(refined_end - refined_start)
                locus_df['GeneTypes'].append(','.join(sorted(set(range_df['GeneType']))))
                locus_df['NumV'].append(len(range_v_df))
                locus_df['NumProdV'].append(len(range_prod_v_df))
                frac_prod = 0
                if len(range_v_df) != 0:
                    frac_prod = len(range_prod_v_df) / len(range_v_df)
                locus_df['FracProdV'].append(frac_prod)
                locus_df['RelStart'].append(refined_start / contig_len_dict[contig])
                locus_df['RelEnd'].append(refined_end / contig_len_dict[contig])
                locus_id += 1
    locus_df = pd.DataFrame(locus_df)
    return locus_df

def VisualizeSummary(locus_df, output_fname):
    if len(locus_df) == 0:
        return
    locus_colors = {'IGH' : '#9367BD', 'IGK' : 'orange', 'IGL' : '#2AA02B', 'TRA' : '#D62727', 'TRB' : '#1F77B4', 'TRG' : '#E377C1'}
    fig, axes = plt.subplots(nrows = 3, figsize = (15, 10))
    x = np.array(range(len(locus_df)))
    colors = [locus_colors[locus_df['Locus'][i]] for i in range(len(locus_df))]
    # lengths
    sns.barplot(x = x, y = locus_df['Length'], palette = colors, ax = axes[0])
    plt.sca(axes[0])
    plt.xticks([], [])
    # fraction productive Vs
    sns.barplot(x = x, y = locus_df['NumProdV'], palette = colors, ax = axes[1])
    plt.sca(axes[1])
    plt.xticks([], [])
    # positions
    rel_pos_df = {'Index' : [], 'RelPos' : []}
    for i in range(len(locus_df)):
        rel_pos_df['Index'].append(i)
        rel_pos_df['RelPos'].append(locus_df['RelStart'][i])
        rel_pos_df['Index'].append(i)
        rel_pos_df['RelPos'].append(locus_df['RelEnd'][i])
    rel_pos_df = pd.DataFrame(rel_pos_df)
    sns.swarmplot(x = 'Index', y = 'RelPos', data = rel_pos_df, palette = colors, ax = axes[2])
    plt.sca(axes[2])
    plt.xticks(x, locus_df['Contig'], rotation = 90, fontsize = 6)
    plt.ylim(0, 1.01)
    plt.xlabel('')
    # output
    plt.tight_layout()
    plt.savefig(output_fname, dpi = 300)
    plt.clf()

def GetRangeBasename(summary_df, idx):
    return summary_df['Locus'][idx] + '_' + summary_df['Contig'][idx] + '_' + str(summary_df['NumV'][idx]) + 'Vs'

def OutputLociToFastaFiles(summary_df, contig_seq_dict, output_dir):
    for i in range(len(summary_df)):
        contig_seq = contig_seq_dict[summary_df['Contig'][i]]
        contig_len = len(contig_seq)
        start_pos = summary_df['StartPos'][i]
        end_pos = summary_df['EndPos'][i]
        fragment = contig_seq[start_pos : end_pos]
        fname = os.path.join(output_dir, GetRangeBasename(summary_df, i) + '.fasta')
        fh = open(fname, 'w')
        fh.write('>' + summary_df['Contig'][i] + '_' + str(summary_df['LocusID'][i]) + '_' + summary_df['Locus'][i] + '\n' + fragment + '\n')
        fh.close()

def VisualizeGenePositions(gene_df, locus_df, output_dir):
    gene_color = {'V' : '#1F77B4', 'D' : '#FF7F0F', 'J' : '#2AA02B'}
    for i in range(len(locus_df)):
        contig = locus_df['Contig'][i]
        start_pos = locus_df['StartPos'][i]
        end_pos = locus_df['EndPos'][i]
        sub_df = gene_df.loc[(gene_df['Contig'] == contig) & (gene_df['Pos'] >= start_pos) & (gene_df['Pos'] <= end_pos)]
        if len(sub_df) == 0:
            continue
        sub_df = sub_df.sort_values(by = 'Pos').reset_index()
        scale = 300
        scaled_pos = []
        color_list = []
        for j in range(len(sub_df)):
            scaled_pos.append((sub_df['Pos'][j] - start_pos) / (end_pos - start_pos) * scale)
            color_list.append(gene_color[sub_df['GeneType'][j]])
        plt.figure(figsize = (6, 3))
        plt.bar(scaled_pos, [1] * len(scaled_pos), color = color_list)
        plt.xlim(0, scale)
        plt.xticks([], [])
        plt.yticks([], [])
        plt.savefig(os.path.join(output_dir, GetRangeBasename(locus_df, i) + '.png'), dpi = 300)
        plt.clf()
        plt.close()

def main(genome_fasta, input_dir, output_dir):
    files = ['combined_genes_IGH.txt', 'combined_genes_IGK.txt', 'combined_genes_IGL.txt', 'combined_genes_TRA.txt', 'combined_genes_TRB.txt', 'combined_genes_TRG.txt']
    dfs = [pd.read_csv(os.path.join(input_dir, fname), sep = '\t', dtype = {'Contig' : str}) for fname in files]
    df = pd.concat(dfs)
    #### reading contigs and contig lengths
    contig_set = set(df['Contig'])
    contig_len_dict = dict()
    contig_seq_dict = dict()
    if genome_fasta.endswith('.gz'):
        genome_fasta_handle = gzip.open(genome_fasta, 'rt')
    else:
        genome_fasta_handle = open(genome_fasta, 'r')
    for r in SeqIO.parse(genome_fasta_handle, 'fasta'):
        contig_id = r.id.replace('|', '_')
        if contig_id not in contig_set:
            continue
        contig_len_dict[contig_id] = len(r.seq)
        contig_seq_dict[contig_id] = str(r.seq)
    genome_fasta_handle.close()

    #### creating summary dataframe
    shift = 10000
    locus_df = ComputeSummaryDF(dfs, contig_len_dict, shift)
    locus_df = locus_df.sort_values(by = ['Contig', 'RelStart'], ascending = True).reset_index()
    locus_df.to_csv(os.path.join(output_dir, 'summary.csv'), index = False, columns = ['LocusID', 'Locus', 'Contig', 'StartPos', 'EndPos', 'Length', 'GeneTypes', 'NumV', 'NumProdV', 'FracProdV', 'RelStart', 'RelEnd'])
    VisualizeSummary(locus_df, os.path.join(output_dir, 'summary.png'))
    
    #### output IG loci into fasta
    igloci_fasta_dir = os.path.join(output_dir, 'igloci_fasta')
    if not os.path.exists(igloci_fasta_dir):
        os.mkdir(igloci_fasta_dir)
    OutputLociToFastaFiles(locus_df, contig_seq_dict, igloci_fasta_dir)

    #### visualize positions of genes with ranges
    gene_pos_dir = os.path.join(output_dir, 'gene_pos_plots')
    if not os.path.exists(gene_pos_dir):
        os.mkdir(gene_pos_dir)
    VisualizeGenePositions(df, locus_df, gene_pos_dir)

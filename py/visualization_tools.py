import os
import pandas as pd
import numpy as np

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def OutputHeatmap(filenames, output_fname):
    dfs = [pd.read_csv(fname, sep = '\t') for fname in filenames if os.path.exists(fname)]
    genes = ['IGHV', 'IGHD', 'IGHJ', 'IGKV', 'IGKJ', 'IGLV', 'IGLJ']
    contigs = []
    for df in dfs:
        contig_set = set(df['Contig'])
        for contig in contig_set:
            if contig not in contigs:
                contigs.append(contig)
    matrix = []
    annot_matrix = []
    for contig in contigs:
        matrix.append([0] * len(genes))
        annot_matrix.append([''] * len(genes))
    df = pd.concat(dfs).reset_index()
    df['LongGeneType'] = [df['Locus'][i] + df['GeneType'][i] for i in range(len(df))]
    for x_idx, contig in enumerate(contigs):
        for y_idx, gene in enumerate(genes):
            sub_df = df.loc[(df['LongGeneType'] == gene) & (df['Contig'] == contig)]
            if len(sub_df) != 0:
                matrix[x_idx][y_idx] = len(sub_df)
                annot_matrix[x_idx][y_idx] = str(len(sub_df))
    plt.figure(figsize = (12, 8))
    sns.heatmap(matrix, annot = np.array(annot_matrix), yticklabels = contigs, xticklabels = genes, cmap = 'coolwarm', robust = True, fmt = '', cbar = False)
    plt.yticks(fontsize = 6)
    plt.savefig(output_fname, dpi = 300)
    plt.clf()

def OutputPositionsPerContig(filename, locus, output_dir):
    if not os.path.exists(filename):
        return
    df = pd.read_csv(filename, sep = '\t')
    contigs = set(df['Contig'])
    color_dict = {('V', True) : '#1F77B4', ('V', False) : '#AEC7E8', 'D' : '#FF7F0F', 'J' : '#2AA02B'}
    strand_dict = {'+' : 1, '-' : -1}
    scale = 500
    for contig in contigs:
        contig_df = df.loc[df['Contig'] == contig].reset_index()
        min_pos = max(0, contig_df['Pos'][0] - 10000)
        max_pos = contig_df['Pos'][len(contig_df) - 1] + 10000
        scaled_pos_list = []
        color_list = []
        y_list = []
        for i in range(len(contig_df)):
            scaled_pos_list.append((contig_df['Pos'][i] - min_pos) / (max_pos - min_pos) * scale)
            gene_key = contig_df['GeneType'][i]
            if gene_key == 'V':
                gene_key = (gene_key, contig_df['Productive'][i])
            color_list.append(color_dict[gene_key])
            y_list.append(1) # * strand_dict[contig_df['Strand'][i]])
        plt.figure(figsize = (6, 3))
        plt.bar(scaled_pos_list, y_list, color = color_list)
        plt.ylim(0, 1)
        plt.xlim(0, scale)
        plt.xticks([], [])
        plt.yticks([], [])
        plt.savefig(os.path.join(output_dir, locus + '_' + contig + '.png'), dpi = 300)
        plt.clf()

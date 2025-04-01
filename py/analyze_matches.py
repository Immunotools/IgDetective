import os
import sys
import gzip
import shutil
import re
import numpy as np

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO

class Match:
    def __init__(self, gene_id, cigar):
        self.gene_id = gene_id
        self.cigar = cigar
        self.cigar_splits, self.cigar_delims = GetCigarSplits(cigar)

    def GeneId(self):
        return self.gene_id

    def MatchLength(self):
        match_length = 0
        for s, d in zip(self.cigar_splits, self.cigar_delims):
            if d == 'M':
                match_length += s
        return match_length

    def __repr__(self):
        return self.gene_id

def GetCigarSplits(cigar):
    delims = ['M', 'S', 'H', 'I', 'D']
    delim_str = '|'.join(delims)
    splits = re.split(delim_str, cigar)
    if splits[-1] == '':
        splits = splits[ : len(splits) - 1]
    cigar_pos = 0
    delims = []
    for s in splits:
        delim = cigar[cigar_pos + len(s) : cigar_pos + len(s) + 1]
        delims.append(delim)
        cigar_pos += len(s) + 1
    splits = [int(s) for s in splits]
    return splits, delims

def AnalyzeMatches(sam_file):
    contig_match_dict = dict()
    lines = open(sam_file).readlines()
    for l in lines:
        if l[0] == '@':
            continue
        line_splits = l.strip().split()
        read_id = line_splits[0]
        ref_id = line_splits[2]
        cigar = line_splits[5]
        if cigar == '*':
            continue
        start_pos = int(line_splits[3])
        if ref_id not in contig_match_dict:
            contig_match_dict[ref_id] = dict()
        if start_pos not in contig_match_dict[ref_id]:
            contig_match_dict[ref_id][start_pos] = []
        contig_match_dict[ref_id][start_pos].append(Match(read_id, cigar))
    return contig_match_dict

def CompressMatches(match_dict, gene_type):
    # match dict: start pos -> [Match]
    window_dict = {'V' : 200, 'J' : 30, 'C' : 500}
    window_size = window_dict[gene_type.split('-')[0]]
    cur_match_index = 0
    compressed_positions = [] # positions
    best_hits = [] # Matches with longest lengths
    processed_matches = set()
    sorted_start_positions = sorted(match_dict.keys())
    for pos in sorted_start_positions: #range(sorted_start_positions[0], sorted_start_positions[-1] + 1):
        if pos in processed_matches:
            continue
        bounds = (pos, pos + window_size)
        num_bound_matches = 0 # the number of matches within bounds
        bound_matches = [] # matches within bounds
        for match_ind in range(cur_match_index, len(sorted_start_positions)):
            match_pos = sorted_start_positions[match_ind]
            if match_pos >= bounds[0] and match_pos < bounds[1]:
                num_bound_matches += 1
                processed_matches.add(match_pos)
                for m in match_dict[match_pos]:
                    bound_matches.append(m)
            else:
                cur_match_index = match_ind
                break
        if num_bound_matches > 0:
#            compressed_positions.append((pos, num_bound_matches))
            best_match = ''
            best_length = 0
            for m in bound_matches:
                if m.MatchLength() > best_length:
                    best_length = m.MatchLength()
                    best_match = m
            best_hits.append(best_match)
            compressed_positions.append((pos, best_match))
    return compressed_positions

def OutputLoci(contig_id, contig_seq, loci_bounds, output_dir):
    contig_id = contig_id.replace('|', '_')
    shift_length = 1000
    total_bounds = (len(contig_seq), 0)
    for locus, gene in loci_bounds:
        gene_locus_bounds = loci_bounds[(locus, gene)]
        real_bounds = (max(0, gene_locus_bounds[0] - 1000), min(len(contig_seq), gene_locus_bounds[1] + 1000))
        contig_subseq = contig_seq[real_bounds[0] : real_bounds[1]]
        output_fname = os.path.join(output_dir, contig_id + '_' + locus + gene + '.fasta')
        output_fh = open(output_fname, 'w')
        output_fh.write('>CONTIG:' + contig_id + '|START_POS:' + str(real_bounds[0]) + '|END_POS:' + str(real_bounds[1]) + '|LOCUS:' + locus + '|GENE:' + gene + '\n')
        output_fh.write(contig_subseq + '\n')
        total_bounds = (min(total_bounds[0], real_bounds[0]), max(total_bounds[1], real_bounds[1]))
    # output the whole contig
    locus = ','.join(sorted(set([l[0] for l in loci_bounds])))
    genes = ','.join(sorted(set([l[1] for l in loci_bounds])))
    fh = open(os.path.join(output_dir, locus + '_' + contig_id + '.fasta'), 'w')
    fh.write('>CONTIG:' + contig_id + '|GENES:' + genes + '\n')
    fh.write(contig_seq + '\n')
    fh.close()
    # output overall locus
    fh = open(os.path.join(output_dir, contig_id + '_IG.fasta'), 'w')
    fh.write('>CONTIG:' + contig_id + '|START_POS:' + str(total_bounds[0]) + '|END_POS:' + str(total_bounds[1]) + '\n')
    fh.write(contig_seq[total_bounds[0] : total_bounds[1]] + '\n')
    fh.close()
    # output IGHD locus if possible
    if ('IGH', 'V') in loci_bounds and ('IGH', 'J') in loci_bounds:
        output_fname = os.path.join(output_dir, contig_id + '_IGHD.fasta')
        ighd_bounds = (len(contig_seq), 0)
        reverse = 'False'
        if loci_bounds[('IGH', 'V')][1] < loci_bounds[('IGH', 'J')][0]:
            ighd_bounds = (loci_bounds[('IGH', 'V')][1], loci_bounds[('IGH', 'J')][0])
        elif loci_bounds[('IGH', 'V')][0] > loci_bounds[('IGH', 'J')][1]:
            ighd_bounds = (loci_bounds[('IGH', 'J')][1], loci_bounds[('IGH', 'V')][0])
            reverse = 'True'
        if ighd_bounds[0] > ighd_bounds[1]:
            print('ERROR: IGHD locus is impossible to localize', contig_id, loci_bounds[('IGH', 'V')], loci_bounds[('IGH', 'J')])
            sys.exit(1)
        output_fh = open(output_fname, 'w')
        output_fh.write('>CONTIG:' + contig_id + '|START_POS:' + str(ighd_bounds[0]) + '|END_POS:' + str(ighd_bounds[1]) + '|LOCUS:IGH|GENE:D|REVERSE:' + reverse + '\n')
        output_fh.write(contig_seq[ighd_bounds[0] : ighd_bounds[1]] + '\n')
        output_fh.close()

######################################################
input_dir = sys.argv[1]
output_dir = sys.argv[2]
contig_file = sys.argv[3]

if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.mkdir(output_dir)

loci = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG']
genes = ['V', 'J', 'C']

input_files = os.listdir(input_dir)

gene_sam_dict = dict()
for f in input_files:
    for l in loci:
        for g in genes:
            gene = l + g
            if f.find(gene) != -1 and f.find('sam') != -1:
                gene_sam_dict[(l, g)] = os.path.join(input_dir, f)

contig_matches = dict() # exact positions of matches
combined_matches = dict() # contig -> locus, gene -> compressed positions
for locus, gene in gene_sam_dict:
    matches = AnalyzeMatches(gene_sam_dict[(locus, gene)]) # config -> start pos -> list of alignments
    for contig in matches:
        if contig not in contig_matches:
            contig_matches[contig] = dict()
            combined_matches[contig] = dict()
        contig_matches[contig][(locus, gene)] = matches[contig]
        compressed_matches = CompressMatches(matches[contig], gene)
        print(locus, gene, contig, compressed_matches)
        combined_matches[contig][(locus, gene)] = compressed_matches

contig_seqs = dict() # contig ID -> seq
if contig_file.endswith('.gz'):
    contig_handle = gzip.open(contig_file, 'rt')
else:
    contig_handle = open(contig_file, 'r')
for r in SeqIO.parse(contig_handle, 'fasta'):
    if r.id not in combined_matches:
        continue
    contig_seqs[r.id] = str(r.seq)
contig_handle.close()

matrix = []
annot_matrix = []
ylabels = []
for contig in contig_matches:
    num_matches = []
    for l in loci:
        for g in genes:
            if (l, g) in combined_matches[contig]:
                num_matches.append(len(combined_matches[contig][(l, g)]))
            else:
                num_matches.append(0)
    matrix.append(num_matches)
    annot_row = [''] * len(num_matches)
    for i in range(len(num_matches)):
        if num_matches[i] != 0:
            annot_row[i] = str(num_matches[i])
    annot_matrix.append(annot_row)
    ylabels.append(contig)

xlabels = []
for l in loci:
    for g in genes:
        xlabels.append(l + g)   

if len(matrix) != 0:
    plt.figure(figsize = (12, 8))
    plt.title(input_dir) 
    sns.heatmap(matrix, annot = np.array(annot_matrix), cmap = 'coolwarm', robust = True, fmt = '', xticklabels = xlabels, yticklabels = ylabels, cbar = False)
    plt.yticks(fontsize = 6)
    plt.savefig(os.path.join(output_dir, '__summary.png'), dpi = 300)
    plt.clf()

summary_txt = os.path.join(output_dir, '__summary.txt')
summary_fh = open(summary_txt, 'w')
summary_fh.write('ContigID\tContigLength\tLocus\tGeneType\tPosition\tGeneName\n')
for contig in contig_matches:
    contig_str = contig.replace('|', '_')
#    print(contig)
    locus_gene_matches = combined_matches[contig]
    max_pos = 0
    min_pos = sys.maxsize
    for locus, gene in locus_gene_matches:
        matches_pos = [p[0] for p in locus_gene_matches[(locus, gene)]]
        max_pos = max(max_pos, max(matches_pos))
        min_pos = min(min_pos, min(matches_pos))
    gene_bounds = dict() # geneType -> bounds
    for locus, gene in locus_gene_matches:
        locus_gene_bounds = (sys.maxsize, 0)
        # txt writing
        for pos, gene_name in sorted(locus_gene_matches[(locus, gene)], key = lambda x : x[0]):
            summary_fh.write(contig_str + '\t' + str('-') + '\t' + locus + '\t' + gene + '\t' + str(pos) + '\t' + str(gene_name) + '\n')
            locus_gene_bounds = (min(locus_gene_bounds[0], pos), max(locus_gene_bounds[1], pos))
        # gene bounds
        gene_bounds[(locus, gene)] = locus_gene_bounds
    # extracting loci and subloci
    OutputLoci(contig_str, contig_seqs[contig], gene_bounds, output_dir)

summary_fh.close()

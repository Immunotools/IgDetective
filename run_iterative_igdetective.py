import os
import sys
import subprocess
import shutil
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

sys.path.append('py')
import extract_aligned_genes as extr_genes

ref_gene_dir = 'datafiles/human_reference_genes'

def CheckPythonVersionFatal():
    if sys.version_info.major != 3:
        print('ERROR: Python 3 was not found. Please install Python 3 and rerun IgDetective')
        sys.exit(1)

def CheckMinimapFatal():
    try:
        subprocess.call(['which', 'minimap2'])
    except FileNotFoundError:
        print("ERROR: minimap2 was not found. Please install minimap2 and rerun IgDetective")
        sys.exit(1)

def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        print('WARN: output directory ' + output_dir + ' exists and will be overwritten!')
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

def GetRange(min_pos, max_pos, seq_len, max_len = 10000000):
    prefix_len = max_pos
    suffix_len = seq_len - min_pos
    gap = 1000000
    if prefix_len > max_len and suffix_len > max_len:
        return (max(min_pos - gap, 0), min(prefix_len + gap, seq_len))
    if prefix_len < suffix_len:
        return (0, min(prefix_len + gap, seq_len))
    return (max(min_pos - gap, 0), seq_len)

def AlignIgGenes(genome_fasta, ig_gene_fasta, sam_file):
    os.system('minimap2 -a ' + genome_fasta + ' ' + ig_gene_fasta + ' -o ' + sam_file + ' > /dev/null 2>&1')

def AlignReferenceGenes(align_dir, genome_fasta, output_dir):
    ref_gene_dict = dict()
    for f in os.listdir(ref_gene_dir):
        gene_type = f.split('.')[0]
        if gene_type not in ref_gene_dict:
            ref_gene_dict[gene_type] = os.path.join(ref_gene_dir, f)
    for gene_type in ref_gene_dict:
        print('Aligning ' + gene_type + ' genes (' + ref_gene_dict[gene_type] + ')...')
        AlignIgGenes(genome_fasta, ref_gene_dict[gene_type], os.path.join(align_dir, gene_type + '.sam'))

def IdentifyIGContigs(igcontig_dir, alignment_dir, output_dir, genome_fasta):
    match_log = igcontig_dir + '.out'
    os.system('python py/analyze_matches.py ' + alignment_dir + ' ' + igcontig_dir + ' ' + genome_fasta + ' > ' + match_log)

def GetPositionRange(sorted_positions):
    if len(sorted_positions) == 1:
        return sorted_positions[0], sorted_positions[0]
    distances = [sorted_positions[i] - sorted_positions[i - 1] for i in range(1, len(sorted_positions))]
    median = np.median(distances)
    start_pos = sorted_positions[0]
    end_pos = sorted_positions[-1]
    if distances[0] / median > 100:
        start_pos = sorted_positions[1]
    if distances[-1] / median > 100:
        end_pos = sorted_positions[-2]
    return start_pos, end_pos

def RunIgDetective(igcontig_dir, output_dir, locus = 'IGH'):
    print('==== Running RSS-based IgDetective for ' + locus + '...')
    txt = os.path.join(igcontig_dir, '__summary.txt')
    if not os.path.exists(txt):
        return
    df = pd.read_csv(txt, sep = '\t')
    igh_df = df.loc[df['Locus'] == locus]
    if len(igh_df) == 0:
        return
    fasta = os.path.join(output_dir, 'combined_contigs_' + locus + '.fasta')
    fh = open(fasta, 'w')
    contigs = set(igh_df['ContigID'])
    for c in contigs:
        c_df = igh_df.loc[igh_df['ContigID'] == c]
        contig_fasta = os.path.join(igcontig_dir, locus + '_' + str(c) + '.fasta')
        if not os.path.exists(contig_fasta):
            print('WARN: ' + contig_fasta + ' does not exist')
            continue
        seq = ''
        seq_id = ''
        for r in SeqIO.parse(contig_fasta, 'fasta'):
            seq = str(r.seq)
            seq_id = r.id
        #min_pos = min(c_df['Position'])
        #max_pos = max(c_df['Position'])
        positions = sorted(c_df['Position'])
        min_pos, max_pos = GetPositionRange(positions)
        contig_range = GetRange(min_pos, max_pos, len(seq))
        print('Contig: ' + c + ', contig range: ' + str(contig_range) + ', approx locus length: ' + str(contig_range[1] - contig_range[0]))
        fh.write('>' + seq_id + '|START:' + str(contig_range[0]) + '|END:' + str(contig_range[1]) + '\n')
        fh.write(seq[contig_range[0] : contig_range[1]] + '\n')
    fh.close()
    # running IgDetective
    igdetective_dir = os.path.join(output_dir, 'predicted_genes_' + locus)
    command_line = 'python py/IGDetective.py -i ' + fasta + ' -o ' + igdetective_dir + ' -m 1 -l ' + locus
    print('Running: ' + command_line)
    os.system(command_line + ' > ' + os.path.join(output_dir, 'predicted_genes_' + locus + '.out'))

def CombineIGGenes(genes_fasta, igdetective_tsv, output_fasta):
    nucl_seqs = set()
    if os.path.exists(genes_fasta):
        for r in SeqIO.parse(genes_fasta, 'fasta'):
            nucl_seqs.add(str(r.seq))
    df = pd.read_csv(igdetective_tsv, sep = '\t')
    for i in range(len(df)):
#        aa_seq = str(Seq(df['gene sequence'][i]).translate())
#        if aa_seq.find('*') == -1:
        nucl_seqs.add(df['gene sequence'][i])
    fh = open(output_fasta, 'w')
    for seq_idx, seq in enumerate(nucl_seqs):
        fh.write('>seq_' + str(seq_idx) + '\n' + seq + '\n')
    fh.close()

def AlignGenesIteratively(ref_gene_fasta, igdetective_tsv, genome_fasta, output_dir, gene_type, num_iter = 5):
    # aligning reference genes
    iter0_dir = os.path.join(output_dir, gene_type + '_iter0')
    extr_genes.main(genome_fasta, ref_gene_fasta, iter0_dir)
    iter0_fasta = os.path.join(iter0_dir, 'genes.fasta')
    # combining genes
    combined_fasta = os.path.join(output_dir, gene_type + '_combined.fasta')
    CombineIGGenes(iter0_fasta, igdetective_tsv, combined_fasta)
    # iterative alignment
    prev_iter_seqs = [r for r in SeqIO.parse(combined_fasta, 'fasta')]
    if len(prev_iter_seqs) == 0:
        return
    print('# combined genes: ' + str(len(prev_iter_seqs)))
    prev_fasta = combined_fasta
    for i in range(num_iter):
        print('== Iteration ' + str(i + 1) + '...')
        iter_dir = os.path.join(output_dir, gene_type + '_iter' + str(i + 1))
        extr_genes.main(genome_fasta, prev_fasta, iter_dir)
        curr_iter_fasta = os.path.join(iter_dir, 'genes.fasta')
        if not os.path.exists(curr_iter_fasta):
            print('gene file does not exist')
            break
        curr_iter_seqs = [r for r in SeqIO.parse(curr_iter_fasta, 'fasta')]
        print('# current genes: ' + str(len(curr_iter_seqs)) + ', # previous genes: ' + str(len(prev_iter_seqs)))
        if len(prev_iter_seqs) >= len(curr_iter_seqs) and i != 0:
            print('no new genes were detected, stopping the iterative search')
#            shutil.rmtree(iter_dir)
            break
        prev_iter_seqs = curr_iter_seqs
        prev_fasta = curr_iter_fasta

def ReadGeneDir(ig_gene_dir):
    files = os.listdir(ig_gene_dir)
    gene_dict = dict()
    for f in files:
        gene_type = f.split('.')[0]
        gene_dict[gene_type] = os.path.join(ig_gene_dir, f)
    return gene_dict

def main(genome_fasta, output_dir, ig_gene_dir):
    #### preparation
    CheckPythonVersionFatal()
    CheckMinimapFatal()
    PrepareOutputDir(output_dir)

    #### running IG gene alignments
    print('==== Aligning human IG genes...')
    alignment_dir = os.path.join(output_dir, 'alignments')
    os.mkdir(alignment_dir)
    AlignReferenceGenes(alignment_dir, genome_fasta, output_dir)
    
    #### identifying IG contigs
    print('==== Identifying IG contigs...')
    igcontig_dir = os.path.join(output_dir, 'ig_contigs')
    IdentifyIGContigs(igcontig_dir, alignment_dir, output_dir, genome_fasta)

    #### running IgDetective
    loci = ['IGH', 'IGK', 'IGL']
    for locus in loci:
        RunIgDetective(igcontig_dir, output_dir, locus)

    #### aligning IG genes
    ig_genes = ReadGeneDir(ig_gene_dir)
    gene_types = ['V', 'J']
    for locus in loci:
        for gene_type in gene_types:
            gene = locus + gene_type
            print('==== Iterative processing ' + gene + ' genes...')
            if gene not in ig_genes:
                continue
            ref_gene_fasta = ig_genes[gene]
            igdetective_tsv = os.path.join(os.path.join(output_dir, 'predicted_genes_' + locus), 'genes_' + gene_type + '.tsv')
            AlignGenesIteratively(ref_gene_fasta, igdetective_tsv, genome_fasta, output_dir, gene)

    #### the end
    print('Thank you for using IgDetective!')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('python run_iterative_igdetective.py genome.fasta output_dir')
        sys.exit(1)
    genome_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    ig_gene_dir = 'datafiles/combined_reference_genes' #sys.argv[3]
    main(genome_fasta, output_dir, ig_gene_dir)

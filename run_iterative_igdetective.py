import os
import sys
import subprocess
import shutil
import pandas as pd
from Bio import SeqIO

sys.path.append('py')

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

def GetRange(min_pos, max_pos, seq_len, max_len = 30000000):
    prefix_len = max_pos
    suffix_len = seq_len - min_pos
    gap = 1000000
    print(min_pos, max_pos, seq_len, prefix_len, suffix_len)
    if prefix_len > max_len and suffix_len > max_len:
        return (max(min_pos - gap, 0), min(prefix_len + gap, seq_len))
    if prefix_len < suffix_len:
        return (0, min(prefix_len + gap, seq_len))
    return (max(min_pos - gap, 0), seq_len)


def AlignIgGenes(genome_fasta, ig_gene_fasta, sam_file):
    os.system('minimap2 -a ' + genome_fasta + ' ' + ig_gene_fasta + ' -o ' + sam_file + ' > /dev/null')

def AlignReferenceGenes(genome_fasta, output_dir):
    print('==== Aligning human IG genes...')
    align_dir = os.path.join(output_dir, 'alignments')
    os.mkdir(align_dir)
    ref_gene_dict = dict()
    for f in os.listdir(ref_gene_dir):
        gene_type = f.split('.')[0]
        if gene_type not in ref_gene_dict:
            ref_gene_dict[gene_type] = os.path.join(ref_gene_dir, f)
    for gene_type in ref_gene_dict:
        print('Aligning ' + gene_type + ' genes (' + ref_gene_dict[gene_type] + ')...')
        AlignIgGenes(genome_fasta, ref_gene_dict[gene_type], os.path.join(align_dir, gene_type + '.sam'))
    return align_dir 

def IdentifyIGContigs(alignment_dir, output_dir, genome_fasta):
    igcontig_dir = os.path.join(output_dir, 'ig_contigs')
    print('==== Identifying IG contigs...')
    os.system('python py/analyze_matches.py ' + alignment_dir + ' ' + igcontig_dir + ' ' + genome_fasta)
    return igcontig_dir

def RunIgDetective(igcontig_dir, output_dir, locus = 'IGH'):
    txt = os.path.join(igcontig_dir, '__summary.txt')
    df = pd.read_csv(txt, sep = '\t')
    igh_df = df.loc[df['Locus'] == locus]
    fasta = os.path.join(output_dir, 'combined_contigs_' + locus + '.fasta')
    fh = open(fasta, 'w')
    contigs = set(igh_df['ContigID'])
    for c in contigs:
        c_df = igh_df.loc[igh_df['ContigID'] == c]
        contig_fasta = os.path.join(igcontig_dir, locus + '_' + c + '.fasta')
        if not os.path.exists(contig_fasta):
            print('WARN: ' + contig_fasta + ' does not exist')
            continue
        seq = ''
        seq_id = ''
        for r in SeqIO.parse(contig_fasta, 'fasta'):
            seq = str(r.seq)
            seq_id = r.id
        min_pos = min(c_df['Position'])
        max_pos = max(c_df['Position'])
        contig_range = GetRange(min_pos, max_pos, len(seq))
        print(contig_range, contig_range[1] - contig_range[0])
        fh.write('>' + seq_id + '|START:' + str(contig_range[0]) + '|END:' + str(contig_range[1]) + '\n')
        fh.write(seq[contig_range[0] : contig_range[1]] + '\n')
    fh.close()
    # running IgDetective
    igdetective_dir = os.path.join(output_dir, 'predicted_genes')
    command_line = 'python IGDetective.py -i ' + fasta + ' -o ' + igdetective_dir + ' -m 1'
    print('Running: ' + command_line)
    os.system(command_line)

def main(genome_fasta, output_dir):
    #### preparation
    CheckPythonVersionFatal()
    CheckMinimapFatal()
    PrepareOutputDir(output_dir)

    #### running IG gene alignments
    alignment_dir = AlignReferenceGenes(genome_fasta, output_dir)
    
    #### identifying IG contigs
    igcontig_dir = IdentifyIGContigs(alignment_dir, output_dir, genome_fasta)

    #### running IgDetective
    RunIgDetective(igcontig_dir, output_dir)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('python run_iterative_igdetective.py genome.fasta output_dir')
        sys.exit(1)
    genome_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    main(genome_fasta, output_dir)

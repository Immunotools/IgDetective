import os
import sys
import shutil
import pandas as pd
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq

def FindAlignmentRange(gene_alignment):
    first_pos = ''
    for pos, c in enumerate(gene_alignment):
        if c not in ['.', ' ']:
            first_pos = pos
            break
    last_pos = ''
    for i in range(len(gene_alignment)):
        pos = len(gene_alignment) - i - 1
        if gene_alignment[pos] not in ['.', ' ']:
            last_pos = pos + 1
            break
    return [first_pos, last_pos]

def GetNumMatches(alignment):
    splits = str(alignment).split('\n')
    query_alignment = splits[0].upper()
    gene_alignment = splits[2].upper()
#    print(query_alignment + '\n' + gene_alignment + '\n===')
    alignment_range = FindAlignmentRange(gene_alignment)
    query_subalignment = query_alignment[alignment_range[0] : alignment_range[1]].upper()
    matches = [i for i in range(min(len(gene_alignment), len(query_alignment))) if query_alignment[i] == gene_alignment[i]]
    return len(matches)

class Alignment:
    def __init__(self):
        self.gene_seq = ''
        self.gene_id = ''
        self.pi = 0

    def Initiate(self, gene_seq, gene_id, pi):
        self.gene_seq = gene_seq
        self.gene_id = gene_id
        self.pi = pi

    def Empty(self):
        return self.gene_seq == '' or self.gene_seq[0] == ' '

def ComputeAlignment(aligner, query_list, gene_seqs):
    best_alignment = ''
    best_num_matches = 0
    best_gene = ''
    for query in query_list:
        for gene in gene_seqs:
            alignments = aligner.align(query, gene.seq)
            if len(alignments) == 0:
                continue
            alignment = alignments[0]
            num_matches = GetNumMatches(alignment)
            if num_matches >= best_num_matches:
                best_num_matches = num_matches
                best_alignment = alignment
                best_gene = gene.id
    splits = str(best_alignment).split('\n')
    if len(splits) == 1:
        return Alignment()
    query_alignment = splits[0].upper()
    gene_alignment = splits[2].upper()
    align_range = FindAlignmentRange(gene_alignment)
    alignment = Alignment()
    alignment.Initiate(''.join([aa for aa in query_alignment[align_range[0] : align_range[1]] if aa != '-']), best_gene, float(best_num_matches) / (align_range[1] - align_range[0]))
    return alignment

def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

def ProcessSamFile(sam_file):
    position_dict = dict()
    for l in open(sam_file).readlines():
        l = l.strip()
        if l[0] == '@':
            continue
        splits = l.split()
        contig_id = splits[2]
        pos = int(splits[3])
        if contig_id == '*':
            continue
        if contig_id not in position_dict:
            position_dict[contig_id] = []
        if pos not in position_dict[contig_id]:
            position_dict[contig_id].append(pos)
    return position_dict

def main(genome_fasta, gene_fasta, output_dir):
    PrepareOutputDir(output_dir)

    print('Running minimap...')
    print('Alignment of IG genes ' + gene_fasta + ' to ' + genome_fasta)
    sam_file = os.path.join(output_dir, 'alignment.sam')
    os.system('minimap2 -a ' + genome_fasta + ' ' + gene_fasta + ' -o ' + sam_file + '> /dev/null 2>&1')

    print('Processing SAM file...')
    position_dict = ProcessSamFile(sam_file)
    if len(position_dict) == 0:
        print('no matches were found')
        return

    contig_dict = dict()
    for r in SeqIO.parse(genome_fasta, 'fasta'):
        if r.id in position_dict:
            contig_dict[r.id] = str(r.seq)

    genes = []
    processed_families = set()
    for r in SeqIO.parse(gene_fasta, 'fasta'):
        family = r.id.split('*')[0].split('-')[0].split('S')[0]
        if family not in processed_families:
            r.seq = str(r.seq).upper()
            genes.append(r)
            processed_families.add(family)

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = 0
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1

    gene_len = 400 #max([len(gene) for gene in genes])
    df = {'Contig' : [], 'Pos' : [], 'Seq' : [], 'AASeq' : [], 'PI' : [], 'BestHit' : [], 'Productive' : []}
    for c_id in position_dict:
        contig_seq = contig_dict[c_id]
        prev_pos = -1
        for pos in sorted(position_dict[c_id]):
#            print(pos)
            if pos - prev_pos <= gene_len:
#                print('redundant position: ' + str(pos))
                continue
            fragment = contig_seq[max(0, pos - gene_len) : min(len(contig_seq), pos + gene_len)]
            fragment_rc = str(Seq(fragment).reverse_complement())
            alignment = ComputeAlignment(aligner, [fragment, fragment_rc], genes)
            if alignment.Empty():
#                print('empty alignment')
                continue
            aa_seq = str(Seq(alignment.gene_seq).translate())
#            if aa_seq.find('*') != -1:
#                print('non-productive gene')
#                continue
#            print('==== ' + c_id + ', pos: ' + str(pos) + ', gene: ' + alignment.gene_id + ', PI: ' + str(alignment.pi))
#            print(alignment.gene_seq)
#            print(aa_seq)
            df['Contig'].append(c_id)
            df['Pos'].append(pos)
            df['Seq'].append(alignment.gene_seq)
            df['AASeq'].append(aa_seq)
            df['PI'].append(alignment.pi)
            df['BestHit'].append(alignment.gene_id)
            df['Productive'].append(aa_seq.find('*') == -1)
            prev_pos = pos

    df = pd.DataFrame(df)
    df.to_csv(os.path.join(output_dir, 'genes.tsv'), index = False, sep = '\t')

    print('# detected genes: ' + str(len(df)))

    fh = open(os.path.join(output_dir, 'genes.fasta'), 'w')
    for i in range(len(df)):
        fh.write('>Contig:' + df['Contig'][i] + '|Pos:' + str(df['Pos'][i]) + '\n' + df['Seq'][i] + '\n')
    fh.close()

if __name__ == '__main__':
    if len(sys.argv) != 4:
       print('Invalid arguments')
       print('python extract_aligned_genes.py genome.fasta reference_IG_genes.fasta output_dir')
       sys.exit(1)
    genome_fasta = sys.argv[1]
    gene_fasta = sys.argv[2]
    output_dir = sys.argv[3]
    main(genome_fasta, gene_fasta, output_dir)

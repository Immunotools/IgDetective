import os
import sys
import gzip
import shutil
import pandas as pd
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq

class BioAlign:
    def __init__(self, alignment):
        self.alignment = alignment
        self.query_align = self.alignment[0].upper()
        self.gene_align = self.alignment[1].upper()
        self.gene_range = self._ComputeGeneRange()
        ####
        self.num_matches = -1
        self.query_seq = ''

    def _ComputeGeneRange(self):
        start_pos = 0
        for i in range(len(self.gene_align)):
            if self.gene_align[i] != '-':
                start_pos = i
                break
        end_pos = len(self.gene_align)
        for i in range(len(self.gene_align)):
            pos = len(self.gene_align) - i - 1
            if self.gene_align[pos] != '-':
                end_pos = pos + 1
                break
        return start_pos, end_pos

    def NumMatches(self):
        if self.num_matches == -1:
            for i in range(self.gene_range[0], self.gene_range[1]):
                if self.query_align[i] == self.gene_align[i]:
                    self.num_matches += 1
        return self.num_matches

    def PI(self):
        return self.NumMatches() / len(self) * 100

    def QuerySeq(self):
        if self.query_seq == '':
            nucls = []
            for n in self.query_align[self.gene_range[0] : self.gene_range[1]]:
                if n != '-':
                    nucls.append(n)
            self.query_seq = ''.join(nucls)
        return self.query_seq

    def __len__(self):
        return self.gene_range[1] - self.gene_range[0]

    def AlignmentRange(self):
        return self.gene_range

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

def SetupAligner(aligner):
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = 0
    aligner.target_internal_open_gap_score = -2
    aligner.target_internal_extend_gap_score = -1
    aligner.target_end_open_gap_score = -2
    aligner.target_end_extend_gap_score = -1
    aligner.query_internal_open_gap_score = -2
    aligner.query_internal_extend_gap_score = -1
    aligner.query_end_open_gap_score = 0
    aligner.query_end_extend_gap_score = 0

def ComputeAlignment(aligner, query_list, strand_list, gene_seqs):
    best_alignment = ''
    best_pi = 0
    best_gene = ''
    best_strand = ''
    for query, strand in zip(query_list, strand_list):
        for gene in gene_seqs:
            alignments = aligner.align(query, gene.seq)
            if len(alignments) == 0:
                continue
            alignment = BioAlign(alignments[0])
            if alignment.PI() >= best_pi:
                best_pi = alignment.PI()
                best_alignment = alignment
                best_gene = gene.id
                best_strand = strand
    if len(best_alignment) == 0:
        return Alignment(), ''
    alignment = Alignment()
    alignment.Initiate(best_alignment.QuerySeq(), best_gene, best_pi)
    return alignment, best_strand

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
    
    if genome_fasta.endswith('.gz'):
        genome_handle = gzip.open(genome_fasta, 'rt')
    else:
        genome_handle = open(genome_fasta, 'r')
    contig_dict = dict()
    for r in SeqIO.parse(genome_handle, 'fasta'):
        if r.id in position_dict:
            contig_dict[r.id] = str(r.seq)
    genome_handle.close()

    genes = []
    for r in SeqIO.parse(gene_fasta, 'fasta'):
        r.seq = str(r.seq).upper()
        genes.append(r)

    aligner = Align.PairwiseAligner()
    SetupAligner(aligner)

    gene_len = 400 #max([len(gene) for gene in genes])
    df = {'Contig' : [], 'Pos' : [], 'Seq' : [], 'AASeq' : [], 'PI' : [], 'BestHit' : [], 'Productive' : [], 'Strand' : []}
    for c_id in position_dict:
        contig_seq = contig_dict[c_id]
        prev_pos = -1
        for pos in sorted(position_dict[c_id]):
            if pos - prev_pos <= gene_len:
                continue
            fragment = contig_seq[max(0, pos - gene_len) : min(len(contig_seq), pos + gene_len)]
            fragment_rc = str(Seq(fragment).reverse_complement())
            alignment, strand = ComputeAlignment(aligner, [fragment, fragment_rc], ['+', '-'], genes)
            if alignment.Empty():
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
            df['Strand'].append(strand)
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

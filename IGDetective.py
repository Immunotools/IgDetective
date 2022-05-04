import os,getopt
import sys
import csv
import itertools
import pickle

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align

import numpy as np

from multiprocessing import Pool

#INITIALIZE VARIABLES
V = 'V'
D = 'D'
J = 'J'
DR = 'D_right'
DL = 'D_left'
FWD = '+'
REV = '-'
GENE_TYPES = [V,D,J]
GENE_TYPES_TOFIND = [V,D,J]
SPACER_LENGTH = {V:23 , DL:12 , DR:12 , J:23}
GENE_LENGTH = {V: 350 , J: 70}
ALIGNMENT_EXTENSION = {V:REV , J:FWD, D:None}

ALIGNER = Align.PairwiseAligner()
ALIGNER.mode = 'local'




#READ DATAFILES
try:
    with open('datafiles/motifs', 'rb') as f:   
        VALID_MOTIFS = pickle.load(f)
except:
    print("Error: could not find the input data files. Please make sure the IGDetective.py file and datafiles folder are in the same directory")

#PARSE COMMAND LINE ARGUMENTS
argumentList = sys.argv[1:]
options = "hi:o:m:rg:"
long_options = ["help","input_file=", "output_directory=", "multi_process=", "rss_only" , "genes_type="]
force_output = True
received_input = False
RSS_MODE = False
help_flag = False
NUM_THREADS = 1
try:
    arguments, values = getopt.getopt(argumentList, options, long_options)
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-h", "--help"):
            print ("Diplaying Help")
            print("Flags and their usage :")
            print("-h , --help : Get this message")
            print("-i, --input_file : provide a fasta file for gene detection")
            print("-o, --output_directory : (optional) provide an output directory for generated results. Default location is in the parent directory of the input file")
            print("-m, --multi_process : (optional) provide number of parallel processing units if available. Default is 1")
            print("-r, --rss_only : (optional) switch to RSS finding mode")
            print("-g, --genes_type : (optional) specify which genes (v,d,j) to find. Eg: vdj, d, vj, jv. Default is vdj")
            help_flag = True
            
        elif currentArgument in ("-i", "--input_file"):
            INPUT_PATH = str(currentValue)
            received_input = True
        
        elif currentArgument in ("-o", "--output_directory"):
            OUTPUT_PATH = str(currentValue)
            force_output = False

        elif currentArgument in ("-m", "--multi_process"):
            NUM_THREADS = int(currentValue)

        elif currentArgument in ("-r", "--rss_only"):
            RSS_MODE = True

        elif currentArgument in ("-g", "--genes_type"):
            GENE_TYPES_TOFIND = set(list(str(currentValue).upper()))
            for g in GENE_TYPES_TOFIND:
                if g not in GENE_TYPES:
                    raise NameError('gene types must be from v,d or j')


    if not received_input and not help_flag:
        raise NameError('no input file was given')
                
except getopt.error as err:
    print (str(err))
    sys.exit(0)

if force_output == True:
    OUTPUT_PATH = ".".join(INPUT_PATH.split('.')[:-1])
if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)
    
#create signal types from gene types
SIGNAL_TYPES = []
if V in GENE_TYPES_TOFIND:
    SIGNAL_TYPES.append(V)
if J in GENE_TYPES_TOFIND:
    SIGNAL_TYPES.append(J)
if D in GENE_TYPES_TOFIND:
    SIGNAL_TYPES.extend([DL,DR])


#READ INPUT FASTA FILE
input_seq_dict= {rec.id : rec.seq for rec in SeqIO.parse(INPUT_PATH, "fasta")}

canonical_genes = {V : {} , J : {}}
for gene in (V,J):
    file_path = 'datafiles/human_{}.fasta'.format(gene)
    canonical_genes[gene] = {rec.id : rec.seq for rec in SeqIO.parse(file_path, "fasta")}

#DEFINE RSS FINDING METHODS
#Find indexes of valid motifs
def find_valid_motif_idx(locus,motifs,k):
    motif_idx = []
    for i in range(0,len(locus)- int(k) +1):
        candidate = locus[i:i+int(k)].upper()
        if candidate in motifs:
            motif_idx.append(i)
    return set(motif_idx)

#return idx of heptamer and nonamer
def find_valid_rss(heptamer_idx, nonamer_idx, sig_type, strand, seq_length):
    rss_idx = []
    spacer = SPACER_LENGTH[sig_type]
    
    #set the 5' appearing k-mer
    if sig_type == V or sig_type == DR:
        k_first = 7
        first_set = heptamer_idx
        second_set = nonamer_idx
    elif sig_type == J or sig_type == DL:
        k_first = 9
        first_set = nonamer_idx
        second_set = heptamer_idx
    
    #search for spacer separation between heptamer and nonamer    
    for idx in first_set:
        if spacer + idx + k_first in second_set:
            rss_idx.append((idx ,spacer  + idx + k_first))
        elif spacer - 1 + idx + k_first in second_set:
            rss_idx.append((idx , spacer - 1 + idx + k_first))
        elif spacer + 1 + idx + k_first in second_set:
            rss_idx.append((idx , spacer + 1 + idx + k_first))
    
    #set tuple to start with heptamer only
    if sig_type == J or sig_type == DL:
        rss_idx = [(x[1],x[0]) for x in rss_idx]

    if strand == REV:
        rss_idx = [(seq_length-x[0]-7 , seq_length-x[1]-9) for x in rss_idx]

    return rss_idx

#combine data of heptamer and nonamer indexes
def get_contigwise_rss(sig_type,strand,parent_seq):
    parallel_heptamers = []
    parallel_nonamers = []
    parallel_rss = []
    
    #find valid heptamer and nonamers motifs
    for i,contigs in enumerate(list(parent_seq.keys())):
        if strand == FWD:
            sequence = str(parent_seq[contigs])
        elif strand == REV:
            sequence = str(parent_seq[contigs].reverse_complement())
            
        parallel_heptamers.append((sequence, VALID_MOTIFS[sig_type]['7'], 7))
        parallel_nonamers.append((sequence, VALID_MOTIFS[sig_type]['9'], 9))

    p = Pool(NUM_THREADS)    
    heptamer_resultset = p.starmap(find_valid_motif_idx,parallel_heptamers)
    nonamer_resultset = p.starmap(find_valid_motif_idx,parallel_nonamers)

    #combine valid heptamer and nonamer motifs
    for i,contig in enumerate(list(parent_seq.keys())):
        L = len(parent_seq[contig])
        parallel_rss.append((heptamer_resultset[i], nonamer_resultset[i], sig_type, strand, L)) 
    result = p.starmap(find_valid_rss , parallel_rss)
    rss_resultset = {contig : result[i] for i,contig in enumerate(list(parent_seq.keys()))}

    return rss_resultset

#D_left(D_right) idx is of the form "input_rss_info['D_left(D_right)']"
def combine_D_RSS(D_left_idx , D_right_idx, input_seq_dict, strand, Dgene_len = 150):
    rss_resultset = {contigs : [] for contigs in input_seq_dict.keys()}
    for contig in input_seq_dict:
        for dr in D_right_idx[contig]:
            for dl in D_left_idx[contig]:
                if strand == FWD:
                    right_heptamer = dr[0]
                    left_heptamer = dl[0]
                elif strand == REV:
                    right_heptamer = dl[0]
                    left_heptamer = dr[0]
       
                if right_heptamer-(left_heptamer+7) <= Dgene_len and left_heptamer < right_heptamer:
                    rss_resultset[contig].append((dl[0] , dl[1], dr[0], dr[1]))
    return rss_resultset

#write RSS details to file
def write_rss_to_file(filepath, rss_idx , parent_seq, min_heptamer = None , max_heptamer = None):
    detected_RSS_info = [['7-mer index' , '9-mer index' , '7-mer' , '9-mer', 'reference contig', 'strand']]
    for strand in [FWD , REV]:
        for contigs in parent_seq:
            for pair in rss_idx[strand][contigs]:
                hepta_index = pair[0]
                nona_index = pair[1]

                if strand == FWD:
                    hepta = parent_seq[contigs][hepta_index:hepta_index+7].upper()
                    nona = parent_seq[contigs][nona_index:nona_index+9].upper()
                
                elif strand == REV:
                    hepta = parent_seq[contigs][hepta_index:hepta_index+7].reverse_complement().upper()
                    nona = parent_seq[contigs][nona_index:nona_index+9].reverse_complement().upper()

                accept = True
                if min_heptamer and hepta_index< min_heptamer:
                    accept = False
                if max_heptamer and L-hepta_index < max_heptamer:
                    accept = False

                if accept:
                    rss_info = [hepta_index , nona_index, hepta, nona, contigs, strand]
                    detected_RSS_info.append(rss_info)
                    
    with open(filepath, "w", newline="") as f:
        writer =csv.writer(f , delimiter = '\t')
        writer.writerows(detected_RSS_info)

#METHODS FOR SEQUENCE ALIGNMENT

# Extract S fragmnt
def extract_s_fragment(index, extract_dir, length, parent_seq):
    #assumes index is one index after end / before start of S fragment
    if extract_dir == FWD:
        fragment = parent_seq[index:index+length]
    elif extract_dir == REV:
        fragment = parent_seq[index-length:index]
    return fragment

def get_s_fragment_from_RSS(gene, strand):
    s_fragments = {contig : [] for contig in input_rss_info[gene][strand]}
    if (gene == V and strand == FWD) or (gene == J and strand == REV):
        for contig in input_rss_info[gene][strand]:
            for rss in input_rss_info[gene][strand][contig]:
                fragment = extract_s_fragment(rss[0], REV, GENE_LENGTH[gene], input_seq_dict[contig])
                s_fragments[contig].append(fragment)
            if strand == REV:
                s_fragments[contig] = [x.reverse_complement() for x in s_fragments[contig]]
                
    elif (gene == J and strand == FWD) or (gene == V and strand == REV):
        for contig in input_rss_info[gene][strand]:
            for rss in input_rss_info[gene][strand][contig]:
                fragment = extract_s_fragment(rss[0] + 7, FWD, GENE_LENGTH[gene], input_seq_dict[contig])
                s_fragments[contig].append(fragment)
            if strand == REV:
                s_fragments[contig] = [x.reverse_complement() for x in s_fragments[contig]]
    
    elif gene == D:
        for contig in input_rss_info[gene][strand]:
            for rss in input_rss_info[gene][strand][contig]:
                if strand == FWD:
                    index = rss[0]+7
                    gene_len = rss[2] - index
                elif strand == REV:
                    index = rss[2]+7
                    gene_len = rss[0] - index
                fragment = extract_s_fragment(index, FWD,gene_len, input_seq_dict[contig])
                s_fragments[contig].append(fragment)
            if strand == REV:
                s_fragments[contig] = [x.reverse_complement() for x in s_fragments[contig]]

    return s_fragments

#define alignment score based on scheme
def set_aligner(scheme):
    aligner_scoring = {'MAXK' :[1,-10000,-10000,-10000],
                   'AFFINE' :[1,0,-1,-0.5]}
    ALIGNER.match_score = aligner_scoring[scheme][0]
    ALIGNER.mismatch_score = aligner_scoring[scheme][1]
    ALIGNER.open_gap_score = aligner_scoring[scheme][2]
    ALIGNER.extend_gap_score = aligner_scoring[scheme][3]


#count length of alignment
def count_alignment_len(al, gene_al, extend_alignment):
    start,end = FindAlignmentRange(al, gene_al,  extend_alignment)
    return 1+end-start

#find start and end of alignment
def FindAlignmentRange(al, gene_al, extend_alignment = None):
    if '-' in al:
        start,end = min(al.find('|'),al.find('-')) , max(al.rfind('|'),al[1].rfind('-'))
    else:
        start,end = al.find('|') , al.rfind('|')
        
    if extend_alignment == REV:
        residual = len(gene_al[:start]) - gene_al[:start].count(' ')
        start -= residual
    elif extend_alignment == FWD:
        residual = len(gene_al[end:]) - gene_al[end:].count(' ')
        end += residual
        
    return start,end

#get numnber of matches in an alignment
def GetNumMatches(alignment, extend_alignment = None):
    splits = str(alignment).split('\n')
    seq_A_alignemnt = splits[0].upper()
    seq_B_alignment = splits[2].upper()
    coded_alignment = splits[1].upper()
    start,end = FindAlignmentRange(coded_alignment, seq_B_alignment,  extend_alignment)
    matches = (coded_alignment[start : end+1]).count('|')

    return matches

#align 2 strings
def ComputeAlignment(seq_A, seq_B, extend_alignment = None):
    query = seq_A
    query_rc = query.reverse_complement()
    alignment = ALIGNER.align(query, seq_B)[0]
    alignment_rc = ALIGNER.align(query_rc, seq_B)[0]
    fwd_matches  = GetNumMatches(alignment, extend_alignment)
    rev_matches = GetNumMatches(alignment_rc, extend_alignment)
    if  fwd_matches > rev_matches:
        return alignment, '+', fwd_matches
    else:
        return alignment_rc, '-', rev_matches

def align_fragment_to_genes(fragments, canon_genes, scoring_scheme, gene):
    set_aligner(scoring_scheme)
    parallel_alignments = []
    pis = []
    alignment_lens = []
    
    for i in range(0,len(fragments)):
        for j in range(0,len(canon_genes)):
            seq_A = fragments[i]
            seq_B = canon_genes[j]
            parallel_alignments.append((seq_A,seq_B,ALIGNMENT_EXTENSION[gene]))
            
    p = Pool(NUM_THREADS)
    alignment_results = p.starmap(ComputeAlignment , parallel_alignments)
    for a in alignment_results:
        alignment = a[0]
        strand = a[1]
        num_matches = a[2]
        splits = str(alignment).split('\n')
        alig_len = count_alignment_len(splits[1],splits[2], ALIGNMENT_EXTENSION[gene])
        perc_identity = num_matches* 100 / alig_len
        pis.append(perc_identity)
        alignment_lens.append(alig_len)
        
    pi_mat = np.zeros((len(fragments), len(canon_genes)))
    alig_len_mat = np.zeros((len(fragments), len(canon_genes)))
    
    for i in range(0,len(fragments)):
        for j in range(0,len(canon_genes)):
            pi_mat[i][j] = pis.pop(0)
            alig_len_mat[i][j] = alignment_lens.pop(0)
            
    return pi_mat , alig_len_mat


#Evaluate and print genes
def extract_genes(parent_seq, gene, rss_idx, fragments, fragment_alignments):
    final_genes = []
    if gene == D:
        for strand in rss_idx:
            for contig in rss_idx[strand]:
                for e in rss_idx[strand][contig]:
                    if strand == FWD:
                        lh , ln , rh, rn = parent_seq[contig][e[0]:e[0]+7].upper() , parent_seq[contig][e[1]:e[1]+9].upper(), \
                        parent_seq[contig][e[2]:e[2]+7].upper(), parent_seq[contig][e[3]:e[3]+9].upper()
                        gs, ge, predicted_gene = e[0]+7, e[2]-1, parent_seq[contig][e[0]+7:e[2]].upper()
                    elif strand == REV:
                        lh , ln , rh, rn = parent_seq[contig][e[0]:e[0]+7].reverse_complement().upper() , parent_seq[contig][e[1]:e[1]+9].reverse_complement().upper(), \
                        parent_seq[contig][e[2]:e[2]+7].reverse_complement().upper(), parent_seq[contig][e[3]:e[3]+9].reverse_complement().upper()
                        gs, ge ,predicted_gene =  e[2]+7, e[0]-1, parent_seq[contig][e[2]+7:e[0]].reverse_complement().upper()      
                    final_genes.append([contig, strand, e[0],e[1], lh, ln, e[2],e[3],rh,rn,gs, ge, predicted_gene])
                    
    return final_genes

def print_predicted_genes(filepath, gene, predictions):
    if gene == D:
        detected_gene_info = [['reference contig', 'strand', 'left heptamer index' , 'left nonamer index' ,\
                               'left heptamer' , 'left nonamer', 'right heptamer index' , 'right nonamer index' ,\
                               'right heptamer' , 'right nonamer', 'start of gene', 'end of gene', 'gene sequence']]
        detected_gene_info.extend(predictions)
    with open(filepath, "w", newline="") as f:
        writer =csv.writer(f , delimiter = '\t')
        writer.writerows(detected_gene_info)


#FIND RSS IN INPUT FASTA
print("Finding candidate RSS...",end =" ")
input_rss_info = {st : {strand : get_contigwise_rss(st,strand, input_seq_dict) for strand in (FWD , REV)} for st in SIGNAL_TYPES}
if D in GENE_TYPES_TOFIND:
    input_rss_info[D] = { strand: combine_D_RSS(input_rss_info[DL][strand] , input_rss_info[DR][strand], input_seq_dict , strand)\
                       for strand in (FWD, REV)}
print("Done")
if RSS_MODE:
    for st in SIGNAL_TYPES:
        write_rss_to_file('{}/rss_{}.csv'.format(OUTPUT_PATH, st), input_rss_info[st], input_seq_dict)
    sys.exit(0)


#create and alingn S fragments
s_fragments = {g : {strand : get_s_fragment_from_RSS(g,strand) for strand in (FWD, REV)} for g in GENE_TYPES_TOFIND}
fragments_to_align = {gene : [] for gene in GENE_TYPES_TOFIND}
for gene in GENE_TYPES_TOFIND:
    for strand in (FWD,REV):
        for contig in s_fragments[gene][strand]:
            fragments_to_align[gene].extend(s_fragments[gene][strand][contig])

print("Aligning candidate genes...",end =" ")          
s_fragment_alignment = {gene : { strand : {contig : [] for contig in s_fragments[gene][strand]} for strand in (FWD,REV)} for gene in GENE_TYPES_TOFIND}
for gene in GENE_TYPES_TOFIND:
    if gene == D:
        continue
    pi_mat , _ = align_fragment_to_genes(fragments_to_align[gene], list(canonical_genes[gene].values()) , 'AFFINE', gene)
    _ , maxk_mat = align_fragment_to_genes(fragments_to_align[gene], list(canonical_genes[gene].values()) , 'MAXK', gene)
    k = 0
    for strand in (FWD,REV):
        for contig in s_fragments[gene][strand]:
            for sequence in s_fragments[gene][strand][contig]:
                best_alignment_index = np.argmax(pi_mat[k])
                pi = pi_mat[k][best_alignment_index]
                maxk = maxk_mat[k][best_alignment_index]
                k+=1
                s_fragment_alignment[gene][strand][contig].append((best_alignment_index,pi,maxk))
print("Done")

#Print genes to tsv file
print_predicted_genes('{}/genes_{}.csv'.format(OUTPUT_PATH, D) , D, extract_genes(input_seq_dict, D, input_rss_info[D], None, None))
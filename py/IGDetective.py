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
from multiprocessing import get_context

import extract_aligned_genes as align_utils


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
#INITIALIZE VARIABLES, default = IGH
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
GENE_LENGTH = {V: 350 , J: 70, D:150}
ALIGNMENT_EXTENSION = {V:REV , J:FWD, D:None}
PI_CUTOFF = {'strict' : {V: 70, J: 70} , 'relax': {V: 60 ,J: 65}}
MAXK_CUTOFF = {V: 15 , J: 11}

ALIGNER = Align.PairwiseAligner()

def InitializeVariables(locus):
    if locus == 'IGK':
        SPACER_LENGTH = {V:12, DL:23, DR:12, J:23}
        GENE_TYPES_TOFIND = [V,J]
    if locus == 'IGL':
        SPACER_LENGTH = {V:23, DL:12, DR:23, J:12}    
        GENE_TYPES_TOFIND = [V,J]
    if locus == 'TRA':
        SPACER_LENGTH = {V:23, J:12}
        GENE_TYPES_TOFIND = [V,J]
    if locus == 'TRB':
        SPACER_LENGTH = {V:23, J:12, DL:12, DR:23}
        GENE_TYPES_TOFIND = [V,D,J]
    if locus == 'TRD':
        SPACER_LENGTH = {V:23, J:12, DL:12, DR:23}
        GENE_TYPES_TOFIND = [V,D,J]

#READ DATAFILES
try:
    motifs_file = os.path.join(SCRIPT_DIR[:-3], 'datafiles', 'motifs')
    with open(motifs_file, 'rb') as f:   
        VALID_MOTIFS = pickle.load(f)
except:
    print("Error: could not find the input data files. Please make sure the IGDetective.py file and datafiles folder are in the same directory")

#PARSE COMMAND LINE ARGUMENTS
argumentList = sys.argv[1:]
options = "hi:o:m:rg:l:"
long_options = ["help","input_file=", "output_directory=", "multi_process=", "rss_only" , "genes_type=", "locus="]
force_output = True
received_input = False
LOCUS = 'IGH'
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
            print("-l, --locus : immunoglobulin locus IGH, IGK, or IGL. Default is IGH")
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

        elif currentArgument in ("-l", "--locus"):
            if currentValue not in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG']:
                print('Incorrect locus argument: ' + currentValue)
                sys.exit(1)
            LOCUS = currentValue
            InitializeVariables(LOCUS)

        elif currentArgument in ("-m", "--multi_process"):
            NUM_THREADS = int(currentValue)

        elif currentArgument in ("-r", "--rss_only"):
            RSS_MODE = True

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
    file_path = os.path.join(SCRIPT_DIR[:-3], 'datafiles', 'combined_reference_genes', LOCUS + gene + '.fa') #'datafiles/human_{}.fasta'.format(gene)
    canonical_genes[gene] = {rec.id : rec.seq.upper() for rec in SeqIO.parse(file_path, "fasta")}

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

#    p = Pool(NUM_THREADS)    
    p = get_context("fork").Pool(NUM_THREADS)
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
def combine_D_RSS(D_left_idx , D_right_idx, input_seq_dict, strand, Dgene_len = GENE_LENGTH[D]):
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
    align_utils.SetupAligner(ALIGNER)

#align 2 strings
def ComputeAlignment(seq_A, seq_B, extend_alignment = None):
    query = Seq(str(seq_A).upper())
    query_rc = query.reverse_complement()
    alignment = ALIGNER.align(query, seq_B)[0]
    alignment_rc = ALIGNER.align(query_rc, seq_B)[0]
    fwd_matches  = align_utils.BioAlign(alignment).NumMatches() 
    rev_matches = align_utils.BioAlign(alignment_rc).NumMatches() 
    if fwd_matches > rev_matches:
        return alignment, '+', fwd_matches
    else:
        return alignment_rc, '-', rev_matches

def align_fragment_to_genes(fragments, canon_genes, scoring_scheme, gene):
    set_aligner(scoring_scheme)
    parallel_alignments = []
    pis = []
    alignment_lens = []

# non-parallel version
#    alignment_results = []
#    for i in range(0, len(fragments)):
#        for j in range(0, len(canon_genes)):
#            seq_A = fragments[i]
#            seq_B = canon_genes[j]
#            if len(seq_A) == 0:
#                seq_A = 'A' * len(seq_B)
#            alignment_results.append(ComputeAlignment(seq_A, seq_B))    

    for i in range(0,len(fragments)):
        for j in range(0,len(canon_genes)):
            seq_A = fragments[i]
            seq_B = canon_genes[j]
            if len(seq_A) == 0:
                seq_A = 'A' * len(seq_B) #canon_genes[j]
            parallel_alignments.append((seq_A,seq_B,ALIGNMENT_EXTENSION[gene]))
            
    p = get_context("fork").Pool(NUM_THREADS)
    alignment_results = p.starmap(ComputeAlignment, parallel_alignments)
    for a in alignment_results:
        alignment = a[0]
        strand = a[1]
        num_matches = a[2]
        bio_align = align_utils.BioAlign(alignment)
        pis.append(bio_align.PI())
        alignment_lens.append(len(bio_align))
        
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
        
    elif gene == V or gene == J:
        for strand in fragment_alignments:
            for contig in fragment_alignments[strand]:
                r = rss_idx[strand][contig]
                f = fragments[strand][contig]
                c = list(canonical_genes[gene].keys()) 
                for i,e in enumerate(fragment_alignments[strand][contig]):
                    if e[1] >= PI_CUTOFF['strict'][gene] or (e[1] >= PI_CUTOFF['relax'][gene] and e[2] >= MAXK_CUTOFF[gene]):
                        a = ComputeAlignment(f[i], canonical_genes[gene][c[e[0]]], ALIGNMENT_EXTENSION)
                        alignment = align_utils.BioAlign(a[0])
                        alig_direction = a[1]
                        start,end = alignment.AlignmentRange()
                        predicted_gene = alignment.QuerySeq() 

                        if strand == FWD:
                            if gene == V:
                                ge = r[i][0] -1
                                gs = r[i][0]- GENE_LENGTH[gene] + start
                            elif gene == J:
                                gs = r[i][0] + 7
                                ge = gs + end
                            #    predicted_gene = parent_seq[contig][gs:ge].upper()
                            h = parent_seq[contig][r[i][0]:r[i][0]+7].upper()
                            n = parent_seq[contig][r[i][1]:r[i][1]+9].upper()

                        elif strand == REV:
                            if gene == V:
                                gs = r[i][0] + 7
                                ge = gs + GENE_LENGTH[gene] - start -1
                            #    predicted_gene = parent_seq[contig][gs:ge+1].reverse_complement().upper()
                            elif gene == J:
                                ge = r[i][0] -1
                                gs = r[i][0] - end
                            #    predicted_gene = parent_seq[contig][gs:ge+1].reverse_complement().upper()
                            h = parent_seq[contig][r[i][0]:r[i][0]+7].reverse_complement().upper()
                            n = parent_seq[contig][r[i][1]:r[i][1]+9].reverse_complement().upper()
                        hi = r[i][0]
                        ni = r[i][1]
                        human_neighbour = c[e[0]]
                        
                        final_genes.append([contig, strand, hi, ni, h, n, gs, ge, human_neighbour, alig_direction, int(round(e[1])), e[2], predicted_gene])
   
    return final_genes

def print_predicted_genes(filepath, gene, predictions):
    if gene == D:
        detected_gene_info = [['reference contig', 'strand', 'left heptamer index' , 'left nonamer index' ,\
                               'left heptamer' , 'left nonamer', 'right heptamer index' , 'right nonamer index' ,\
                               'right heptamer' , 'right nonamer', 'start of gene', 'end of gene', 'gene sequence']]
    else:
         detected_gene_info = [['reference contig', 'strand', 'heptamer index' , 'nonamer index' ,\
                               'heptamer' , 'nonamer', 'start of gene', 'end of gene','best aligned human gene', \
                                'alignment direction', 'alignment PI', 'longest common k-mer','gene sequence']]
    detected_gene_info.extend(predictions)
    with open(filepath, "w", newline="") as f:
        writer =csv.writer(f , delimiter = '\t')
        writer.writerows(detected_gene_info)


print("Finding immunoglobulin genes for locus " + LOCUS + '...')

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
for gene in GENE_TYPES_TOFIND:
#    print(gene, input_rss_info[gene], s_fragments[gene], s_fragment_alignment[gene])
    if gene == D:
        print_predicted_genes('{}/genes_{}.tsv'.format(OUTPUT_PATH, D) , D, extract_genes(input_seq_dict, D, input_rss_info[D], None, None))
    else:
        predictions = extract_genes(input_seq_dict, gene, input_rss_info[gene], s_fragments[gene], s_fragment_alignment[gene])
        print_predicted_genes('{}/genes_{}.tsv'.format(OUTPUT_PATH, gene) , gene, predictions)

print("Please see {}/ for gene predictions".format(OUTPUT_PATH))  

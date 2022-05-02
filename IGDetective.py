import os,getopt
import sys
import csv
import itertools
import pickle

from Bio.Seq import Seq
from Bio import SeqIO

from multiprocessing import Pool

#INITIALIZE VARIABLES
V = 'V'
D = 'D'
J = 'J'
DR = 'D_right'
DL = 'D_left'
GENE_TYPES = [V,D,J]
GENE_TYPES_TOFIND = [V,D,J]
SPACER_LENGTH = {V:23 , DL:12 , DR:12 , J:23}
GENE_LENGTH = {V: 350 , J: 70}
FWD = '+'
REV = '-'




#READ DATAFILES############################################
try:
    with open('datafiles/motifs', 'rb') as f:   
        VALID_MOTIFS = pickle.load(f)
except:
    print("Error: could not find the input data files. Please make sure the IGDetective.py file and datafiles folder are in the same directory")

#PARSE COMMAND LINE ARGUMENTS##############################
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


#READ INPUT FASTA FILE#####################################
input_seq_dict= {rec.id : rec.seq for rec in SeqIO.parse(INPUT_PATH, "fasta")}

#DEFINE RSS FINDING METHODS################################


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

#FIND RSS IN INPUT FASTA################################
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


##ALIGN V-FRAGMENTS TO HUMAN V GENES##############################
import os
import sys,getopt
import pandas as pd
import math
from statistics import *
import numpy as np
from numpy import array, all
import random
import csv
import copy
import re
from multiprocessing import Pool

from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align

#filenames and constants
genetypes = ['V','D','J']
signaltypes = ['V','D_left','D_right','J']
n2a= {0:'A' , 1:'C' , 2:'G' , 3:'T'}
a2n= {'A':0 , 'C':1 , 'G':2 , 'T':3}
VJ_spacer = 23
D_spacer = 12
max_gene_length = {'V' : 350 , 'D' : 150 , 'J':70}
max_common_kmer = {'V': 15 , 'J' : 11}
human_genes = {'V' : None , 'J' : None}

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner_scoring = {'maxk' :[1,-10000,-10000,-10000],
                   'fin' :[1,0,-1,-0.5]
                  }
try:
    with open("datafiles/canonical_rss.txt","r") as f:
        canonical_rss = eval(f.read())
        
    with open("datafiles/allsignals.txt", "r") as f:
        signals = eval(f.read())
        
    with open("datafiles/top_allprofile.txt",'r') as f:
        allprofile = eval(f.read())

    with open("datafiles/top_allconsensus.txt",'r') as f:
        allconsensus = eval(f.read())
    
    with open("datafiles/Lmin.txt",'r') as f:
        L_all_vals = eval(f.read())
    
    for k in human_genes:
        human_genes[k]= {str(rec.id) : str(rec.seq) for rec in SeqIO.parse('datafiles/human_{}.fasta'.format(k), "fasta")}
            
except:
    print("Error: could not find the input data files. Please make sure the IGDetective.py file and datafiles folder are in the same directory")
    
#parse args
argumentList = sys.argv[1:]
options = "h:i:o:r:m:"
long_options = ["help=","input_file=", "output_directory=","reference_matrix=", "multi_process="]
force_output = True
received_input = False
num_threads = 1
ref_spec = 'combined'
try:
    arguments, values = getopt.getopt(argumentList, options, long_options)
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-h", "--help"):
            print ("Diplaying Help")
            print("Flags and their usage :")
            print("-h , --help : Get this message")
            print("-i, --input_file : provide a fasta file for gene detection")
            print("-o, --output_directory : (optional) provide an output directory for generated results. Default location is in the parent directory of the input file")
            print("-r, --reference_matrix : (optional) provide alternate reference. Default is combined matrix")
            print("-m, --multi_process : (optional) provide number of parallel processing units if available. Default is 1")
            
        elif currentArgument in ("-i", "--input_file"):
            input_path = str(currentValue)
            received_input = True
        
        elif currentArgument in ("-o", "--output_directory"):
            output_path = str(currentValue)
            force_output = False
        
        elif currentArgument in ("-r", "--reference_matrix"):
            ref_spec = str(currentValue)

        elif currentArgument in ("-m", "--multi_process"):
            num_threads = int(currentValue)
    if not received_input:
    	raise NameError('no input file was given')
    			
except getopt.error as err:
    print (str(err))
    sys.exit(0)

if force_output == True:
    output_path = ".".join(input_path.split('.')[:-1])
if not os.path.exists(output_path):
    os.makedirs(output_path)
    
input_seq_dict= {str(rec.id) : str(rec.seq).upper() for rec in SeqIO.parse(input_path, "fasta")}
contig_list = list(input_seq_dict.keys())
for contig in contig_list:
    input_seq_dict[contig+"_reverse_complement"] = str(Seq(input_seq_dict[contig]).reverse_complement()).upper()
 
#############################################################################methods to find candidate RSSs##########################################################################################

#prob(string | profilematrix)
def prob_string_given_profile(sequence , prof):
    a_sequence = [a2n[x] for x in sequence]
    prob = 1
    for idx,s in enumerate(a_sequence):
        prob = prob* prof[s,idx]
    return prob

#Find L values of all k-mers in input

def find_L_value(kmer):
    if(bool(re.match('^[ACGTacgt]+$', kmer))):
        L = prob_string_given_profile(kmer,find_all_L_values.k_profile)/find_all_L_values.k_den
    else:
        L = 0
    
    return L
    
#Find L values of all k-mers in input contig
def find_all_L_values(input_contig,ref_matrix):
    probability_dict = {g: {k: None for k in ('7','9')} for g in signaltypes}
    for g in signaltypes:
        for k in ['7','9']:
                L_all = []
                k_all = []
                find_all_L_values.k_profile = allprofile[ref_matrix][g][k]
                find_all_L_values.k_den = prob_string_given_profile(allconsensus[ref_matrix][g][k] ,find_all_L_values.k_profile)
                S = input_seq_dict[input_contig]
                numk = int(k)
                
                for i in range(0,len(S)-numk+1):
                    k_all.append(S[i:i+numk])
                
                p = Pool(num_threads)
                L_all = p.map(find_L_value,k_all)

                probability_dict[g][k] = L_all
    return probability_dict


#first pass of IGDetective
#finds all 7-mer and 9-mers over L min threshold seperately
def find_denovo_kmer(p_dict_vector,L):
    return [idx for idx,x in enumerate(p_dict_vector) if x>L]

def find_pass1_counts(probability_dict,L_vals):
    pass1_mapped_indices = {g:{k: [] for k in ('7','9')} for g in signaltypes}
    L_experiment = {g:{k: [] for k in ('7','9')} for g in signaltypes}
    for g in signaltypes:
        L_experiment[g]['7'] = L_vals['7'][g]
        L_experiment[g]['9'] = L_vals['9'][g]
        for l7 in L_vals['7'][g]:
            pass1_mapped_indices[g]['7'].append(find_denovo_kmer(probability_dict[g]['7'],l7))
        for l9 in L_vals['9'][g]:
            pass1_mapped_indices[g]['9'].append(find_denovo_kmer(probability_dict[g]['9'],l9))

    return pass1_mapped_indices , L_experiment

#second pass of IGDetective
#finds all candidate 7-mer and 9-mer pairs forming candidate RSS.
def find_pass2_counts(pass1_mapped_indices,L_experiment):
    pass2_mapped_indices = {g: [] for g in signaltypes}
    for g in signaltypes:
        if g == 'V':
            spacer = 7 + VJ_spacer 
        elif g == 'J':
            spacer = -VJ_spacer - 9
        elif g == 'D_left':
            spacer = -D_spacer - 9
        elif g == 'D_right':
            spacer = 7 +D_spacer

        for l_7 in range(0, len(L_experiment[g]['7'])):
            for l_9 in range(0, len(L_experiment[g]['9'])):
                pass2_indices = []
                seven_indices = set(pass1_mapped_indices[g]['7'][l_7])
                nine_indices = set(pass1_mapped_indices[g]['9'][l_9])
                for si in seven_indices:
                    if (si+spacer in nine_indices):
                        pass2_indices.append((si,si+spacer))
                    elif (si+spacer-1 in nine_indices):
                        pass2_indices.append((si,si+spacer-1))
                    elif (si+spacer+1 in nine_indices):
                        pass2_indices.append((si,si+spacer+1))
                        
                pass2_mapped_indices[g].append(pass2_indices)
                
    return pass2_mapped_indices

#third pass, where we ensure D_left and D_right are separated by a length of |gene| ~ 150 BP
def find_pass3_counts(pass2_mapped_indices,L_experiment):
    pass3_mapped_indices = copy.deepcopy(pass2_mapped_indices)
    pass3_mapped_indices['D'] = []

    for l_7_l in range(0,len(L_experiment['D_left']['7'])):
        for l_9_l in range(0,len(L_experiment['D_left']['9'])):
            l_l  = l_7_l * len(L_experiment['D_left']['7']) + l_9_l
            for l_7_r in range(0,len(L_experiment['D_right']['7'])):
                for l_9_r in range(0,len(L_experiment['D_right']['9'])):
                    l_r  = l_7_r * len(L_experiment['D_right']['7']) + l_9_r

                    Dlist = []
                    for x_ in pass3_mapped_indices['D_left'][l_l]:
                        for y_ in pass3_mapped_indices['D_right'][l_r]:
                            if y_[0] - x_[0] - 7 <= max_gene_length['D'] and y_[0] > x_[0]:
                                Dlist.append((x_,y_))

                    pass3_mapped_indices['D'].append(Dlist)
                    
    return pass3_mapped_indices


####################################################################methods to find candidate genes##################################################################

#define alignment score based on scheme
def set_aligner(scheme):
    aligner.match_score = aligner_scoring[scheme][0]
    aligner.mismatch_score = aligner_scoring[scheme][1]
    aligner.open_gap_score = aligner_scoring[scheme][2]
    aligner.extend_gap_score = aligner_scoring[scheme][3]

#count length of alignment 
def count_alignment_len(al,gene_al):
    start,end = FindAlignmentRange(al,gene_al)
    return 1+end-start

#find start and end of alignment
def FindAlignmentRange(al,gene_al):
    if '-' in al:
        start,end = min(al.find('|'),al.find('-')) , max(al.rfind('|'),al[1].rfind('-'))
    else:
        start,end = al.find('|') , al.rfind('|')
        
    residual = len(gene_al[:start]) - gene_al[:start].count(' ')
    start -= residual
    return start,end

#get numnber of matches in an alignment
def GetNumMatches(alignment):
    splits = str(alignment).split('\n')
    query_alignment = splits[0].upper()
    gene_alignment = splits[2].upper()
    coded_alignment = splits[1].upper()
    start,end = FindAlignmentRange(coded_alignment,gene_alignment)
    matches = (coded_alignment[start : end+1]).count('|')
    
    return matches

#align 2 strings
def ComputeAlignment(aligner, locus_fragment, gene_seq):
    query = Seq(locus_fragment)
    query_rc = query.reverse_complement()
    alignment = aligner.align(query, gene_seq)[0]
    alignment_rc = aligner.align(query_rc, gene_seq)[0]
    if GetNumMatches(alignment) > GetNumMatches(alignment_rc):
        return alignment, '+', GetNumMatches(alignment)
    return alignment_rc, '-', GetNumMatches(alignment_rc)

#compute best alignment between any string and all human genes
def ComputeBestAlignment(aligner, locus_fragment, gene_dict):
    best_gene = ''
    best_PI = 0
    best_alignment = ''
    best_strand = ''
    for gene in gene_dict:
        alignment, strand, num_matches = ComputeAlignment(aligner, locus_fragment.upper(), gene_dict[gene].upper())
        splits = str(alignment).split('\n')
        alig_len = count_alignment_len(splits[1],splits[2])
        PI = num_matches *100 / alig_len
        if PI > best_PI:
            best_alignment = alignment
            best_strand = strand
            best_gene = gene
            best_PI = PI
    return best_alignment, best_strand, best_gene

#method to get alignment and all related statistics between a set of candidate RSS and all human genes
def get_alignments(contig,gene_type,alig_scheme,pass3):
    #values to return
    genes = [] 
    perc_identities = []
    candidate_7 = []
    candidate_9 = []
    local_alig_len = []
    alignments = []
    gene_strands = []
    
    #logic starts here
    if gene_type in ['V','J']:
        spec_genes = human_genes[gene_type]
        
        for i in range(0,len(pass3)):
            heptamer_pos = pass3[i][0]
            nonamer_pos = pass3[i][1]
            
            if heptamer_pos <= max_gene_length[gene_type] or (len(input_seq_dict[contig]) - heptamer_pos)<= max_gene_length[gene_type]+7 :
                local_alig_len.append(-1)
                perc_identities.append(-1)
                alignments.append("maxlen_elim")
                gene_strands.append("maxlen_elim")
                genes.append("maxlen_elim")
                candidate_7.append('maxlen_elim')
                candidate_9.append('maxlen_elim')
                continue
        
    
            heptamer = input_seq_dict[contig][heptamer_pos : heptamer_pos+7]
            nonamer = input_seq_dict[contig][nonamer_pos : nonamer_pos + 9]
            
            if gene_type == 'V':
                gene_superstr = input_seq_dict[contig][heptamer_pos - max_gene_length[gene_type] : heptamer_pos]
            elif gene_type == 'J':
                gene_superstr = input_seq_dict[contig][heptamer_pos+7:heptamer_pos+7+max_gene_length[gene_type]]
                
            if 'N' in gene_superstr:
                local_alig_len.append(-1)
                perc_identities.append(-1)
                #perc_identities_2.append(-1)
                alignments.append("N_elim")
                gene_strands.append("N_elim")
                genes.append("N_elim")
                candidate_7.append('N_elim')
                candidate_9.append('N_elim')
                continue
            
            set_aligner(alig_scheme)
            alignment, gene_strand, gene  = ComputeBestAlignment(aligner, gene_superstr, spec_genes)
            num_matches = float(GetNumMatches(alignment))
            splits = str(alignment).split('\n')
            alig_len = count_alignment_len(splits[1] , splits[2])
            perc_identity = num_matches* 100 / alig_len
            
            
            candidate_7.append(heptamer_pos)
            candidate_9.append(nonamer_pos)
            local_alig_len.append(alig_len)
            perc_identities.append(perc_identity)
            alignments.append(str(alignment).split('\n')[:-1])
            gene_strands.append(gene_strand)
            genes.append(gene)

            
    return(perc_identities,local_alig_len,alignments,gene_strands,genes,candidate_7,candidate_9)
    
# method to predict genes from alignments
def predict_genes(gene_type, contig, args):
    if contig.endswith('_reverse_complement'):
        strand = '-'
        contig_print = contig[:-len('_reverse_complement')]
    else:
        strand ='+'
        contig_print = contig
        
    perc_identities = args[0]
    local_alig_len = args[1]
    alignments = args[2]
    gene_strands = args[3]
    genes = args[4]
    candidate_7 = args[5]
    candidate_9 = args[6]
    maxks = args[7]
    
    passlist = []
    
    for i in range(0,len(args[0])):
        gene_pass = False

        if perc_identities[i]>=70:
            gene_pass = True
        elif perc_identities[i] >=60 and maxks[i]>=max_common_kmer[gene_type]:
            gene_pass = True
            
        if gene_pass == True:
            start,end = FindAlignmentRange(alignments[i][1],alignments[i][2] )
            end = end - alignments[i][0].count('-') +1
            heptamer = input_seq_dict[contig][candidate_7[i]:candidate_7[i]+7]
            nonamer = input_seq_dict[contig][candidate_9[i]:candidate_9[i]+9]
            
            if gene_type == 'V':
                detected_start = candidate_7[i]-max_gene_length['V']+start
                detected_end = candidate_7[i] -1
            elif gene_type == 'J':
                detected_start = candidate_7[i]+7
                detected_end = candidate_7[i]+7+end
            detected = input_seq_dict[contig][detected_start:detected_end].upper()
            
            if strand == '+':
                print_start = detected_start
                print_end = detected_end
                print_heptamer_pos = candidate_7[i]
                print_nonamer_pos = candidate_9[i]
            elif strand == '-':
                print_start = len(input_seq_dict[contig]) - 1 - detected_end
                print_end = len(input_seq_dict[contig]) - 1 - detected_start
                print_heptamer_pos =len(input_seq_dict[contig]) - 1 - (candidate_7[i]+6)
                print_nonamer_pos = len(input_seq_dict[contig]) - 1 - (candidate_9[i]+8)
                
            
            
            gene_print = [contig_print,strand,print_heptamer_pos,print_nonamer_pos,heptamer,nonamer,\
                         genes[i],gene_strands[i],perc_identities[i],maxks[i],print_start,print_end,\
                         detected]
            passlist.append(gene_print)
        
                
    return passlist    


#method to identify D genes
def identify_Dgenes(gene_type,contig,pass3):
    if contig.endswith('_reverse_complement'):
        strand = '-'
        contig_print = contig[:-len('_reverse_complement')]
    else:
        strand ='+'
        contig_print = contig

    passlist = []
    for i in pass3:
        left_heptamer = input_seq_dict[contig][i[0][0]:i[0][0]+7]
        left_nonamer = input_seq_dict[contig][i[0][1]:i[0][1]+9]
        right_heptamer = input_seq_dict[contig][i[1][0]:i[1][0]+7]
        right_nonamer = input_seq_dict[contig][i[1][1]:i[1][1]+9]
        detected_start = i[0][0]+7
        detected_end = i[1][0]
        detected = input_seq_dict[contig][detected_start:detected_end]
        
        if strand == '+':
            l7print,l9print,r7print,r9print = i[0][0],i[0][1],i[1][0],i[1][1]
            print_start = detected_start
            print_end = detected_end
        elif strand =='-':
            l7print,l9print,r7print,r9print = len(input_seq_dict[contig]) - 1 - (i[0][0]+6),\
                                                len(input_seq_dict[contig]) - 1 - (i[0][1]+8),\
                                                len(input_seq_dict[contig]) - 1 - (i[1][0]+6),\
                                                len(input_seq_dict[contig]) - 1 - (i[1][1]+8)
            print_start = len(input_seq_dict[contig]) - 1 - detected_end
            print_end = len(input_seq_dict[contig]) - 1 - detected_start
            
        
        gene_print = [contig_print,strand,l7print,l9print,left_heptamer,left_nonamer,\
                      r7print,r9print,right_heptamer,right_nonamer,print_start,print_end,detected] 
        passlist.append(gene_print)
        
    return passlist  
####################################################################################Run IGDetective###################################################################################
#find Candidate RSS 
pass3_contigs = {}
print("Finding candidate RSS...",end =" ")
for k in input_seq_dict:
    probability_dict = find_all_L_values(k,ref_spec)
    L_vals = L_all_vals[ref_spec]
    (pass1_mapped_indices , L_experiment) = find_pass1_counts(probability_dict,L_vals)
    pass2_mapped_indices = find_pass2_counts(pass1_mapped_indices,L_experiment)
    pass3_mapped_indices = find_pass3_counts(pass2_mapped_indices,L_experiment)
    pass3_contigs[k] = pass3_mapped_indices
print("Done\n")

#find V and J genes
JV_header = ['contig','reference strand', 'heptamer index', 'nonamer index', 'heptamer', 'nonamer',\
            'best aligned human gene','alignment direction','alignment PI','longest common k-mer','start of gene',\
            'end of gene','gene sequence']

for g_type in ['V','J']:
    print("Identifying {} genes...".format(g_type),end =" ")
    all_print = []
    f= open("{}/{}_predicted_genes.csv".format(output_path,g_type), "w")
    writer = csv.writer(f)
    writer.writerow(JV_header)
    for k in input_seq_dict:
        pass3 = pass3_contigs[k][g_type][0]
        alig_result = get_alignments(k,g_type,'maxk',pass3)
        maxk_len = alig_result[1]
        alig_result = get_alignments(k,g_type,'fin',pass3)
        args = list(alig_result)
        args.append(maxk_len)
        contig_print = predict_genes(g_type, k, tuple(args))
        all_print.extend(contig_print) 
    writer.writerows(all_print)
    f.close()
    print("Done\n")

#find D genes
D_header = ['contig','reference strand', 'left heptamer index', 'left nonamer index', 'left heptamer', 'left nonamer',\
        'right heptamer index', 'right nonamer index', 'right heptamer', 'right nonamer','start of gene',\
        'end of gene','gene sequence']

print("Identifying D genes...",end =" ")
all_print = []
g_type = 'D'
f= open("{}/{}_predicted_genes.csv".format(output_path,g_type), "w")
writer = csv.writer(f)
writer.writerow(D_header)
for k in input_seq_dict:
    pass3 = pass3_contigs[k][g_type][0]
    contig_print = identify_Dgenes(g_type,k,pass3)
    all_print.extend(contig_print) 
writer.writerows(all_print)
f.close()
print("Done\n")


import os,getopt
import sys
import csv
import itertools
import pickle
import math
from statistics import median,mean
import copy



from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align

import numpy as np
from numpy import array, all
import networkx as nx
from scipy.cluster.hierarchy import single, fcluster
from scipy.spatial.distance import pdist

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
GENE_TYPES_TOFIND = [V]
SPACER_LENGTH = {V:23 , DL:12 , DR:12 , J:23}
GENE_LENGTH = {V: 350 , J: 70, D:150}
ALIGNMENT_EXTENSION = {V:REV , J:FWD, D:None}
PI_CUTOFF = {'strict' : {V: 70, J: 70} , 'relax': {V: 60 ,J: 65}}
MAXK_CUTOFF = {V: 15 , J: 11}
CLUSTER_SPAN_LIMIT = 1000000
COLOCALIZATION_SPAN = 500000
COMPONENT_SIZE = {'GIANT': 500, 'LARGE' : 3}
ANNOTATION_THRESHOLD = {'human': 80, 'species': 90}
GRANDPARENT_FOLDER = "/".join(os.path.abspath(os.getcwd()).split('/')[:-1])


ALIGNER = Align.PairwiseAligner()
ALIGNER.mode = 'local'


#PARSE COMMAND LINE ARGUMENTS
argumentList = sys.argv[1:]
options = "hi:a:o:m:"
long_options = ["help","input_file=", "annotated_genes=" ,"output_directory=", "multi_process="]
force_output = True
received_input = False
help_flag = False
NUM_THREADS = 1
BLIND_ANNOTATED_FILE = None
try:
    arguments, values = getopt.getopt(argumentList, options, long_options)
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-h", "--help"):
            print ("Diplaying Help")
            print("Flags and their usage :")
            print("-h , --help : Get this message")
            print("-i, --input_file : provide a fasta file for gene detection")
            print("-a, --annotated_genes : provide a directory containing annotated_genes")
            print("-o, --output_directory : (optional) provide an output directory for generated results. Default location is in the parent directory of the input file")
            print("-m, --multi_process : (optional) provide number of parallel processing units if available. Default is 1")
            help_flag = True
            sys.exit(0)
            
        elif currentArgument in ("-i", "--input_file"):
            INPUT_PATH = str(currentValue)
            received_input = True

        elif currentArgument in ("-a", "--annotated_genes"):
            BLIND_ANNOTATED_FILE = str(currentValue)
        
        elif currentArgument in ("-o", "--output_directory"):
            OUTPUT_PATH = str(currentValue)
            force_output = False

        elif currentArgument in ("-m", "--multi_process"):
            NUM_THREADS = int(currentValue)


    if not received_input and not help_flag:
        raise NameError('no input file was given')
                
except getopt.error as err:
    print (str(err))
    sys.exit(0)
if force_output == True:
    OUTPUT_PATH = ".".join(INPUT_PATH.split('.')[:-1])
if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)
#Read input Fasta file
input_seq_dict= {rec.id : rec.seq.upper() for rec in SeqIO.parse(INPUT_PATH, "fasta")}

#READ DATAFILES
try:
    with open('{}/datafiles/motifs'.format(GRANDPARENT_FOLDER), 'rb') as f:   
        VALID_MOTIFS = pickle.load(f)
except:
    print("Error: could not find the input motif files. Please make sure the run_blind_igdetective.py file and datafiles folder are in the same directory and that motifs is present in datafiles")
    sys.exit(0)

blind_annotated_genes = {spec : {} for spec in ('human','species')}
try:
    for g in GENE_TYPES_TOFIND:
        filename = '{}/datafiles/blind_ig/annotated_human_{}.fa'.format(GRANDPARENT_FOLDER,g)
        blind_annotated_genes['human'][g] = {rec.id : rec.seq.upper() for rec in SeqIO.parse(filename, 'fasta')}

        if BLIND_ANNOTATED_FILE:
            filename = BLIND_ANNOTATED_FILE
        blind_annotated_genes['species'][g] = {rec.id : rec.seq.upper() for rec in SeqIO.parse(filename, 'fasta')}

except:
    print("Error: could not find the input annotated gene files. Please make sure the run_blind_igdetective.py file and datafiles folder are in the same directory and that the annotated_human_{}.fa is present in datafiles".format(g))
    sys.exit(0)
    
#create signal types from gene types
SIGNAL_TYPES = []
if V in GENE_TYPES_TOFIND:
    SIGNAL_TYPES.append(V)
if J in GENE_TYPES_TOFIND:
    SIGNAL_TYPES.append(J)
if D in GENE_TYPES_TOFIND:
    SIGNAL_TYPES.extend([DL,DR])

#temp delete some contigs
for key in list(input_seq_dict.keys()):
    if key!= 'NC_000014.9':
        del(input_seq_dict[key])

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
        if strand == FWD:
            first_set = [x for x in heptamer_idx if x>= GENE_LENGTH[sig_type]]
        elif strand == REV:
            first_set = [x for x in heptamer_idx if seq_length - x >= GENE_LENGTH[sig_type]+7]
        second_set = nonamer_idx
    elif sig_type == J or sig_type == DL:
        k_first = 9
        if strand == FWD:
            second_set = [x for x in heptamer_idx if x>= GENE_LENGTH[sig_type]]
        elif strand == REV:
            second_set = [x for x in heptamer_idx if seq_length - x >= GENE_LENGTH[sig_type]+7]
        first_set = nonamer_idx
    
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
        
# S FRAGMENT EXTRACTION
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
                
#PAIRWISE ALIGNMENT
#define alignment score based on scheme
def set_aligner(scheme):
    aligner_scoring = {'MAXK' :[1,-10000,-10000,-10000],
                   'AFFINE' :[1,0,-1,-0.5]}
    ALIGNER.match_score = aligner_scoring[scheme][0]
    ALIGNER.mismatch_score = aligner_scoring[scheme][1]
    ALIGNER.open_gap_score = aligner_scoring[scheme][2]
    ALIGNER.extend_gap_score = aligner_scoring[scheme][3]

#count length of alignment
def count_alignment_len(al, gene_al, extend_alignment = None):
    start,end = find_alignment_range(al, gene_al,  extend_alignment)
    return 1+end-start

#find start and end of alignment
def find_alignment_range(al, gene_al, extend_alignment = None):
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
def get_num_matches(alignment, extend_alignment = None):
    splits = str(alignment).split('\n')
    seq_A_alignemnt = splits[0].upper()
    seq_B_alignment = splits[2].upper()
    coded_alignment = splits[1].upper()
    start,end = find_alignment_range(coded_alignment, seq_B_alignment,  extend_alignment)
    matches = (coded_alignment[start : end+1]).count('|')

    return matches

#align 2 strings
def compute_alignment(seq_A, seq_B, extend_alignment = None):
    query = seq_A
    query_rc = query.reverse_complement()
    alignment = ALIGNER.align(query, seq_B)[0]
    alignment_rc = ALIGNER.align(query_rc, seq_B)[0]
    fwd_matches  = get_num_matches(alignment, extend_alignment)
    rev_matches = get_num_matches(alignment_rc, extend_alignment)
    if  fwd_matches > rev_matches:
        return alignment, '+', fwd_matches
    else:
        return alignment_rc, '-', rev_matches

def get_PI_matrix(fragments, scoring_scheme, fragments_other = None):
    set_aligner(scoring_scheme)
    parallel_alignments = []
    pis = []
    
    if fragments_other is None:
        pi_mat = np.zeros((len(fragments), len(fragments)))
        for i in range(0,len(fragments)):
            for j in range(0,len(fragments)):
                if j>i:
                    seq_A = fragments[i]
                    seq_B = fragments[j]
                    parallel_alignments.append((seq_A,seq_B))
    else:
        pi_mat = np.zeros((len(fragments), len(fragments_other)))
        for i in range(0,len(fragments)):
            for j in range(0,len(fragments_other)):
                seq_A = fragments[i]
                seq_B = fragments_other[j]
                parallel_alignments.append((seq_A,seq_B))
                
    p = Pool(NUM_THREADS)
    alignment_results = p.starmap(compute_alignment , parallel_alignments)
    for a in alignment_results:
        alignment = a[0]
        strand = a[1]
        num_matches = a[2]
        splits = str(alignment).split('\n')
        alig_len = count_alignment_len(splits[1],splits[2])
        perc_identity = num_matches* 100 / alig_len
        pis.append(perc_identity)
    
    if fragments_other is None:
        for i in range(0,len(fragments)):
            for j in range(0,len(fragments)):
                if j>i:
                    pi_mat[i][j] = pis.pop(0)
    else:
        for i in range(0,len(fragments)):
            for j in range(0,len(fragments_other)):
                pi_mat[i][j] = pis.pop(0)
        
    return pi_mat

#S GRAPH ANALYSIS
def make_S_graph(single_vertex_list,PI_matrix):
    num_rows, num_cols = PI_matrix.shape
    edge_list = []
    
    vertex_list = {n:single_vertex_list[n] for n in range(0,len(single_vertex_list))}
    for i in range(0,num_rows):
        for j in range(0,num_cols):
            if j>i and (PI_matrix[i][j]>=70):
                edge_list.append((i,j , PI_matrix[i][j]))
    G = nx.Graph()
    G.add_nodes_from(vertex_list)
    for e in edge_list:
        G.add_edge(e[0],e[1],weight = e[2])
    
    return G

#helper function to find extended coding length for V genes
#new_start for REV implies the end of gene as per usal notation - Caution!
def find_extended_coding_length(G, node, frag_len):
    strand = G.nodes[node]['strand']
    index = G.nodes[node]['index']
    frame = G.nodes[node]['aa_frame']
    contig = G.nodes[node]['contig']
    coding_length = G.nodes[node]['coding_length']
    added_codons = 0
    
    if strand == FWD:
        coding_start_index = index + frame
        while(True):
            codon = input_seq_dict[contig][coding_start_index-3:coding_start_index]
            if codon.translate() == '*':
                break
            coding_start_index -=3
            added_codons+=1
        real_start_coordinate =  index + frame - (added_codons*3)
        new_coding_leng = coding_length + (added_codons*3)
    elif strand == REV:
        coding_start_index = index + frag_len - frame
        while(True):
            codon = input_seq_dict[contig][coding_start_index:coding_start_index+3]
            if (codon.reverse_complement()).translate() == '*':
                break
            coding_start_index +=3
            added_codons+=1
        real_start_coordinate = coding_start_index
        new_coding_leng = coding_length + (added_codons*3)
    return real_start_coordinate , new_coding_leng
    
#add node meta data
def add_node_data(g, fragments, G, contig, rss_info):
    for node in list(G.nodes):
        G.nodes[node]['contig'] = contig
        incident_edges = list(G.edges([node]))
        incident_weights = [G.get_edge_data(x[0],x[1])['weight'] for x in incident_edges]
        G.nodes[node]['weight'] = max(incident_weights)
        if node < len(fragments[FWD][contig]):
            fragment_index = node
            strand = FWD
            G.nodes[node]['index'] = rss_info[g][strand][contig][fragment_index][0] - GENE_LENGTH[g]
        else:
            fragment_index = node - len(fragments[FWD][contig])
            strand = REV
            G.nodes[node]['index']  = rss_info[g][strand][contig][fragment_index][0] +7
        
        G.nodes[node]['dna_seq'] = fragments[strand][contig][fragment_index]
        G.nodes[node]['strand'] = strand  
        aa1 = fragments[strand][contig][fragment_index].translate()
        aa2 = fragments[strand][contig][fragment_index][1:].translate()
        aa3 = fragments[strand][contig][fragment_index][2:].translate()
        aas = [len(aa1.split('*')[-1]) , len(aa2.split('*')[-1]) , len(aa3.split('*')[-1])]
        coding_leng = 3 * max(aas)
        aa = [aa1,aa2,aa3][aas.index(max(aas))]
        frame = aas.index(max(aas))
        G.nodes[node]['coding_length'] = coding_leng
        G.nodes[node]['aa_seq'] = aa
        G.nodes[node]['aa_list'] = [aa1,aa2,aa3]
        G.nodes[node]['aa_frame'] = frame        
        if coding_leng == math.floor(len(G.nodes[node]['dna_seq']) / 3 )*3:
            real_start_coordinate, new_coding_length = find_extended_coding_length(G, node, len(G.nodes[node]['dna_seq']))
            G.nodes[node]['coding_length'] = new_coding_length
            
    return G

def find_connected_components(graph):
    return nx.connected_components(graph)


#CLUSTER CREATION AND ANALYSIS
def get_group_colocalized_nodes(G,node_set,node):
    colocalized = []
    for other_node in node_set:
        if other_node != node and \
        G.nodes[node]['contig'] == G.nodes[other_node]['contig'] and \
        abs(G.nodes[node]['index'] - G.nodes[other_node]['index']) <= COLOCALIZATION_SPAN :
            colocalized.append(other_node)
    return colocalized

def get_group_center_node(G,group):
    power = {x : len(get_group_colocalized_nodes(G,group,x)) for x in group}
    center = max(power, key=power.get)
    return center, power[center]

def get_clumps_from_component(G,component_id, components):
    clumps = {}
    numclmp = 0
    component_nodes = copy.deepcopy(components[component_id])
    while component_nodes:
        center, power = get_group_center_node(G,component_nodes)
        center_colocals = get_group_colocalized_nodes(G,component_nodes,center)
        center_colocals = center_colocals + [center]
        
        #remove isolated nodes
        k = G.subgraph(center_colocals)
        non_isolates = set(center_colocals).difference(set(nx.isolates(k)))
        if non_isolates:
            clumps[numclmp] = sorted(list(non_isolates))    
            numclmp+=1
        for clmp in center_colocals:
            component_nodes.remove(clmp)
            
    return clumps

def get_group_span(G,group):
    contig_set = set([G.nodes[x]['contig'] for x in group])
    if len(contig_set)>1:
        return -1
    indexes = [G.nodes[x]['index'] for x in group]
    span = max(indexes) - min(indexes)
    return span

def get_clusters_from_clumps(G , contig_clumps, contig_connected_components):
    selected_clumps = []
    #remove short clumps under 1k NT
    for ccid in contig_connected_components:
        for cmpid in contig_clumps[ccid]:
            span = get_group_span(G,contig_clumps[ccid][cmpid])
            if span > 1000:
                selected_clumps.append((ccid,cmpid))
    
    contig_clump_centers_index = {ccid : {cmpid : None for cmpid in contig_clumps[ccid]} for ccid in contig_connected_components}
    centers_index = []
    for (ccid,cmpid) in selected_clumps:
        clump = contig_clumps[ccid][cmpid]
        contig_clump_centers_index[ccid][cmpid] =  G.nodes[get_group_center_node(G, clump)[0]]['index']
        centers_index.append([G.nodes[get_group_center_node(G, clump)[0]]['index']])
        
    clusters_vertices = {}
    clusters_clumps = {}
    if len(centers_index) > 1:
        vertex_set = []
        cluster_set = []
        dist_p = pdist(centers_index)
        Z = single(dist_p)
        clustering = fcluster(Z, CLUSTER_SPAN_LIMIT, criterion='distance')
        for cluster_id in set(clustering):
            combined = []
            clusters_vertices[cluster_id-1]= []
            clusters_clumps[cluster_id-1] = []
            for i in range(0,len(clustering)):
                if clustering[i] == cluster_id:
                    (ccid , cmpid) = selected_clumps[i]
                    combined.extend(contig_clumps[ccid][cmpid])
                    clusters_clumps[cluster_id - 1].append((ccid , cmpid))
            clusters_vertices[cluster_id - 1] = combined
  
    else:
        ccid = selected_clumps[0]
        cmpid = selected_clumps[1]
        clusters_vertices[0] = contig_clumps[ccid][cmpid]
        clusters_clumps[0] = (ccid,cmpid)
    
    return clusters_vertices,clusters_clumps

def get_group_CL(G, group):
    return median([G.nodes[node]['coding_length'] for node in group])

def get_group_density(G, group):
    k = G.subgraph(group)
    num_edges = len(k.edges(group))
    n = len(group)
    max_edges = n*(n-1)/2
    if max_edges == 0:
        return -1 
    group_density =  num_edges*100 / max_edges
    return group_density

def get_group_PI(G, group):
    k = G.subgraph(group)
    group_PI = []
    for node in group:
        incident_edges = list(k.edges([node]))
        if incident_edges:
            incident_weights = [G.get_edge_data(x[0],x[1])['weight'] for x in incident_edges]
            group_PI.append(max(incident_weights))
    if group_PI:
        return median(group_PI)
    else:
        return -1
    
def get_group_span(G, group):
    contig_set = set([G.nodes[x]['contig'] for x in group])
    if len(contig_set)>1:
        return -1
    indexes = [G.nodes[x]['index'] for x in group]
    span = max(indexes) - min(indexes)
    return span

def get_group_conservation(G, group , conservation_type):
    conservations = [G.nodes[node][conservation_type] for node in group]
    return conservations

def get_group_AI(G, group , AI_type):
    AI_type = 'is_annotated_{}'.format(AI_type)
    return [G.nodes[node][AI_type] for node in group].count(True)/len(group)

def get_group_family(G, group):
    if get_group_AI(G, group , 'human') == 0 and get_group_AI(G, group , 'species') ==0:
        return '-'
    else:
        annotated_nodes = sorted([node for node in group if G.nodes[node]['is_annotated_species'] == True or G.nodes[node]['is_annotated_human']])
        locus = G.nodes[annotated_nodes[int(len(annotated_nodes)/2)]]['locus_human']
    return locus

#OUTPUT FORMATING
def create_summary_table(g, cluster_information):
    summary_headers = ['Cluster ID', 'Cluster size', 'No. of clumps/|Largest clump|', 'PI', 'CL', 'Density', 'Span (Mb)', 'Center contig', 'Center coordinate (Mb)', 'AI/AI80',\
                     'Min/Median conservation wrt species genes', 'Min/Median conservation wrt human genes' , 'Locus']
    summary_table = []
    clusterid = 0
    
    for contig in clusters[g]:
        for clsid in clusters[g][contig][0]:
            if cluster_information['locus'][g][contig][clsid] != '-':
                cls_size = len(clusters[g][contig][0][clsid])
                num_clumps = len(clusters[g][contig][1][clsid])
                largest_clump = max([len(clumps[g][contig][ccid][cmpid]) for ccid,cmpid in clusters[g][contig][1][clsid]])
                clump_info = '{}/{}'.format(num_clumps,largest_clump)
                pi = '{:.0f}'.format(cluster_information['PI'][g][contig][clsid])
                cl = 3*round(math.floor(cluster_information['CL'][g][contig][clsid])/3)
                density = '{:.0f}'.format(cluster_information['density'][g][contig][clsid])
                span = '{:.3f}'.format(cluster_information['span'][g][contig][clsid]/1000000)
                center_contig = contig
                center_node = get_group_center_node(S_graphs[g][contig],clusters[g][contig][0][clsid])[0]
                center_cood = '{:.3f}'.format(S_graphs[g][contig].nodes[center_node]['index']/1000000)
                ai_ai80 = '{:.2f}/{:.2f}'.format(cluster_information['AI_species'][g][contig][clsid] , cluster_information['AI_human'][g][contig][clsid])
                spec_cons = '{:.0f}/{:.0f}'.format(min(cluster_information['conservation_species'][g][contig][clsid]) , \
                                                   median(cluster_information['conservation_species'][g][contig][clsid]))
                human_cons = '{:.0f}/{:.0f}'.format(min(cluster_information['conservation_human'][g][contig][clsid]) , \
                                                   median(cluster_information['conservation_human'][g][contig][clsid]))
                locus = cluster_information['locus'][g][contig][clsid]
                summary_table.append([cls_size, clump_info, pi, cl, density, span, center_contig, center_cood, ai_ai80, spec_cons, human_cons, locus])
                
    summary_table.sort(key=lambda x: (x[-1], x[0]), reverse = True)
    for idx, row in enumerate(summary_table):
        summary_table[idx] = [idx] + summary_table[idx]
        
    return [summary_headers] + summary_table

def create_fasta_from_dict(name_seq_dict):
    fasta=""
    for x in name_seq_dict:
        fasta+= '>{}\n{}\n'.format(x,name_seq_dict[x])
    return fasta

def create_fasta_files(g, cluster_information):
    fasta_dict = {}
    gene_id = 0
    for contig in clusters[g]:
        for clsid in clusters[g][contig][0]:
            if cluster_information['locus'][g][contig][clsid] != '-':
                for node in clusters[g][contig][0][clsid]:
                    vertex = S_graphs[g][contig].nodes[node]
                    if vertex['is_annotated_species'] == True:
                        annotation_flag = 'as'
                        neighbour = vertex['annotated_neighbour_species']
                        locus = vertex['locus_species']
                        seq = vertex['annotated_aligned_gene_species']
                    elif vertex['is_annotated_human'] == True:
                        annotation_flag = 'ah'
                        neighbour = vertex['annotated_neighbour_human']
                        locus = vertex['locus_human']
                        seq = vertex['annotated_aligned_gene_human']
                    else:
                        annotation_flag = 'un'
                        neighbour = 'UNK'
                        locus = 'UNK'
                        seq = vertex['dna_seq']
                        
                    name = '{}|{}|{}|{}|prediction_{}'.format(locus, contig, annotation_flag, neighbour, gene_id)
                    gene_id+=1
                    sequence = str(seq).upper()
                    fasta_dict[name] = sequence
    
    return create_fasta_from_dict(fasta_dict)

def print_results(cluster_information):
    for g in GENE_TYPES_TOFIND:
        summary_sheet = create_summary_table(V, cluster_information)
        fasta_seqs = create_fasta_files(V, cluster_information)
        with open('{}/{}_cluster_summary.tsv'.format(OUTPUT_PATH,g), "w", newline="") as f:
            writer =csv.writer(f , delimiter = '\t')
            writer.writerows(summary_sheet)
        with open('{}/{}_gene_predictions.fasta'.format(OUTPUT_PATH,g), "w", newline="") as f:
            f.write(fasta_seqs)

#Find RSS in input fasta
print("Finding candidate RSS...",end =" ")
input_rss_info = {st : {strand : get_contigwise_rss(st,strand, input_seq_dict) for strand in (FWD , REV)} for st in SIGNAL_TYPES}
if D in GENE_TYPES_TOFIND:
    input_rss_info[D] = { strand: combine_D_RSS(input_rss_info[DL][strand] , input_rss_info[DR][strand], input_seq_dict , strand)\
                       for strand in (FWD, REV)}

for st in SIGNAL_TYPES:
    write_rss_to_file('{}/rss_{}.tsv'.format(OUTPUT_PATH, st), input_rss_info[st], input_seq_dict)
print("Done")

#Find S fragments
s_fragments = {g : {strand : get_s_fragment_from_RSS(g,strand) for strand in (FWD, REV)} for g in GENE_TYPES_TOFIND}

print("Creating alignment matrices...",end =" ")
#Calculate_PI_matrices
scoring_scheme = 'AFFINE'
similarity_matrices = {g: {contig : None for contig in set(s_fragments[g][FWD].keys()).union(set(s_fragments[g][REV].keys()))} for g in GENE_TYPES_TOFIND}
for g in GENE_TYPES_TOFIND:
    for contig in similarity_matrices[g]:
        contig_fragments = s_fragments[g][FWD][contig] + s_fragments[g][REV][contig]
        similarity_matrices[g][contig] = get_PI_matrix(contig_fragments, scoring_scheme)
print("Done")

print("Creating S-graphs...",end =" ")
#find S-graphs and connected components
S_graphs = {g: {contig : make_S_graph(s_fragments[g][FWD][contig] + s_fragments[g][REV][contig], similarity_matrices[g][contig]) for contig in similarity_matrices[g]} for g in GENE_TYPES_TOFIND}
connected_components = {g: {contig : None for contig in similarity_matrices[g]} for g in GENE_TYPES_TOFIND}

#Remove isolated nodes
for g in GENE_TYPES_TOFIND:
    for contig in  S_graphs[g]:
        isolated_nodes = list(nx.isolates(S_graphs[g][contig]))
        S_graphs[g][contig].remove_nodes_from(isolated_nodes)
        connected_components[g][contig] = list(find_connected_components(S_graphs[g][contig]))
        for cc in connected_components[g][contig]:
            if len(cc) > COMPONENT_SIZE['GIANT']:
                S_graphs[g][contig].remove_nodes_from(cc) 
for g in GENE_TYPES_TOFIND:
    for contig in similarity_matrices[g]:
        S_graphs[g][contig] = add_node_data(g, s_fragments[g],  S_graphs[g][contig] , contig, input_rss_info)
        connected_components[g][contig] = list(find_connected_components(S_graphs[g][contig]))
        connected_components[g][contig].sort(key = lambda x:len(x) , reverse = True)
        cc_dict = {idx:x for idx,x in enumerate(connected_components[g][contig])}
        connected_components[g][contig] = copy.deepcopy(cc_dict)
print("Done")

print("Annotating S-graphs...",end =" ")
#Annotating vertices
for ann_spec in ('human', 'species'):
    annotated_similarity_matrix =  {g: {contig : None for contig in S_graphs[g]} for g in GENE_TYPES_TOFIND}
    for g in GENE_TYPES_TOFIND:
        for contig in S_graphs[g]:        
            fragments_for_annotation = [S_graphs[g][contig].nodes[node]['dna_seq'] for node in S_graphs[g][contig]]
            canonical_annotated_sequences = list(blind_annotated_genes[ann_spec][g].values())
            annotated_similarity_matrix[g][contig] = get_PI_matrix(fragments_for_annotation, scoring_scheme, canonical_annotated_sequences)
            node_conservation = list(np.amax( annotated_similarity_matrix[g][contig], axis=1))
            node_closest_annotated_gene = [list(blind_annotated_genes[ann_spec][g].keys())[x] for x in list(np.argmax(annotated_similarity_matrix[g][contig], axis=1))]

            for idx,node in enumerate(S_graphs[g][contig]):
                S_graphs[g][contig].nodes[node]['conservation_{}'.format(ann_spec)] = node_conservation[idx]
                locus,neighbour = (node_closest_annotated_gene[idx]).split('|',1)
                S_graphs[g][contig].nodes[node]['annotated_neighbour_{}'.format(ann_spec)] = neighbour
                S_graphs[g][contig].nodes[node]['locus_{}'.format(ann_spec)] = locus
                if node_conservation[idx] >= ANNOTATION_THRESHOLD[ann_spec]:
                    S_graphs[g][contig].nodes[node]['is_annotated_{}'.format(ann_spec)] = True
                    if g == V:
                        a = compute_alignment(S_graphs[g][contig].nodes[node]['dna_seq'],blind_annotated_genes[ann_spec][g][node_closest_annotated_gene[idx]], ALIGNMENT_EXTENSION[g])
                        splits = str(a[0]).split('\n')
                        alig_direction = a[1]
                        start,end = find_alignment_range(splits[1].upper(), splits[2].upper(), ALIGNMENT_EXTENSION[g])
                        S_graphs[g][contig].nodes[node]['annotated_aligned_gene_{}'.format(ann_spec)] = S_graphs[g][contig].nodes[node]['dna_seq'][start:]
                        if S_graphs[g][contig].nodes[node]['strand'] == FWD:
                            S_graphs[g][contig].nodes[node]['annotated_aligned_gene_fiveprime_index_{}'.format(ann_spec)] = S_graphs[g][contig].nodes[node]['index'] + start
                        elif S_graphs[g][contig].nodes[node]['strand'] == REV:
                             S_graphs[g][contig].nodes[node]['annotated_aligned_gene_fiveprime_index_{}'.format(ann_spec)] = S_graphs[g][contig].nodes[node]['index'] + GENE_LENGTH[g] - start -1

                else:
                    if g == V:
                        S_graphs[g][contig].nodes[node]['is_annotated_{}'.format(ann_spec)] = False
                        S_graphs[g][contig].nodes[node]['annotated_aligned_gene_{}'.format(ann_spec)] = None
                        S_graphs[g][contig].nodes[node]['annotated_aligned_gene_fiveprime_index_{}'.format(ann_spec)] = None
print("Done")

print("Analyzing Clusters...",end =" ")
#Finding clumps and clusters
clumps = {g: {contig : {ccid : get_clumps_from_component(S_graphs[g][contig],ccid,connected_components[g][contig]) \
                        for ccid in connected_components[g][contig]} \
              for contig in S_graphs[g]} \
          for g in GENE_TYPES_TOFIND}
clusters = {g: {contig : get_clusters_from_clumps(S_graphs[g][contig] , clumps[g][contig], connected_components[g][contig]) for contig in S_graphs[g]} for g in GENE_TYPES_TOFIND}

#Filling cluster information
cluster_information = {}
cluster_information['CL'] = {g: {contig : {clsid: get_group_CL(S_graphs[g][contig],clusters[g][contig][0][clsid]) for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

cluster_information['density'] = {g: {contig : {clsid: get_group_density(S_graphs[g][contig],clusters[g][contig][0][clsid]) for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

cluster_information['PI'] = {g: {contig : {clsid: get_group_PI(S_graphs[g][contig],clusters[g][contig][0][clsid]) for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

cluster_information['span'] = {g: {contig : {clsid: get_group_span(S_graphs[g][contig],clusters[g][contig][0][clsid]) for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

cluster_information['conservation_human'] = {g: {contig : {clsid: get_group_conservation(S_graphs[g][contig],clusters[g][contig][0][clsid],'conservation_human') for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

cluster_information['conservation_species'] = {g: {contig : {clsid: get_group_conservation(S_graphs[g][contig],clusters[g][contig][0][clsid],'conservation_species') for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

cluster_information['AI_human'] = {g: {contig : {clsid: get_group_AI(S_graphs[g][contig],clusters[g][contig][0][clsid],'human') for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

cluster_information['AI_species'] = {g: {contig : {clsid: get_group_AI(S_graphs[g][contig],clusters[g][contig][0][clsid],'species') for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

cluster_information['locus'] = {g: {contig : {clsid: get_group_family(S_graphs[g][contig],clusters[g][contig][0][clsid]) for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

cluster_information['size'] = {g: {contig : {clsid: len(clusters[g][contig][0][clsid]) for \
                                           clsid in clusters[g][contig][0]} \
                                 for contig in S_graphs[g]} for \
                             g in GENE_TYPES_TOFIND}

print_results(cluster_information)
print("Done")

import os
import sys
from Bio import SeqIO

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]

records = [r for r in SeqIO.parse(in_fasta, 'fasta')]

min_length = 50
gene_names = set()
out_f = open(out_fasta, 'w')

for r in records:
    splits = r.id.split('|')
    gene_allele = splits[1]
    gene_name = gene_allele.split('*')[0]
#    if gene_name in gene_names:
#        continue
    r.seq = str(r.seq)
    seq = ''.join([c for c in r.seq if c != '.']).upper()
    print(gene_allele, len(seq))
    if len(seq) <= min_length:
        continue
    out_f.write('>' + gene_allele + '\n' + seq + '\n')
    gene_names.add(gene_name)

out_f.close()

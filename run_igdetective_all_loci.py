import os
import sys
import shutil

def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

def ProcessLocus(genome_fasta, locus, locus_dir):
    # iteration 0 = IgDetective
    # iteration 0 = find genes by alignment
    # combine genes
    for i in range(1, 11):
        # iteration i
        # run find genes by alignment
        # compare genes
        # if there are no new genes, break
        # reassign genes

def main(genome_fasta, output_dir):
    PrepareOutputDir(output_dir)
    loci = ['IGH', 'IGK', 'IGL']
    for locus in loci:
        locus_dir = os.path.join(output_dir, locus)
        os.mkdir(locus_dir)
        ProcessLocus(genome_fasta, locus, locus_dir)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Invalid argument')
        return 
    genome_fasta = sys.argv[1]
    

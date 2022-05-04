# IgDetective
a tool for annotation of immunoglobulin genes (V, D and J genes) in genome assemblies

## System Requirements
* Anaconda package manager 4.8.4 or higher
* Bioconda channel for Anaconda. See [bioconda setup](https://bioconda.github.io/user/install.html#install-conda) for setup guide

## Usage

IGDetective operates in two stages. It first establishes candidate RSSs based on similarity to RSSs' of a known reference species. After this step it detects genes flanking candidate RSSs. Both stages are run through the python script in IGDetective.py :  `python3 IGDetective.py -h -i -o -r -m -g` . The flags are explained below:

* `-h , --help` : Get information on usage and flags
* `-i, --input_file` : provide an input fasta file containing a genome assembly for gene detection
* `-o, --output_directory` : (optional) provide an output directory for generated results. Default location is in the parent directory of the input file
* `-r, --rss_only` : (optional) switch to RSS finding mode
* `-m, --multi_process` : (optional) provide number of parallel processing units if available. Default is 1
* `-g, --genes_type` : (optional) specify which genes (v,d,j) to find. Eg: vdj, d, vj, jv. Default is vdj

## Examples

We have provided 3 sample input species - human (Homo Sapiens), cow (Bos Taurus) and mouse(Mus Musculus) IGH assemblies  in the *example/* directory. We have run IGDetective on the human instance using the command `python3 IGDetective.py -i examples/human.fasta -g vdj-m 20` and generated an annotation of the V,D and J genes in the similarly named directory. 

## Generated result

The predicted V(D,J) genes are listed in a comma seperated file, `genes_V(D,J).tsv`  in the output directory. 
The headers are explained below. Note that all indexing is done with respect the the forward strand inputted by the user : 

For V and J genes
1. `contig` : Name of the input fasta contig in which the gene is detected
1. `reference strand` : Direction of input contig. `+` indicates standard, `-` indicates reverse complement
1. `heptamer index` : starting index of the RSS heptamer
1. `nonamer index` : starting index of the RSS nonamer
1. `heptamer ` :  RSS heptamer sequence
1. `nonamer ` :  RSS nonamer sequence
1. `best aligned human gene ` :  The human gene with which the predicted gene best aligns with, measured as percent identity. The human V(J) genes are listed in `datafiles/human_V(human_J).fasta`.
1. `alignment direction` : Direction of the above mentioned alignment
1. `algnment PI` : Percent identity of the above mentioned alignment
1. `longest common k-mer` : Longest common substring shared by both the best aligned human gene and the predicted gene
1. `start of gene` : Starting index of gene
1. `end of gene` : Ending index of gene
1. `gene sequence` : Nucleotide sequence of predicted gene

For D genes
1. `contig` : Name of the input fasta contig in which the gene is detected
1. `reference strand` : Direction of input contig. `+` indicates standard, `-` indicates reverse complement
1. `left(right) heptamer index` : starting index of the RSS left(right) heptamer
1. `left(right) nonamer index` : starting index of the RSS left(right) nonamer
1. `left(right) heptamer ` :  RSS left(right) heptamer sequence
1. `left(right) nonamer ` :  RSS left(right) nonamer sequence
1. `start of gene` : Starting index of gene
1. `end of gene` : Ending index of gene
1. `gene sequence` : Nucleotide sequence of predicted gene

## Development

IGDetective is still in its nascent stages with room for improvement. Please send your suggestions or any bugs detected to vsirupur@ucsd.edu

# IgDetective
a tool for annotation of variable (V), diversity (D), and joining (J) immunoglobulin genes in genomes.

## System Requirements
* Anaconda package manager 4.8.4 or higher. To install Bioconda, see the [setup guide](https://bioconda.github.io/user/install.html#install-conda).
* Minimap2. To install Minimap2, see the [setup guide](https://anaconda.org/bioconda/minimap2).

## Usage
IGDetective takes a genome in FASTA format as an input and operates in three stages:
* Finding contigs containing IG gene matches. This step is performed using minimap2 and usually takes several minutes. 
* Detecting candidate RSSs based on similarity to RSSs of known reference species. 
* Detects candidates of IG genes flanking candidate RSSs. 

To run IGDetective, type:
```
python run_iterative_igdetective.py genome.fasta output_dir
```
Please note that **IGDetective overwrites the output directory**, so make sure that it does not contain important files.

## Output format

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
IGDetective is still in its nascent stages with room for improvement. Please report any bugs to GitHub. We welcome your comments and suggestions on IGDetective development. Please feel free to send it to Vikram Sirupurapu (vsirupur@ucsd.edu) and/or Yana Safonova (ysafonova@cs.jhu.edu).

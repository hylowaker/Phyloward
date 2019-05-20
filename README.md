Phyloward (#Working in Progress!)
================

Phyloward is a tool for phylogenetic analysis using single-copy prokaryotic core genes.

The SSU rRNA (16S) gene has played an essential role in prokaryotic taxonomy.
However, the short sequence of this single gene is not informative enough
  to acquire stability and resolution of the phylogeny analysis.
Here, we present the set of bacterial and archaeal core genes covering all phyla.

The package provides the following features:
* Extraction of prokaryotic core genes from genome assemblies
* Multiple alignment of core genes
* Concatenation of core gene sequences
* Customization of the pipeline by alternative model or parameters
* Phylogenetic analysis using FastTree
* Calculation of Gene Support Index (GSI) which indicates how many genes
  support the branch in the concatenated phylogenetic tree

## Getting Started

### Prerequisites

**OS**:  
Linux, macOS  
(It will also work on Windows in theory, but we haven't tested. For Windows users, we recommend using [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10) instead)

**Python Interpreter**:  
Python 3.5 or higher from [official Python website](https://www.python.org)
  or [Anaconda Distribution](https://www.anaconda.com/download/).  
*Note: Python 2 is __not__ supported.*

**Third party dependencies**:  
The package depends on several external programs.
The following programs should be installed in advance.
* [Prodigal](https://github.com/hyattpd/prodigal/releases/) (>= 2.6.3) : Gene prediction
* [HMMER](http://hmmer.org/download.html) (>= 3.1) : Core gene identification
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) (>= 7.402) : Multiple alignment
* [FastTree](http://www.microbesonline.org/fasttree/#Install) (>= 2.1.10) : Phylogeny tree construction
* ~~RaxML~~ (Not implemented yet)

Make sure that each of `prodigal`, `hmmsearch`, `mafft`, 
  and `FastTree` command is executable in shell.

Alternatively, if you are using Linux or macOS and have `conda` installed, you can install dependencies by this one-liner:
```
$ conda install -c bioconda prodigal hmmer mafft fasttree
```

### Download and Installation 
Clone this project with `git`.

    $ git clone https://github.com/hylowaker/Phyloward.git

The package can be installed by `pip3`, a package manager for Python 3. (If your default pip is Python 3 version, you may use `pip` instead of `pip3`)

    $ cd Phyloward
    $ pip3 install .

### Basic command line use
Try invoking the program in shell.

    $ phyloward

The above command will simply print help message and exit. 
The actual pipeline can be executed with two subcommands: `extract` and `align`

1.  First, you have to identify core genes from a genome by `phyloward extract` command.
    The following command lets the core genes information from _such.fasta_ file
    be stored in _doge_ directory as JSON format.
     
         $ phyloward extract such.fasta doge
     
    Alternatively, a batch of multiple fasta files in a same directory can be processed.
    The command below will allow the program to take all FASTA files in _much_ directory as input.
         
         $ phyloward extract much doge
    
    Note: If you are dealing with archaeal genome, add `--archaea` option.

2.  After several genomes were processed with `phyloward extract`,
    you can align those core gene sequences to infer phylogenetic trees using `phyloward align` command.

         $ phyloward align --tree doge wow
    
    (If you want only alignments, not trees, remove `--tree` option)  
    Using the extracted core gene data stored in _doge_ directory,
      the alignments and the trees will be saved in _wow_ directory
      as FASTA and Newick format, respectively.
    Another single JSON file with all results and metadata joined together will also be created.
   


## Import as Python package
It is also possible to import the package and run the pipeline in your Python code.

An example code to extract core genes from single genome file is shown below.
```python
from phyloward import extractor

file = 'some_bacteria.fasta'
extracted = extractor.extract_core_genes(file)
print(extracted.get_best_hit_sequence('rpsC'))
```
    
For further information, please check documentation. _(TODO)_


## Documentation
_(We're working on it)_


## Citation
_We're working on a manuscript for Phyloward now._


## Copyright and License
Copyright (c) 2018-2019 JaeHeung Han  
Phyloward is available under the terms of [GPLv3](LICENSE).

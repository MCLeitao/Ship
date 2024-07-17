# SHIP • Safe Harbor Identification Program
SHIP is Python program for the identification of Genomic Safe Harbors (GSHs) in eukaryotic organisms.

SHIP software was written in the Python 3.9 programming language using standard libraries, which implies a simplified installation process and transparent use on multiple platforms, such as Linux and Windows.

The input consists of (1) the genome annotation file of the target organism, in the GFF3 format (Generic Feature Format); and the Features and [Track](https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html) files in JSON format. 

For the identification of GSHs based on biological premises, a series of annotation attributes were selected to direct a search for regions with a minimum size and devoid of genes, but also that can harbor regulatory contexts, such as binding sites for transcription factors, DNA hypermethylation and state of chromatin and histone alterations, which are evaluated using information collected in external databases through the [Entrez Programming Utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) (E-utilities) Application Programming Interface (API) National Center for Biotechnology Information (NCBI), [Ensembl REST API](https://rest.ensembl.org/) and the [UCSC Genome Browser](https://genome.ucsc.edu/).

The software performs a combinatorial analysis to identify **intergenic regions** and the direction of their flanking genes. Subsequently, intergenic regions presenting the genetic elements described in the track.json are excluded.

## Contents


## Description


## Features


## Usage



## Data

## Usage

## Reference
Please cite the following preprint when referencing Ship.


## Copyright



## License

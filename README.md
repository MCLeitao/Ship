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

### Inputs

### Navigating

```
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   Welcome on Board Sailor!   
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Enter your email to access NCBI: user-name@email.com
Enter file name <include the path>: /Users/user1/Documents/target_organism.gff3
Enter the interval for genomic analysis (default 500): (enter)
```
```
------------------------------
ANALYZING ...
------------------------------
Performing the following steps:
    [1] Selecting and sorting the genes by their start within each chromosome.
    [2] Disregarding completely overlapping genes.
    [3] Filtering out convergent, divergent and tandem genes.
------------------------------
features : ['gene', 'ncRNA_gene', 'pseudogene', 'centromere', 'telomere', 'long_terminal_repeat', 'mobile_genetic_element', 'origin_of_replication', 'transposable_element_gene', 'meiotic_recombination_region', 'sequence_feature']
------------------------------
BK006935.2, BK006935.2 BK006936.2, BK006936.2 BK006937.2, BK006937.2 BK006938.2, BK006938.2 BK006942.2, 
BK006942.2 AJ011856.1, AJ011856.1 BK006939.2, BK006939.2 BK006940.2, BK006940.2 BK006941.2, BK006941.2 BK006934.2, 
BK006934.2 BK006943.2, BK006943.2 BK006944.2, BK006944.2 BK006945.2, BK006945.2 BK006946.2, BK006946.2 BK006947.3, 
BK006947.3 BK006948.2, BK006948.2 BK006949.2, BK006949.2 

---------- Types that will be processed ----------
gene = 6600        transposable_element_gene = 91
ncRNA_gene = 424   pseudogene = 12

-------- Types that will not be processed --------
chromosome = 17   CDS = 6913   snoRNA = 77                   snRNA = 6  
mRNA = 6600       ncRNA = 18   transposable_element = 91     five_prime_UTR = 4
exon = 7507       tRNA = 299   pseudogenic_transcript = 12   rRNA = 24
--------------------------------------------------
   Interval (bp)                   Tandem            Divergent           Convergent
         Overlapping                 142                 264                 184
            1-500 bp                2077                 842                1445
        501-1,000 bp                 698                 518                 154
      1,001-1,500 bp                 169                 124                  46
      1,501-2,000 bp                  74                  56                  10
      2,001-2,500 bp                  26                  22                   4
      2,501-3,000 bp                  11                  13                   1
      3,001-3,500 bp                   6                   7                   1
      3,501-4,000 bp                   8                   4                   1

Quantity of Chromosomes = 17
Number of genes located = 7127
Total intergenic intervals = 6925
Total Intervals flanked by Tandem genes     = 3218
Total Intervals flanked by Divergent genes  = 1855
Total Intervals flanked by Convergent genes = 1852
Number of complete overlaps between genes   = 191
```
```
Would you like to plot the Genes chart? (Y/N): Y
Types of flanking genes: 
[1] CONVERGENT (+ -) 
[2] DIVERGENT (- +) 
[3] TANDEM (+ +/- -)
> Enter your choice: 1
> Minimum size of the intergenic region (Standard 1500): 1300
> Maximum size of the intergenic region (Standard 2000): 2200
```
```
Do you want to perform analysis on the UCSC Database? (Y/N) Y
Do you want to perform Cross References? (Y/N) Y
Do you want to perform regulatory analysis? (Y/N) Y
Enter the species to be analyzed: _Target organism_
Enter the regulatory build Path and File name: /Users/user1/Documents/target_organism.Regulatory_Build.regulatory_features.gff
```
```
------------------------------
PROCESSING ...
------------------------------
[1] Selecting the intergenic regions with size between 1500 to 2000 bp.
[2] Searching for sequences at NCBI.
[3] Crossing information with other databases.
Wait for processing... Go have a cup of coffee.
☕︎︎ ☕︎︎ ☕︎︎ ☕︎︎ ☕︎︎ ☕︎︎ ☕︎︎ ☕︎︎
```
```
------------------------------
RESULT
------------------------------
File generated with the results: /Users/user1/Documents/target_organism_result_SHIP.txt
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Do you want to analyze another GFF file? ...
```

### Outputs

Genomic Overview
![Saccharomyces_cerevisiae R64-1-1 49 _plot](https://github.com/user-attachments/assets/ee3f9a87-6d5d-4fbe-9ff5-0cd766eafcdd)

pGSH List (without Regulatory Analysis)
<img width="1002" alt="Screenshot 2024-07-18 at 16 07 22" src="https://github.com/user-attachments/assets/53abd587-9cec-4e4a-bdf8-7d0b2ad367f5">

pGSH List (with Regulatory Analysis)




## Data

## Usage

## Reference
Please cite the following preprint when referencing Ship.


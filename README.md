# :anchor: SHIP: Safe Harbor Identification Program
_Software for the identification of Genomic Safe Harbors (GSHs) in eukaryotic organisms_<br/>

## Description

A Software with a graphical user interface written in Python 3.9 using standard libraries, which implies a simplified installation process and a possibility of being used on multiple platforms, such as Linux and Windows.

The input consists of (1) the genome annotation file of the target organism, in the GFF3 format (Generic Feature Format); and (2) the Features and [Track](https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html) files in JSON format.

For the identification of GSHs based on biological premises, a series of annotation attributes were selected to direct a search for regions with a minimum size and devoid of genes, but also that do not harbor regulatory elements, such as binding sites for transcription factors, DNA hypermethylation and state of chromatin and histone alterations, which are evaluated using information collected in external databases through the [Entrez Programming Utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) (E-utilities) Application Programming Interface (API) National Center for Biotechnology Information (NCBI), [Ensembl REST API](https://rest.ensembl.org/) and the [UCSC Genome Browser](https://genome.ucsc.edu/).

The software performs a combinatorial analysis to identify **intergenic regions** and the direction of their flanking genes. Subsequently, intergenic regions presenting the genetic elements described in the track.json are excluded.


## Contents
- [Dependencies](https://github.com/MCLeitao/Ship#dependencies)
- [Usage](https://github.com/MCLeitao/Ship#usage)
  - [Inputs](https://github.com/MCLeitao/Ship#inputs)
  - [Outputs](https://github.com/MCLeitao/Ship#outputs)
- [Navigating](https://github.com/MCLeitao/Ship#navigating)
- [Reference](https://github.com/MCLeitao/Ship#reference)
- [License](https://github.com/MCLeitao/Ship#license)


## Dependencies

* Python 3.9
* Packages (version): bio (1.5.9); biopython (1.81); biothings-client (0.3.0; certifi (2023.5.7); charset-normalizer (3.1.0); contourpy (1.0.7); cycler (0.11.0); fonttools (4.39.4); gprofiler-official (1.0.0); idna (3.4); kiwisolver (1.4.4); matplotlib (3.7.1); mygene (3.2.2); numpy (1.24.3); packaging (23.1); pandas (2.0.2); Pillow (9.5.0); platformdirs (3.5.1); pooch (1.7.0); pyparsing (3.0.9); python-dateutil (2.8.2); pytz (2023.3); requests (2.31.0); six (1.16.0); tqdm (4.65.0); tzdata (2023.3); urllib3 (2.0.3).


## Usage

### Inputs

To use the program, you must have an NCBI account due to the Database API Usage Guidelines and Requirements ([Entrez](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen)). The email registered in this account is required.


<br/><ins>Before running the program you will need all the following files in the same directory</ins>

1. The genome annotation file of the target organism, in the GFF or [GFF3](https://gmod.org/wiki/GFF3) format. (Ex. [`Saccharomyces_cerevisiae.R64-1-1.49.gff3`](examples/yeast/Saccharomyces_cerevisiae.R64-1-1.49.gff3))
   * If it is a multicellular organism, there is the possibility of adding a GFF file with the [Regulatory Features](https://www.ensembl.org/info/genome/funcgen/index.html) of each cell type (Ex. `homo_sapiens-regulatory_features.gff`)

2. The Feature list that SHIP will use to analyze the GFF file in JSON format. (Ex. [`features.json`](ship/features.json))
   * At the beginning, the software will show you all the features present in the GFF file

3. The [Track](https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html) list that SHIP will use to check the presence in each intergenic region in UCSC Genome Browser Database in JSON format. (Ex. [`tracks.json`](examples/yeast/tracks.json))

> [!WARNING]
> Do not change the name of the JSON files.


<br/><ins>During running the user needs to choose the following key parameters</ins>

1. Neighbour gene orientation
   * Convergent (+ -)
   * Divergent  (- +)
   * Tandem     (+ + / - -)

2. Intergenic size range
   * Maximum size of the intergenic region (Standard 2000)
   * Minimum size of the intergenic region (Standard 1500)

The user also needs to indicate the analyses that should be performed

1. UCSC Database
2. Ensembl (Cross References)
3. Regulatory analysis (Multicellular organisms)

---
### Outputs
All files will be generated in the same directory as the input files.

1. Genomic Overview
  ![Saccharomyces_cerevisiae R64-1-1 49 _plot](https://github.com/user-attachments/assets/ee3f9a87-6d5d-4fbe-9ff5-0cd766eafcdd)

2. pGSH List
   * Description of neighboring genes
   * Fast sequence of the intergenic region
   * Description of Regulatory
   * Additional information for each neighboring gene

  _Saccharomyces cerevisiae_
  <img width="1002" alt="Screenshot 2024-07-18 at 16 07 22" src="https://github.com/user-attachments/assets/53abd587-9cec-4e4a-bdf8-7d0b2ad367f5">

  _Homo sapiens_ (with Regulatory Analysis)
  ![Screenshot 2024-07-19 at 19 24 25](https://github.com/user-attachments/assets/0c40effb-26db-4847-80ee-b7ad23b38342)


## Navigating
_Practical example running the software_

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

## Reference
Leitão et al. SHIP: A Computational Tool for the Systematic Identification of Genomic Safe Harbors for Stable Gene Insertion in Eukaryotes. In Review

##  License
* MIT License
* Copyright (c) 2023 Matheus de Castro Leitão
* Software registered under process number BR512023002017-6

# WGS Analysis and Detection of Molecular Markers (WADE)

## What is WADE?

WADE provides a flexible and customizable method to extract specific genes from a large number of genomes at once using BLAST to interrogate assembled genomes, current molecular analyses include antimicrobial resistance, toxins, virulence profiles and several multi-locus sequence typing (MLST) schemes. The Virulence Factor DataBase (VFDB) and antimicrobial resistance factor databases CARD, ARG-ANNOT and ResFinder have also been made available.

Tabular results are output in a format that is compatible for LabWare uploads. These results can consist of simple “Positive-Negative” results corresponding to presence or absence of a queried gene. Curated multi-fasta lookup files can be provided for molecular determinants to create molecular profiles of affective mutations. Fasta file outputs of gene sequences extracted from the genomes can readily be loaded into sequence aligners to correlate nucleotide differences to phenotypic observations.

## Getting Started

This tool can be run using RStudio (available at https://www.rstudio.com/)

### Prerequisites

This tool requires the use of R packages: plyr, dplyr, tidyverse, tidyselect, stringr, Biostrings which can be loaded using:

  ```sh
  library(plyr)
  library(dplyr)
  library(tidyverse)
  library(tidyselect)
  library(stringr)
  library(Biostrings)
  library(shiny)
  library(DT)
  library(here)
  library(wade)  
  ```
and the use of the BLAST+ executable from NCBI: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download


### Installation

1. Install R from https://www.r-project.org/
2. Install RStudio from https://www.rstudio.com
3. Install WADE:  
   a. from github using:
   ```sh
   install.packages("githubinstall")
   library(githubinstall)
   githubinstall("wade")
   ```  
   OR
   b. clone the wade git repository into a directory and run the following:
   ```sh
   install.packages("C:\\path\\to\GitHub\\wade", repos=NULL, type="source") 
   ```  
4. Install required packages
   ```sh
   install.packages("plyr")
   install.packages("dplyr")
   install.packages("tidyverse")
   install.packages("tidyselect")
   install.packages("stringr")
   install.packages("shiny")
   install.packages("DT")
   install.packages("here")
   ```
5. Install Biostrings (https://bioconductor.org/install/)
   ```sh
   if(!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
   BiocManager::install(c("ggtree", "Biostrings"))
   ```
6. Install the BLAST+ executable from NCBI: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/  
Download and install the *.exe file from the FTP directory.  

## Setup

Copy the [WADE](WADE) directory from the github/wade folder into the C:/ drive  
OR  
Set up the following directory
```sh
C:/WADE
```
with two subdirectories:  
```sh
C:/WADE/Output  
C:/WADE/temp  
```

This molecular analysis tool queries pre-assembled fasta files. The location of the contig files, vcf files, and wade data files need to be listed in "DirectoryLocations.csv" in the line for the corresponding organism ID:  

*   **OrgID** = The organism being queried. *NOTE: Do not edit this column*  
*   **LocalDir** = C:\\WADE\\   
*   **SystemDir** = The directory in which the local wade-data files are saved  
*   **ContigsDir** = The directory in which contigs are saved  
*   **VCFDir** = The directory in which vcf files are saved (GONO and PNEUMO Only)  

The contig files must have the file extension".fasta" (eg. MySampleNo_contig.fasta).

To use the multiple sample list option, a sample list file must be located at: C:/WADE/list.csv
list.csv must have the following structure:

|SampleNo	|Variable |
|---------|---------|
| 12345	  | 4 ug/ml |   
| 12346	  | 8 ug/ml |  

## Usage  
  
1. Open the WADE.R file in RStudio. WADE.R is an RShiny UI interface to facilitate usage of WADE.
2. change the current working directory in WADE.R to the folder where your "DirectoryLocations.csv" file and Output and temp subdirectories are located.
```sh
Line 20: curr_work_dir <- "C:\\WADE\\"
```
3. Click on the "Run App" button
4. Select the desired organism  
5. Select the desired analysis  
6. choose the query locus. The default is "list" to query all loci in the chosen category, otherwise type in the desired locus.  
7. If this is the first time you have run this particular analysis or if the allele lookup files have been recently updated, the BLAST databases must be indexed by pressing the "MakeBlastdb" button.
8. Enter sample number, or "list" to query multiple samples  
9. Click the "Go" button. The Output button will display a spreadsheet containing the results of the analysis. These results can also be found in the Output folder.  
  
## Troubleshooting

When running this program on some Windows machines, the MakeBlastdb program can give an error. If this happens, the environmental variables setting will need to be changed as follows:  

1. Go to Windows Settings and search for "Environmental Variables"  

2. In the System Properties dialogue box, click on the "Environmental Variables" button  

3. In the "User Variables for..." box, click "New..." button  

4. Input the following:    
   ```sh  
   Variable Name: BLASTDB_LMDB_MAP_SIZE
   Variable Value: 1000000
   ```

## Legal

Copyright Government of Canada 2022

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Contact

**Walter Demczuk:** Walter.Demczuk@phac-aspc.gc.ca

**Shelley Peterson:** Shelley.Peterson@phac-aspc.gc.ca

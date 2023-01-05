## What is MasterBlastR?

MasterBlastR provides a flexible and customizable method to extract specific genes from assembled contigs.

Using BLAST to interrogate assembled genomes, MasterBlastR is able to extract gene sequences for molecular analyses including genes for antimicrobial resistance, toxins, virulence profiles and several multi-locus sequence typing (MLST) schemes. Gene mutations and alleles can be identified by comparison with a multi-fasta file containing sequences for each allele (i.e. for MLST allele sequences) or using a file listing mutations (amino acid locations) or both. Fasta file outputs of gene sequences extracted from the genomes can readily be loaded into sequence aligners to correlate nucleotide differences to phenotypic observations.

## Getting Started

This tool can be run using RStudio (available at <https://www.rstudio.com/>)

### Prerequisites

This tool requires the use of R packages: plyr, dplyr, tidyverse, tidyselect, stringr, Biostrings which can be loaded using:

``` sh
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyselect)
library(stringr)
library(Biostrings)
```

and the use of the BLAST+ executable from NCBI: <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>

### Installation

1.  Install R from <https://www.r-project.org/>

2.  Install RStudio from <https://www.rstudio.com>

3.  Install required packages

    ``` sh
    install.packages("plyr")
    install.packages("dplyr")
    install.packages("tidyverse")
    install.packages("tidyselect")
    install.packages("stringr")
    ```

4.  Install Biostrings (<https://bioconductor.org/install/>)

    ``` sh
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(c("ggtree", "Biostrings"))
    ```

5.  Install the BLAST+ executable from NCBI: <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>

## Usage

1.  Ensure MasterBlastR is located in a folder with subfolders "output", "temp", "wildgenes" and "allele_lkup_dna".

2.  Ensure blast_evalues.csv is located in the same folder as MasterBlastR\
    blast_evalues.csv must have the following structure:

    ``` sh
    contig     allele      wt_id 
    1.00E-50   1.00E-98    10
    ```

3.  Set the working directory where MasterBlastR is located.

    ``` sh
    line 20: curr_work_dir <- "C:\\MasterBlastR\\"
    ```

4.  This molecular analysis tool queries pre-assembled fasta files. The location of the contig files needs to be assigned to ContigsDir with the file extension ".fasta" (eg. MySampleNo_contig.fasta).

    ``` sh
    line 21: ContigsDir <- "C:\\MasterBlastR\\contigs\\"
    ```

5.  Put reference gene sequences into subfolders. File names need to be the same for both files (i.e. parC.fasta):

    a.  wild-type into "wildgenes"\
    b.  allele lookup multi-fasta files into "allele_lkup_dna" (optional)

6.  Fasta headers for reference gene sequences should be in the following format:

    ``` sh
    >Locusname_AlleleNum_Mutation_Comment_Variable
    ex. >folA_28_I100L_Sample5_NA
    ```

7.  To query genes for which an allele lookup multi-fasta file is unavailable, ensure the "wildgenes" folder contains a .fasta of the wild-type gene and create a loci mutations list in the directory. (eg. C:/MasterBlastR/loci_mutations.csv).\
    loci_mutations.csv must have the following structure:

    ``` sh
    Locus_id   Name      Posn_1    Posn_2    WildType  
    gene1      S81       81        81        S
    gene2      ARA100    100       102       ARA
    ```

8.  Gene names must be consistent between wildgenes file name, lookup table file name, Locus_id in loci_mutations.csv, and query.

9.  To use the multiple sample list option, a multiple sample list file must be located in the directory. (eg. C:/MasterBlastR/list.csv)\
    list.csv must have the following structure:

    ``` sh
    SampleNo     Variable   
    12345          4 ug/ml   
    12346          8 ug/ml   
    ```

10. To use the multiple locus list option, a multiple locus list file must be located in the directory. (eg. C:/MasterBlastR/loci.csv). loci.csv must have the following structure:

``` sh
Locus_id
gene1
gene2
gene3
```

## Troubleshooting

When running this program on some Windows machines, the makeblastdb program can give an error. If this happens, the environmental variables setting will need to be changed as follows:

1.  Go to Windows Settings and search for "Environmental Variables"

2.  In the System Properties dialogue box, click on the "Environmental Variables" button

3.  In the "User Variables for..." box, click "New..." button

4.  Input the following:

    ``` sh
    Variable Name: BLASTDB_LMDB_MAP_SIZE
    Variable Value: 1000000
    ```

## Legal

Copyright Government of Canada 2022.

Written by: National Microbiology Laboratory, Public Health Agency of Canada.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Contact

**Walter Demczuk:** [Walter.Demczuk\@phac-aspc.gc.ca](mailto:Walter.Demczuk@phac-aspc.gc.ca){.email}

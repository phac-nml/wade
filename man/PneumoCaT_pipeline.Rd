% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PneumoCaT.R
\name{PneumoCaT_pipeline}
\alias{PneumoCaT_pipeline}
\title{PneumoCaT_Batch (Waffles version) :
Log in to Waffles
Download files from IRIDA using ngsArchiveLinker using folder structures
Rename files to remove _001.fastq using: prename 's/_001.fastq/.fastq/' */*.fastq
Run command line to send all samples in the "fastqs" folder :
$ for sample in `ls fastqs` ; do  sbatch -c 1 -p high --mem 4G --output $sample-%j.out --wrap="PneumoCaT.py -i fastqs/$sample -o Output/$sample -s /share/apps/samtools-0.1.19/samtools" ; done
$ watch sq
CTRL-C to exit watch.}
\usage{
PneumoCaT_pipeline(sample_num)
}
\arguments{
\item{sample_num}{Sample number associated with contig.fasta file}
}
\value{
A table frame containing the results of the query
}
\description{
PneumoCaT_Batch (Waffles version) :
Log in to Waffles
Download files from IRIDA using ngsArchiveLinker using folder structures
Rename files to remove _001.fastq using: prename 's/_001.fastq/.fastq/' */*.fastq
Run command line to send all samples in the "fastqs" folder :
$ for sample in `ls fastqs` ; do  sbatch -c 1 -p high --mem 4G --output $sample-%j.out --wrap="PneumoCaT.py -i fastqs/$sample -o Output/$sample -s /share/apps/samtools-0.1.19/samtools" ; done
$ watch sq
CTRL-C to exit watch.
}

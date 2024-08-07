% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rename_Submitted.R
\name{rename_submitted}
\alias{rename_submitted}
\title{Renames submitted contig files to NML LabWare numbers
February 5 2024, Walter Demczuk & Shelley Peterson}
\usage{
rename_submitted(Org_id, curr_work_dir)
}
\arguments{
\item{curr_work_dir}{Start up directory from pipeline project to locate system file structure}

\item{Org_code}{Organism to query: GAS, PNEUMO or GONO}
}
\value{
A table frame containing the results of the query
}
\description{
This script makes a new LabWareUpload_metrics_SUBMITTED.csv with submitted lab nos replaced with LabWare nos.
}
\details{
Copies Contig Assemblies and VCF files from submitted IRIDA sample no to LabWare No.
Reads in a list of new LabWare Nos. and submitted genome sample numbers from list.csv
Makes a new LabWareUpload_metrics_SUBMITTED.csv with subm.lab.nos replaced with LabWare Nos.

Input:
SampleNo    Variable
LabWare id  Submitted id in IRIDA

SampleNo is new LabWare No
Variable is the submitted IRIDA sample id.
}

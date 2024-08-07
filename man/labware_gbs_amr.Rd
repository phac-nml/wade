% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LabwareUpload_GBS_AMR.R
\name{labware_gbs_amr}
\alias{labware_gbs_amr}
\title{Labware Upload Formatter for GBS AMR
February 5 2024, Walter Demczuk & Shelley Peterson
Run AMR first
Then run this analysis to combine data the full amr profile to upload to LabWare.}
\usage{
labware_gbs_amr(Org_id, curr_work_dir)
}
\arguments{
\item{Org_id}{Organism to query: GAS, PNEUMO or GONO}

\item{curr_work_dir}{Start up directory from pipeline project to locate system file structure}
}
\value{
A table frame containing the results of the query
}
\description{
Labware Upload Formatter for GBS AMR
February 5 2024, Walter Demczuk & Shelley Peterson
Run AMR first
Then run this analysis to combine data the full amr profile to upload to LabWare.
}

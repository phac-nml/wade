## What is SHINY_Trees?

SHINY_Trees is an RShiny-Based UI for viewing and annotating phylogenetic trees.

Using a Newick tree file and associated sample metadata as input, SHINY_Trees can generate phylogenetic trees with colour coded tree tips in radial or rectangular format, and can also generate temporal plots to complement the trees.

## Getting Started

This tool can be run using RStudio (available at <https://www.rstudio.com/>)

### Prerequisites

This tool requires the use of R packages: tidyverse, scales, shiny, shinyWidgets, ggtree, phytools, and gridExtra which can be loaded using:

``` sh
library(tidyverse)
library(scales)
library(shiny)
library(shinyWidgets)
library(ggtree)
library(phytools)
library(gridExtra)
```

### Installation

1.  Install R from <https://www.r-project.org/>

2.  Install RStudio from <https://www.rstudio.com>

3.  Install required packages

``` sh
install.packages("tidyverse")
install.packages("scales")
install.packages("shiny")
install.packages("shinyWidgets")
install.packages("phytools")
install.packages("gridExtra")
```

4.  Install Biostrings and ggtree (<https://bioconductor.org/install/>)

``` sh
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ggtree", "Biostrings"))
```

## Usage

1.  Set the working directory for SHINY_Trees - This can be the WADE folder or another folder. Ensure the working directory folder has a subfolder called "Output".

``` sh
line 34: curr_work_dir <- "C:/SHINY_Trees/"
```

2.  Ensure ColourList.csv and ShapeList.csv are located in the working directory folder. These files define tip shapes and preferred custom colours and can be altered.

3.  Run the SHINY app by clicking the "Run App" button.

4.  Select plot type from the following options:

    -   **Tree** - creates a rectangular phylogenetic tree
    -   **Radial** - creates a radial phylogenetic tree
    -   **Temporal** - creates a dot plot showing temporal distribution of isolates by date collected
    -   **Tree + Tempora**l - displays rectangular tree and temporal distribution of isolates for a side-by-side comparison

5.  Load **Newick Tree File** and **Metadata.csv file** into corresponding file inputs. \
    Metadata file must have the following columns: **Strain_ID, Date Collected, Year**. *Province is not required, however some errors may occur if this column is omitted.*

6.  Select tree parameters:

    -   **Colour Tips By:** - Any metadata column can be used, but if "Custom Plot Colours" is selected, the number of categories must not be greater than the number of custom colours listed.

    -   **Highlight in Red**: - Will highlight a specific category in red. This will override the colour generated from "Colour Tips By:"

    -   **Colour Other Tips:** - Will change any tips other than those highlighted in red to grey. This overrides "Colour Tips By"

    -   **Apply Shape to:** - Similar to "Colour Tips By" but will apply shapes. Please note, only the first 5 shapes in ShapeList.csv are fillable shapes, so metadata columns with more categories will not be coloured if both shapes and colours are selected.

    -   **Apply Triangles to:** - Similar to "Highlight in Red" will highlight a specific category using the triangle shape.

    -   **Node Size:** - Sliding scale to adjust tree node size

**Please note: if both colours and shapes are selected, the colours will not appear in the legend**

7.  Click **Draw Plot** to display the plot in the RShiny app, or **View Plot** to open the plot as a pdf file.

## Legal

Copyright Government of Canada 2025.

Written by: National Microbiology Laboratory, Public Health Agency of Canada.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Contact

**Shelley Peterson:** [Shelley.Peterson\@phac-aspc.gc.ca](mailto:Shelley.Peterson@phac-aspc.gc.ca){.email}

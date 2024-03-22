[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

# Bioinformatics training - Fighting Malaria Across Borders (FiMAB)

This repository houses the course notes for the Bioinformatics Module of the Fighting Malaria Across Borders (FiMAB) international training programme, created by the Institute of Tropical Medicine (ITM) in Antwerp (Belgium) and supported by VLIR-UOS. Its primary goal is to support the implementation of targeted sequencing assays (in particular, AmpliSeq) to strengthen malaria molecular surveillance and help guide national control programmes. In conjunction with laboratory training, this bioinformatics course is intended to allow young academics around the globe to become familiar with molecular surveillance as a key activity to monitor transmission, sources of epidemics and the emergence and spread of drug resistance mutations in the _Plasmodium_ parasite.

Currently, the course notes focus on an introduction to the Unix command line and data manipulation using R. Further content aimed at sequencing analysis pipeline and population genetics will be added in the future.

The course notes are available at the following URL: [https://pmoris.github.io/FiMAB-bioinformatics/](https://pmoris.github.io/FiMAB-bioinformatics/)

An accompanying GitHub codespace environment containing all the training files and tools, can be accessed in two ways:

1. Bioinformatics training environment (default): a VSCode environment containing a Unix terminal and several common (bioinformatics) bash tools. [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/pmoris/FiMAB-bioinformatics?template=false&devcontainer_path=.devcontainer%2Fdevcontainer.json)
2. RStudio Server (must be selected manually): an RStudio server containing tidyverse, based on the [rocker project](https://rocker-project.org/images/devcontainer/images.html). Once the environment has started (a VSCode session should appear), you can browse to the RStudio IDE URL listed under the ports tab (in the same pane as the terminal) or via the forwarded ports "Radio" icon in the status bar at the bottom of the window (it will look a little bit like this `https://ubiquitous-space-eureka-5pxwr7qg96p2vr79-8787.app.github.dev/`/). [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/pmoris/FiMAB-bioinformatics?template=false&devcontainer_path=.devcontainer%2Frstudio%2Fdevcontainer.json)

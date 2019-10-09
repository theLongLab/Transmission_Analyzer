# Transmission Analyzer

**A pipeline that enables the prediction of HIV transmission and Infection date inference. The provided R code can be utilized to analyze a BEAST2 generated MCC tree to predict novel transmission clusters as well as Infection dates of patients using Sanger sequencing data.**

---

### Table of Contents

- [Description](#description)
- [Getting Started](#getting-started)
- [Function of R and Python Scripts](#function_of_r_and_python_scripts)
- [Procedure](#procedure)
- [License](#license)
- [Authors](#author-info)
- [Citations](#citations)
---

## Description

### Validation of a phylogenetic pipeline to examine transmission networks in a Canadian HIV cohort.

The following pipeline was used to analyze Person to Person (P2P) transmission relationships in 139 patients infected with Type 1 HIV.

The pipeline is capable of inferring not only the transmission relationships but also the dates of infection. In the event where there may be a loss of patient, the pipeline is capable of predicting non sampled patients that may be the source of infection.

The pipeline consists of three steps.

#### Procedure of Pipeline

- Multiple Sequence Alignment
- Rooted Phylogenetic Tree Generation
- Transmission and Infection date inference

Custom made R scripts and Python code was written to conduct the analysis. The process is automated and once the functions are loaded into the R-project you can simply execute the entire process.

- The R scripts main goal is to extract the summary statistics about the topology of the inferred transmission trees and extract inferred person-to-person (P2P) relationships.
- The Python script was written to compare the inferred P2P transmission relationships and infection dates to the clinically known information.

These scripts are highly commented and can be easily modified to extract whatever you want from the TransPhylo ctree object.

[Back To The Top](#table-of-contents)

---

## Getting Started

### Prerequisites
- **Required Software**
  - R 3.0 and above

  - RStudio 1.1 and above
  - Python 2.0 and above

- **R Studio Libraries**
  - TransPhylo
  - ape
  - stringr
  - lubridate
  - readr
  - reticulate
  - igraph


[Back To The Top](#table-of-contents)

---

## Function of R and Python Scripts

#### Main_Run

You can run the entire process from the Main_Run script. Simply edit the parameters as desired and run the script.

#### brLenChecker

First use TransPhylo and generate the ptree from a neewick tree.

>ptree<-ptreeFromPhylo(read.tree('sample.nwk'),dateLastSample=2018.222)

Then use brLenChecker function to identify the row and column of the affected nodes, ID of the originating node, patient ID of the affected node.

If there are no affected nodes the function will return a null value.

Though simulated trees will rarely have this problem but it should be noted that maximum-clade credibility phylogenetic (MCC) trees made from real sequence data may have branches that have zero or negative length.

TransPhylo breaks when this happens, so they have to be removed from the Newick MCC tree file by hand. brLenChecker checks the branch lengths of the phylogenetic tree object made by TransPhylo's ptreeFromPhylo and reports the location of the problematic branches. If there were any, the TransPhylo ptree object must be re-made with the cleaned Newick file.

#### autoSummary

The auto Summary function takes over the TransPhylo inferTTree function. The parameters can be edited using the Main_Run script.

- The inputs inlcude:
 - TransPhylo ptree object

 - Number of sampled patients
 - TransPhylo parameters to be customized

- The outputs are:
 - Inferred TransPhylo ctree object

 - summary statistics in the summary_data dataframe

_It should be noted that autoSummary given an output in the form of a modified ctree, this should be assigned to a R variable and this new ctree should be used in the transReport function_

>c<-autoSummary(ptree,pop,0.99,F,1,gshape,gscalep,sshape,scalep)

#### transReport

Generates a Report of the predicted Transmission relationships and the predicted Infection dates in the form of a text file.

The autoSummary generated c tree and the population size of the sample should be input into the function

>transReport(c,pop)

#### predictVSKnown

This is a python script that can be loaded into R Studio and executed from the same project. This function is automated in the Main_Run script.

There is no input file that has to be entered but the transReport generated tsv file should be left in the same project directory for the script to access. It requires the sample_actual_transmission.tsv rename it as actual_transmission_table.tsv before use and place in project folder.

#### plotCTree

This process is completely automated in the Main_Run script.

Essentially it marks the number of correctly and incorrectly predicted or inferred relationships on the generated phylogenetic tree itself. The function overwrites TransPhylo default plotCTree function.

To manually perform this process you have to extract the correct and incorrect transmissions from the predictvsKnown function's txt output file. The extracted data has to be entered as two separate vectors of correct and incorrect vectors as dates in the decimal format.

To convert dates to decimal you may use the lubridate library.

> plotCTree(c,TRUE,NA,NA,correct,incorrect)

**The output will be as follows:**

![alt text](https://github.com/theLongLab/Transmission_Analyzer/blob/master/Rplot.png)

**Figure 1:** _Transmission tree generated by MUSCLE, adjusted BEAST parameters, and adjusted TransPhylo parameters. The green arrows represent correctly inferred transmissions while the red arrows represent incorrectly inferred transmissions. The squares represent novel predicted patients._

#### clusterGen

clusterGen is used to visualize the transmission clusters. It is an acyclic graph that represents the P2P transmission relationships.  For instance which patient infected whom and graphically view the direction of the transmission. It can be used to view the sampled and non-sampled patients in different colours.

Two txt files should be created for clusterGen. They are as follows
 - The first file contains the nodes
  - The node files contains two columns the Patient ID and the Sampled column.
  - 1 is used to refer to sampled patients and 0 to unsampled or predicted patients.

 -The second text file contains the information of who transmitted to whom
  - The second output file should contain three columns.
  - First column dubbed Primary stating the ID of the source of infection.
  - Second column dubbed Secondary stating the ID of the sink of infection.
  - Third column

_Sample files have been provided._

An instance of the clusterGen is as follows:
>clusterGen('nodes.csv','edges.csv',1)

**The output will be as follows:**
![alt text](https://github.com/theLongLab/Transmission_Analyzer/blob/master/clusterGEN.png)

**Figure 2:** _Novel transmission clusters identified from the transmission tree in Figure 1. Each node represents a sampled patient from the SAC database, and each edge a directed transmission relationship. Coloured nodes represent sampled patients and indicate the corresponding coloured boxes in Figure 1. White nodes represented predicted unsampled patients. Solid edges represent transmission relationships already known to the SAC, whereas dashed edges represent novel inferences of transmission._

[Back To The Top](#table-of-contents)

---

## Procedure
_The entire Procedure has been automated upto the clusterGen in Main_Run._

If you wish to manually conduct the analysis the steps are as follows:
 1. Generate ptree from TransPhylo
 2. brLenChecker
 3. autoSummary
 4. transReport
 5. predictVSKnown
 6. plotCTree _(optional)_
 7. clusterGen _(optional)_  


[Back To The Top](#table-of-contents)

---
## License
**This project is licensed under the MIT License**

MIT License

Copyright (c) 2019 Lauren Mak

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


[Back To The Top](#table-of-contents)

---

## Author Info

- Lauren Mak: Developed and tested the scripts.
- Deshan Perera: Automated and tested the scripts

See also the list of contributors (reference publication to be added) who participated in this project. Didelot et al.: Developed and tested the TransPhylo program.

[Back To The Top](#table-of-contents)

## Citations
_If you find these scripts helpful in your own transmission analyses, please cite (reference publication to be added)._

[Back To The Top](#table-of-contents)

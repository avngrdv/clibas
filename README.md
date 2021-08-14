# FastqProcessor: an NGS data processor tailored for mRNA display selection workflows

FastqProcessor is python package to analyze NGS sequencing data. The tool is specifically
tailored toward analysis of screening outcomes for the selections/screens of genetically
encoded libraries (i.e. mRNA display, phage, yeast, DNa encoded libraries etc).  \
\
The primary purpose of FastqProcessor is library design matching. Given a particular design
for a genetically encoded library, parse .fastq files to remove noise, low confidence reads,
sequences that are too short, too long and much more.\
\
Concominantly, some very basic statistics about the data can be collected as well.\
\
# Dependencies
The code was written and tested for:
	python 3.8.5 
	numpy  1.19.5
	pandas 1.2.4
	matplotlib 3.3.2
\
# Usage
The package provides tools for the analysis, and the analysis pipelines should to be built
for individual applications. code/config.py and code/main.py provide examples for a simple
workflow.

#Documentation
##LibraryDesign
<todo>
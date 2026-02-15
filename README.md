Scripts used in data processing, data analyses, and respective source files, for the article:


**"Genomic and physiological signatures of adaptation in pathogenic fungi"**  
  

Marco Alexandre Guerreiro<sup>1,2,*</sup>, Andrey Yurkov<sup>3</sup>, Minou Nowrousian<sup>4</sup>, Kirk Broders<sup>5</sup>, Eva H. Stukenbrock<sup>1,2</sup>

<sup>1</sup> Environmental Genomics Group, Botanical Institute, Christian-Albrechts University of Kiel, Kiel, Germany  
<sup>2</sup> Max Planck Institute for Evolutionary Biology, Plön, Germany  
<sup>3</sup> Leibniz Institute DSMZ-German Collection of Microorganisms and Cell Cultures, Braunschweig, Germany  
<sup>4</sup> Department of Molecular and Cellular Botany, Ruhr University Bochum, Bochum, Germany  
<sup>5</sup> USDA, Agricultural Research Service, National Center for Agricultural Utilization Research, Mycotoxin Prevention and Applied Microbiology Research Unit, 1815 N. University, Peoria, IL. 61604, U.S.A  

<sup>*</sup> **Corresponding author:** Marco Alexandre Guerreiro,  
Max Planck Institute for Evolutionary Biology,  
August-Thienemann-Str. 2, 24306 Plön  
and  
Christian-Albrechts University of Kiel,  
Am Botanischen Garten 1-9,  
24118 Kiel

Phone: +49 (0) 431 880 6366, Fax: +49 (0) 431 880 6369  
Email: mguerreiro[at]evolbio[dot]mpg[dot]de; mguerreiro[at]bot[dot]uni-kiel[dot]de


Guerreiro MA, Yurkov A, Nowrousian M, Broders K, Stukenbrock EH. Genomic and physiological signatures of adaptation in pathogenic fungi. Nature Communications 17, 748 (2026). https://doi.org/10.1038/s41467-026-68330-6



# Genome_processing
This pipeline performs genome annotation and downstream codon optimization analysis for fungal genomes. It integrates gene prediction (Funannotate), functional annotation (EggNOG, CAZymes) and codon usage analysis (codonR/tAI).

## Quick Start

### 1. Files and directories  
- Provide the path for the directories
- Provide the genome assembly (with .fna extension) in the "Genomes" directory.
- Adjust the paths as needed for the tools. Some tools might change the output format between versions.

### 2. Set the target species  
Edit the pipeline and set the "SPECIES" variable to match the genome assembly file name (without .fna).

### 3. Enable modules  
Toggle specific pipeline steps by setting these flags to Yes or No.  
For a first-time run, it is recommended to do the complete pipeline to generate required intermediate files. Some output files are reformatted for follow-up analyses.

#### Main analyses
- Gene prediction (Funannotate)
- Function annotation (eggNOG)
- Genome completeness (BUSCO)
- Secretome prediction (SignalP)
- tRNA prediction (tRNAscan-SE)    
- CAZymes prediction (dbCAN2)
- Effector prediction (EffectorP)
- Genome metrics (e.g. genome size, GC content)
- Codon usage bias
- Translation efficiency (tAI) and codon optimization (S-values)

The file "tRNA_counts_template.txt" is required for some of the calculation related to RSCU, tAI and S. It is used to reformat output files.


# Data_analyses.r
This scripts loads pre-formatted data from the *Source_files* directory. The raw data is as originated from the *Genome_processing* pipeline. Some tables were simplified or reformatted manually. This script will produce the figures published in the article and available on the *Figures* directory.

## Source_files
This directory contains all files required to generate all the  Figures.

## Figures
This directory contains all the output figures from the *Data_analyses.r* script. Some figures (raw) were further manually edited (edited) in external image editors (e.g. for layout) and are labeled accordingly.

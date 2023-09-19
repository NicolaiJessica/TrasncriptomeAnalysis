# Analyzing global gene expression patterns across publicly available transcriptomes

By Janina K. von Dahlen<sup>1,2</sup>, Kerstin Schulz<sup>1,3</sup>, Jessica Nicolai<sup>1</sup> and Laura E. Rose<sup>1,3</sup>  

<sup>1</sup> Institute of Population Genetics, Heinrich-Heine University Duesseldorf, Universitaetsstr. 1, 40225 Duesseldorf, Germany.  
<sup>2</sup> iGRAD-Plant Graduate School, Heinrich-Heine University Duesseldorf, Duesseldorf, Germany.  
<sup>3</sup> Ceplas, Cluster of Excellence in Plant Sciences, Heinrich-Heine University Duesseldorf, Duesseldorf, Germany.  

# Introduction
This pipeline provides you with a guide to analyze gene expression in a large number of publicly available transcriptomes.  
  
A long-standing objective in biology is to understand which factors affect gene expression. Taking a meta-analysis approach, the amplitude of expression variation across genes in different species and the underlying factors associated with expression differences can be evaluated. One of the strengths as well as a limitation of such kind of meta-analysis is the fact that it implements a large diversity of transcriptomes. With a large number of datasets created under lab-specific settings (using for example different host genotypes and pathogens), consistent gene expression patterns can be viewed as robust, despite the general "noisiness" of the data.  
  
This GitHub manual aims in describing how to analyze gene expression patterns on a global level across a variety of publicly available transcriptomes. The manual refers to the publication:  
  
von Dahlen, J.K., Schulz, K., Nicolai, J., Rose, L.E. (2023): Global expression patterns of R-genes in tomato and potato. Frontiers in Plant Science.  
  
# Requirements & Setup
The pipeline is meant to run on a UNIX based server with the exception of PRIMER which only runs natively on Windows XP or later.  
You will need to install the following programs in advance:  
* [Trimmomatic](https://github.com/usadellab/Trimmomatic)
* FastQC
* Kallisto
* R version X
* R package
** Sleuth
* PRIMER

# Preparing your data (set)
## Data set
Obtain the transcriptomic data you want to analyze. The data can either be obtained through your own sequencing or downloading the .fastq files from databases like NCBI.
You must have one .fastq file per transcriptome for single-end sequencing data and two files for paired-end. It might be necessary to split your .fastq file for paired-end sequencing data. Detailed instruction can be found e.g. on the NCBI website.

## Trimmomatic - Trimming the transcriptomes
To remove low-quality reads and adapter sequencess from the transcriptomes (.fastq) we used Trimmomatic (Bolger et al., 2014). We used the standard parameters. Sample command:  
  
`java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

### Adapter Sequences
The adapter sequences that should be removed from the sequencing data need to be adjusted based on the sequencing method used. Illumina adapters are included in the Trimmomatic download. For further reference see https://github.com/usadellab/Trimmomatic
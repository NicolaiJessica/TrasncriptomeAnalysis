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
  
# Requirements & setup
The pipeline is meant to run on a UNIX based server with the exception of PRIMER which only runs natively on Windows XP or later.  
You will need to install the following programs in advance:  
* [Trimmomatic](https://github.com/usadellab/Trimmomatic)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Kallisto](https://pachterlab.github.io/kallisto/manual)
* R version X
* R package
** [Sleuth](https://github.com/pachterlab/sleuth)
* [PRIMER](https://www.primer-e.com/)

# Preparing your data (set)
## Data set
Obtain the transcriptomic data you want to analyze. The data can either be obtained through your own sequencing or downloading the .fastq files from databases like NCBI.
You must have one .fastq file per transcriptome for single-end sequencing data and two files for paired-end. It might be necessary to split your .fastq file for paired-end sequencing data. Detailed instruction can be found e.g. on the NCBI website.

## Trimmomatic - trimming the transcriptomes
To remove low-quality reads and adapter sequencess from the transcriptomes (.fastq) we used Trimmomatic (Bolger et al., 2014). We used the standard parameters.  
Example command:  
  
`java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

### Adapter sequences
The adapter sequences that should be removed from the sequencing data need to be adjusted based on the sequencing method used. Illumina adapters are included in the Trimmomatic download. For further reference see https://github.com/usadellab/Trimmomatic

## FastQC - qualilty control
Quality controls of trimmed transcriptomes were performed using FastQC (Andrews, 2010). Transcriptomes with overall low-quality or remaining adapters were removed from study.  
Example command:  
  
`fastqc seqfile1 seqfile2 .. seqfileN`

## Kallisto - calculation of transcript abundances 
The program Kallisto (v.0.46.0) can be used to estimate the relative expression of genes (Bray et al., 2016). As a first step, the raw sequence reads were compared to the transcript sequences. This step in Kallisto is designated as the pseudoalignment step. As Kallisto requires information on fragment length for single-end sequenced transcriptomes, the fragment length denoted by the authors was used. If this information was not available, the recommended fragment length of the reported RNA isolation kit was used. 

### Building an index
Before using Kallisto to calculate transcript abundance an index has to be build. The index is used to identify/map the transcriptomic data. In our study we used the CDS to build the index. You need one index for each species you want to analyze.  
Example command:  
  
`kallisto index [arguments] FASTA-files`

### TMP calculation
Transcript abundance was calculated as transcripts per million (TPM; Wagner et al., 2012). TPM normalizes the transcript abundance for gene length and sequencing depth, making TPM values comparable across transcriptomes. To increase the robustness of this calculation we recommend using the -b argument to run the analysis with bootstrap replicates. 
Example command:  
  
`kallisto quant [arguments] FASTQ-files`
  
## Sleuth - differential expression analysis
Differentially expressed genes between conditions (e.g. roots vs. shoot or pathogen-treated vs non-treated transcriptomes) were identified using the R package Sleuth (Pimentel et al., 2017).  
Example command:  
  
### Specify the directory with your kallisto results:  
`sample_id <- dir(file.path("..", "results"))`
  
After executing `sample_id` your output should look like this:  
  
`\## [1] "SRA1" "SRA2" "SRA3" "SRA4"`

### Add paths for each result:  
  
`kal_dirs <- file.path("..", "results", sample_id, "kallisto")`

The output of `kal_dirs` should look look like this:  
  
`## [1] "../results/SRA1/kallisto" "../results/SRA2/kallisto"`

### Adding the experimental design information:  
You need to supply a table that specifies which treatment each SRA belongs to. 
Example:  
|sample|condition|
|SRA1|Mock|
|SRA2|Mock|
|SRA3|Treatment1|
|SRA4|Treatment2|

| Sample 	| Conditon  	|
|--------	|-----------	|
| SRA1   	| Mock      	|
| SRA2   	| Mock      	|
| SRA3   	| Treatment 	|
| SRA4   	| Treatment 	|
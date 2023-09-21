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
* R version >= 3.6.3
* R package: [Sleuth](https://github.com/pachterlab/sleuth)
* [PRIMER](https://www.primer-e.com/)

# Preparing your data
## Data set
This pipeline can be run by using publicly available transcriptomes as well as self-generated ones. Publicly available transcriptomes can be analyzed downloading the .fastq files from databases like NCBI. 
You must have one .fastq file per transcriptome for single-end sequencing data and two files for paired-end. For paired-end sequencing data that only consists of one single .fastq file containing both the forward and reverse reads you need to split your .fastq file. Detailed instruction for how to download publicly available transcriptomes and how to split them can be found e.g. on the NCBI website.

## Trimmomatic - trimming the transcriptomes
To remove low-quality reads and adapter sequences from the transcriptomes (.fastq) we used Trimmomatic with standard parameters (Bolger et al., 2014).   
Example command:  
  
`java -jar trimmomatic-0.39.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

### Adapter sequences
The adapter sequences that should be removed from the sequencing data need to be adjusted based on the sequencing method used. Illumina adapters are included in the Trimmomatic download. For further reference see https://github.com/usadellab/Trimmomatic

## FastQC - qualilty control
Quality controls of trimmed transcriptomes were performed using FastQC (Andrews, 2010). Transcriptomes with overall low-quality or remaining adapters were removed from study. For further reference see https://github.com/s-andrews/FastQC  
Example command:  
  
`fastqc seqfile1 seqfile2 .. seqfileN`

## Kallisto - calculation of transcript abundances 
The program Kallisto (v.0.46.0) can be used to estimate the relative expression of genes (Bray et al., 2016). As a first step, the raw sequence reads were compared to the transcript sequences. This step in Kallisto is designated as the pseudoalignment step. As Kallisto requires information on fragment length for single-end sequenced transcriptomes, the fragment length denoted by the authors was used. If this information was not available, the recommended fragment length of the reported RNA isolation kit was used. For further reference see https://github.com/pachterlab/kallisto

### Building an index
Before using Kallisto to calculate transcript abundance an index has to be build. The index is used to pseudoalign the raw transcriptomic data to the transcript sequences. In our study we used the CDS to build the index. One index for each species you want to analyze is necessary.  
Example command:  
  
`kallisto index [arguments] FASTA-files`

### TPM calculation
Transcript abundance was calculated as transcripts per million (TPM; Wagner et al., 2012). TPM normalizes the transcript abundance for gene length and sequencing depth, making TPM values comparable across transcriptomes. To increase the robustness of this calculation we recommend using the -b argument to run the analysis with bootstrap replicates. 
Example command:  
  
`kallisto quant [arguments] FASTQ-files`
  
## Sleuth - differential expression analysis
Differentially expressed genes between conditions (e.g. roots vs. shoot or pathogen-treated vs non-treated transcriptomes) were identified using the R package Sleuth (Pimentel et al., 2017).  
The p-values need to be adjusted using the Benjamini-Hochberg correction (FDR ≤ 0.05; Benjamini & Hochberg, 1995). Since Sleuth relies on replicates within treatments, transcriptomes without replicates need to be removed from the analysis.   
  
### Loading sleuth
  
`library("sleuth")`  
  
### Specify the directory with your kallisto results  
The input for sleuth are the .tsv and .h5 file that Kallisto generates. The results for each transcriptome must be in their own directory and the directory should be named after the SRA.  
`sample_id <- dir(file.path("..", "results"))`
  
After executing `sample_id` your output should look like this:  
  
`## [1] "SRA1" "SRA2" "SRA3" "SRA4" "SRA5" "SRA6`

### Add paths for each result  
  
`kal_dirs <- file.path("..", "results", sample_id, "kallisto")`

The output of `kal_dirs` should look like this:  
  
`## [1] "../results/SRA1/kallisto" "../results/SRA2/kallisto"`

### Adding the experimental design information  
In the next step, you need to supply a table that specifies which treatment each SRA belongs to. 
Example table:  

| **Sample** 	| **Condition**   |
|--------	    |-----------	    |
| SRA1   	    | Mock           	|
| SRA2   	    | Mock          	|
| SRA3   	    | Mock 	          |
| SRA4   	    | Treatment     	|
| SRA5   	    | Treatment     	|
| SRA6   	    | Treatment     	|
  
```
s2c <- read.table(file.path(".."), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
```  

### Appending kal_dirs to table  

`s2c <- dplyr::mutate(s2c, path = kal_dirs)`  
  
The output from `print(s2c)` should look like this:  
  
| **Sample** 	| **Condition** | **Path**                  |
|------------	|--------------	|--------------------------	|
| SRA1       	| Mock         	| ../results/SRA1/kallisto 	|
| SRA2       	| Mock         	| ../results/SRA2/kallisto 	|
| SRA3       	| Mock        	| ../results/SRA3/kallisto 	|
| SRA4       	| Treatment    	| ../results/SRA4/kallisto 	|
| SRA5       	| Treatment    	| ../results/SRA5/kallisto 	|
| SRA6       	| Treatment    	| ../results/SRA6/kallisto 	|  
  
### Optional - Adding gene annotation

`t2g <- read.table("..", header = TRUE, stringsAsFactors = FALSE)`

### Construction of "sleuth object" and fitting a model
  
```
so <- sleuth_prep(s2c, ~condition, target_mapping = t2g, read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE, transformation_function = function(x) log2(x + 0.5))
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
```  
  
If you did not provide an annotation file, remove "target_mapping = t2g" from the first command of this step.  

### Likelihood ratio test

```
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 
sleuth_significant_lrt <- dplyr::filter(sleuth_table_lrt, qval <= 0.05) 
```  
  
### Wald test / foldchange calculation / significant expression differences calculation
First you need to specify for which condition you want to perform the calculation.  
The available conditions can be accessed by executing `models(so)`.  
  
```
so <- sleuth_wt(so, condition, which_model = "full")  
results_Wald_table <- sleuth_results(so, condition, test_type = "wt")  
sleuth_Wald_significant <- dplyr::filter(results_Wald_table, qval <= 0.05)
```  
  
## PRIMER - multivariate analysis of expression differences
Gene expression levels can be affected among other factors by the tissue, the genotype or the treatment. To identify the factors associated with differences in expression of genes across transcriptomes, we performed an ANOSIM in Primer 7.0.13 (PRIMER-e). ANOSIM is a non-parametric statistical test similar to ANOVA. As sample data, an excel sheet table with TPM-values was used (`File > open > choose input data > sample data` (no titles, data type: unknown/other). Missing values were treated as blanks. The sample data should have the following format:  
  
|               	| **Transciptome 1** 	| **Transciptome 2** 	| **Transciptome 3** 	| **Transciptome 4** 	| **Transciptome 5** 	| **...** 	|
|---------------	|--------------------	|--------------------	|--------------------	|--------------------	|--------------------	|---------	|
| **Gene-ID 1** 	| TPM                	| TPM                	| TPM                	| TPM                	| TPM                	| TPM     	|
| **Gene-ID 2** 	| TPM                	| TPM                	| TPM                	| TPM                	| TPM                	| TPM     	|
| **Gene-ID 3** 	| TPM                	| TPM                	| TPM                	| TPM                	| TPM                	| TPM     	|
| **Gene-ID 4** 	| TPM                	| TPM                	| TPM                	| TPM                	| TPM                	| TPM     	|
| **Gene-ID 5** 	| TPM                	| TPM                	| TPM                	| TPM                	| TPM                	| TPM     	|
| **Gene-ID 6** 	| TPM                	| TPM                	| TPM                	| TPM                	| TPM                	| TPM     	|
| **Gene-ID 7** 	| TPM                	| TPM                	| TPM                	| TPM                	| TPM                	| TPM     	|
| **Gene-ID 8** 	| TPM                	| TPM                	| TPM                	| TPM                	| TPM                	| TPM     	|
| **...**       	| TPM                	| TPM                	| TPM                	| TPM                	| TPM                	| TPM     	|  
  
To add factors to the data set, choose: edit>factors. The factor data should have the following format:  
  
|                     	| **Factor 1 -  Tissue** 	| **Factor 2 - Treatment** 	| **Factor 3 - Genotype** 	| **...** 	|
|---------------------	|------------------------	|--------------------------	|-------------------------	|---------	|
| **Transcriptome 1** 	| Root                   	| Mock                     	| AB                      	| ...     	|
| **Transcriptome 1** 	| Root                   	| Treatment                	| AB                      	| ...     	|
| **Transcriptome 1** 	| Shoot                  	| Mock                     	| AB                      	| ...     	|
| **Transcriptome 1** 	| Shoot                  	| Treatment                	| BB                      	| ...     	|
| **Transcriptome 1** 	| Root                   	| Treatment                	| BB                      	| ...     	|
| **...**             	| ...                    	| ...                      	| ...                     	| ...     	|  
  
The starting point of the analysis is a pairwise dissimilarity matrix. In our case, the dissimilarity matrix was computed as follows: First the TPM values for each gene within each transcriptome were LOG (x+1) transformed (Pre-treatment > transform overall > LOG(x+1)) and standardized across libraries (Pre-treatment > standardized: Standardizes samples by total). On the basis of these transformed TPM values, the dissimilarity in gene expression patterns between transcriptomes were calculated based on Euclidean distances (`Analyze > resemblance > Euclidean distance > analyze between samples`). Ranking was applied to the distance matrix.  
  
### ANOSIM
To determine if gene expression is more similar within groups than between groups (for example when groups are defined by infection status or tissue type) the R test statistic value was calculated. The R values can range from -1 to 1, with larger values corresponding to greater differences between groups (Analyze > ANOSIM > Model: one-way-A > choose factor). Statistical significance is calculated through 999x permutations of the group labels and recalculation of the R value for each replicate. 
  
#### R value categories:  
  
| **R value**        	| **Meaning**                                    	|
|--------------------	|-------------------------------------------------	|
| **0.75 < R < 1**   	| highly different                                	|
| **0.5 < R < 0.75** 	| different                                       	|
| **0.25 < R < 0.5** 	| different with some overlap                     	|
| **0.1 < R < 0.25** 	| similar with some differences (or high overlap) 	|
| **R < 0.1**        	| similar                                         	|  
  
### PCA
For visualizing data as PCA pick the standardize dataset and choose: `analyze > PCA`.  
    
## Statistics
Next to ANOSIM, other statistic tests such as the Mann-Whitney-U test (Mann & Whitney, 1947) or the Spearman's rank correlation (Hollander et al., 2013) can be performed. The following statistic tests were performed in R.  

The data should have to following format:
  
| **Data row 1**      | **Data row 2**                    |
|--------------------	|-----------------------------------|
| TPM / DEG    	      | TPM / DEG                         |
| TPM / DEG          	| TPM / DEG                         |
| TPM / DEG	          | TPM / DEG                        	|
| TPM / DEG 	        | TPM / DEG                        	|
| ...               	| ...                             	|  

### Shapiro test - testing for normal distribution
For testing the data for normal distribution (<5000 data points per column), the Shapiro test can be used (Shapiro & Wilk, 1965):
  
`shapiro.test(file name[,x])`
  
Comments: 
* The [,x] is referring to the text column which should be analyzed
* p <0.05: data is non-normal distributed
* p ≥0.05: data is normal distributed

For testing >5000 data points per column for normal distribution, the ad.test within the Nortest package can be used.  
  
### Testing for equal variance  
For testing the data for equal variance, the var.test can be used:  
  
`var.test(file name[,x], file name[,y], alternative="two.sided")`
  
Comments:
* Data can be either tested for “one.sided” or “two.sided”
* p-value ≥ 0.05: equal variance
* p-value < 0.05: unequal variance
  
### Mann-Whitney-U test for non-parametric data
Non-parametric data sets can be analyzed using the Mann-Whitney-U test:

`wilcox.test(file name[,x],file name[,y],alternative="two.sided")`

Comments:
* p-value ≥ 0.05: accept null hypothesis 
* p-value < 0.05: reject null hypothesis (significant difference)

### Student's t-test for parametric data
Parametric data sets can be analyzed using a student’s t-test:

`t.test(file name[,x],file name[,y],alternative="two.sided",var.equal=TRUE)`

Comments:
* in case the variance is not equal choose option: var.equal=FALSE
* p-value ≥ 0.05: accept null hypothesis 
* p-value < 0.05: reject null hypothesis (significant difference)
  
### Spearman's rank correlation
Spearman's rank correlations (Hollander et al., 2013) are a non-parametric measure of statistical dependency between two variables. To perform a Spearman's rank correlation:

`cor.test(file name[,x],file name[,y],alternative="greater",method=c("spearman"), exact =NULL, continuity = TRUE)`

Comments:
* in case there is a negative association use: alternative =”less”
* Exact: a logical indicating whether an exact p-value should be computed or not

To estimate the strength of the association the Rho value can be used:
Rho 0.00-0.19: very weak association
Rho 0.20-0.39: weak association
Rho 0.40-0,59: moderate association
Rho 0.60-0.79: strong association
Rho 0.80-1.00: very strong association

# Test files
If you want to test this pipeline, you may use the following transcriptomes (available on NCBI) and compare your results with our paper.
* SRR7734429
* SRR7734430
* SRR7734431
* SRR7734432
* SRR7734435
* SRR7734436  
  
# How to cite us?

von Dahlen, J.K., Schulz, K., Nicolai, J., Rose, L.E. (2023): Global expression patterns of R-genes in tomato and potato. Frontiers in Plant Science.

# References
* Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.  
* Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc  
* Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature biotechnology, 34(5), 525.  
* Wagner, G. P., Kin, K., & Lynch, V. J. (2012). Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. Theory in biosciences, 131(4), 281-285.  
* Pimentel, H., Bray, N. L., Puente, S., Melsted, P., & Pachter, L. (2017). Differential analysis of RNA-seq incorporating quantification uncertainty. Nature methods, 14(7), 687.  
* Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal statistical society: series B (Methodological), 57(1), 289-300.  
* Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). Biometrika, 52(3/4), 591-611.  
* Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. The annals of mathematical statistics, 50-60.  
* Hollander, M., Wolfe, D. A., & Chicken, E. (2013). Nonparametric statistical methods (Vol. 751). John Wiley & Sons.  

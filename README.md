06.tpm.py
# MT-DeepSea
MT-DeepSea: an optimized strategy for deep-sea metatranscriptomic analysis
## Dependencies
- Python >3.6
- Nextflow >22.10
- Trimmomatic 0.39
- Fastqc 0.11.9
- BBMap 38.87
- samtools 1.11
- Bowtie2 2.4.2
- featureCounts 2.0.1
- Kraken2 2.1.2

##  Installation
Please execute the following command in your terminal to clone the metatranscriptome repository to your own machine. 
```shell
git clone https://github.com/Li-Lab-SJTU/MT-DeepSea.git
```
## Input description 
The input data should contain the following parts:
Input  | Description
------------- | -------------
metatranscriptomics | reads[1\|2].fastq.gz
rRNA database index  | bowtie2 index for rRNA database
trimmomatic | the absolute path of trimmomatic.jar
adapter | the adapter sequences used for filtration
kraken2 database | database built by kraken2 for taxonomy analysis
reference genome sequence | .fasta
reference genome annotation | .gtf

## Implementation
You can implement the code as following:
```shell
nextflow run MT-DeepSea.nf --stdin1 <reads1.fastq.gz> --stdin2 <reads2.fastq.gz> --rrna_db <rrnadb_index/rrna_db> --trimmomatic <path/to/trimmomatic.jar> --adapter <adapter.fa> --kraken2_db <kraken2_db> --ref <genome.fasta> --gtf <genome.gtf> <Options> 
```
**Note:** This step will take about **two hours**. The input data should be stored in the `data` folder. The output results will be stored in the folder named `results` by default.


#### Options arguments
|parameters|descriptions|
|---|---|
|--sample|Specify the prefix of output files. Default: the prefix of the name of stdin1 (exclude suffix).|
|--outdir|Specify an output folder. Default: the results folder in the current path where the command is executed.|
|--thread| Set the thread. Default: 40.|
|--help|Print this help message.|

# MT-DeepSea
MT-DeepSea: an optimized strategy for deep-sea metatranscriptomic analysis
## Dependencies
- Python >3.6
- Nextflow 22.10.4.5836
- Trimmomatic 0.39
- Fastqc 0.11.9
- BBMap 38.87
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
reference genome sequence | .fasta
reference genome annotation | .gtf

## Implementation
You can implement the code as following:
```shell
nextflow run MT-DeepSea.nf --stdin1 <reads1.fastq.gz> --stdin2 <reads2.fastq.gz> --rrna_db <rrnadb_index/rrna_db> --ref <genome.fasta> --gtf <genome.gtf> --outdir <result_fold> <Options> <Functions> 
```
**Note:** This step will take about **two hours**. The input data should be stored in the `data` folder. The output results will be stored in the folder named `results` by default.


#### Options arguments
|parameters|descriptions|
|---|---|
|--sample|Specify the prefix of output files. Default: the prefix of the name of stdin1 (exclude suffix).|
|--outdir|Specify an output folder. Default: the results folder in the current path where the command is executed.|
|--thread| Set the thread. Default: 40.|
|--cache| The mode to store the process results to a local cache, see [nextflow document](https://www.nextflow.io/docs/latest/basic.html) for details. Default: "deep".|
|--publish_mode|The publication mode of output files, see [nextflow document](https://www.nextflow.io/docs/latest/basic.html) for details. Default: "copy".|
|--help|Print this help message.|
#### Functions arguments
These parameters are built-in functions of Nextflow, they can generate some visual graphics, which or show the total time consumption of the pipeline, or show the time consumption, memory occupation, cpu usage of each process. Interested can add these parameters to observe relative information. See [nextflow document](https://www.nextflow.io/docs/latest/basic.html) for details.
|parameters|descriptions|
|---|---|
|-with-timeline|It renders a timeline.html file that records the time, memory consumption of different processes.|
|-with-report|It generates a report.html file that records the single core CPU Usage, execution time, memory occupation and Disk read write information of different processes.|
|-with-trace|It creates an execution tracing file that contains some useful information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used.|
|-with-dag|It outputs the pipeline execution DAG. It creates a file named dag.dot containing a textual representation of the pipeline execution graph in the DOT format.|
|-resume|It means only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead. Also, the pipeline can be restarted by add the parameter when any disconnection of the network or server occurs.|

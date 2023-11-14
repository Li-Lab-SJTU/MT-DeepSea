#! /usr/bin/env nextflow

/*
=====================================================================================
                                MT-Enviro.nf                                      
=====================================================================================
       An optimized strategy for environmental metatranscriptomic analysis                                     
=====================================================================================
*/

//

params.help = false
if (params.help) {
    HelpMessage()
    exit 0
}

def HelpMessage() {
    log.info """
MT-Enviro by Weiyi Li (oliva630@sjtu.edu.cn)
Please run the command: [ nextflow run MT-Enviro.nf --help ] to print the help message.\n
Usage: 
    
    nextflow run MT-Enviro.nf --stdin1 <reads1.fastq.gz> --stdin2 <reads2.fastq.gz> --rrna_db <rrnadb_index/rrna_db> --adapter <adapter.fa> --kraken2_db <kraken2_db> --ref <genome.fasta> --gtf <genome.gtf> <Options>
        
Input:

    --stdin1         Fastq files, gzip'ed (extension: .gz).

    --stdin2         For paired reads in two files, gzip'ed (extension: .gz).

    --rrna_db        The rRNA database for the rRNA removal. Note that this database needs to be the index files for bowtie2. See bowtie2-build for details.

    --adapter        The adapter sequences used for filtration (.fa).

    --kraken2_db     The kraken2 database for kraken2. Note that this parameter needs a folder that contains kraken2 database. See kraken2 -db for details.

    --ref            Specify a reference genome for reads mapping, fasta file.

    --gtf            Specify a GTF file as the annotation of the reference genome.
 
Options:
    
    --sample        Specify the prefix of output files. Default: the prefix of the name of stdin1 (exclude suffix).

    --outdir        Specify an output folder. Default: the results folder in the current path where the command is executed.

    --thread        Set the thread. Default: 40.

    --help          Print this help message.

    """.stripIndent()
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    SET UP CONFIGURATION VARIABLES                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////


params.sample = false
params.stdin1 = null
params.stdin2 = null
params.outdir = "./results"

params.adapter = null

params.rrna_db = null
params.ref = null 
params.gtf = null
params.kraken2_db = null

params.thread = 40
params.cache = "deep"
params.publish_mode = "copy"

if (params.stdin1 == null && params.stdin2 == null && params.adapter == null && params.rrna_db == null && params.ref == null && params.gtf == null && params.kraken2_db == null) {
    HelpMessage()
    exit 0
}
//When the command line lacks the corresponding parameter, give the necessary warning
def Warnings() {
    if (params.stdin1 == null) {
        HelpMessage()
        println "\nError:\n     The standard input file 1 is not found, please specify it by adding the --stdin1 parameter in the command-line!\n"
        exit 0
    }
    if (params.stdin2 == null) {
        HelpMessage()
        println "\nError:\n     The standard input file 2 is not found, please specify it by adding the --stdin2 parameter in the command-line!\n"
        exit 0
    }
    if (params.adapter == null) {
        HelpMessage()
        println "\nError:\n     The adapter sequences are not found, please specify it by adding the --adapter parameter in the command-line!\n"
        exit 0
    }
    if (params.rrna_db == null) {
        HelpMessage()
        println "\nError:\n     The rRNA database is not found, please specify it by adding the --rrna_db parameter in the command-line!\n"
        exit 0
    }
    if (params.ref == null) {
        HelpMessage()
        println "\nError:\n     The reference database is not found, please specify it by adding the --ref parameter in the command-line!\n"
        exit 0
    }
    if (params.gtf == null) {
        HelpMessage()
        println "\nError:\n     The GTF file is not found, please specify it by adding the --gtf parameter in the command-line!\n"
        exit 0
    }
    if (params.kraken2_db == null) {
        HelpMessage()
        println "\nError:\n     The database for kraken2 is not found, please specify it by adding the --kraken2_db parameter in the command-line!\n"
        exit 0
    }
}
Warnings()

def Setsample() {
    if (params.sample == false) {
        infile = file("$params.stdin1")
        sample = infile.simpleName
    }
    else if (params.sample) {
        sample = params.sample
    }
}
Setsample()

def Makedirs() {
    outdir = file("$params.outdir")
    if (outdir.exists() == false) {
        outdir.mkdir()
    }

    qc_outdir = file("${params.outdir}/01_readTrimming")
    qc_outdir.mkdir()

    fastqc_outdir = file("${params.outdir}/02_qualityEvaluation")
    fastqc_outdir.mkdir()
    
    rrnaRemove_outdir = file("${params.outdir}/03_Deep-rRNA")
    rrnaRemove_outdir.mkdir()
    
    readsMapping_outdir = file("${params.outdir}/04_readAlignment")
    readsMapping_outdir.mkdir()
    
    featureCounts_outdir = file("${params.outdir}/05_readSummarization")
    featureCounts_outdir.mkdir()

    kraken2_outdir = file("${params.outdir}/06_taxonomyAnnotation")
    kraken2_outdir.mkdir()
}
Makedirs()

Channel.fromPath(params.stdin1, checkIfExists: true)
       .ifEmpty{exit 1, "The $params.stdin1 file is empty!"}
       .view{"$workflow.start - INFO - Load standard input file 1: $it"}
       .set{stdin1}

Channel.fromPath(params.stdin2, checkIfExists: true)
       .ifEmpty{exit 1, "The $params.stdin2 file is empty!"}
       .view{"$workflow.start - INFO - Load standard input file 2: $it"}
       .set{stdin2}


// read trimming for raw squencing data. Tools: Trimmomatic

process readTrimming {

    cache params.cache
    publishDir "$qc_outdir", mode: params.publish_mode

    input:
    path r1 from stdin1
    path r2 from stdin2

    output:
    path "*"
    path "${sample}_01_trim_1.fastq.gz" into trim_fq1
    path "${sample}_01_trim_2.fastq.gz" into trim_fq2

    script:
    """
    c1=${sample}_01_trim_1.fastq.gz
    c2=${sample}_01_trim_2.fastq.gz
    t1=${sample}_01_trim_unpaired_1.fastq.gz
    t2=${sample}_01_trim_unpaired_2.fastq.gz

    java -jar ${workflow.projectDir}/bin/trimmomatic-0.39.jar PE -phred33 $r1 $r2 \${c1} \${t1} \${c2} \${t2} ILLUMINACLIP:${params.adapter}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 HEADCROP:10
    """
}

// copy trimmed fastq file into multiple channels
trim_fq1.into{trim_fq1_01; trim_fq1_02}
trim_fq2.into{trim_fq2_01; trim_fq2_02}

// quality evaluation of trimmed fastq file
// tools: fastqc

process qualityEvaluate {

    cache params.cache
    publishDir "$fastqc_outdir", mode: params.publish_mode

    input:
    path c1 from trim_fq1_01
    path c2 from trim_fq2_01

    output:
    path "*"

    script:
    """
    fastqc -f fastq $c1 -o ./ -t $params.thread 
    fastqc -f fastq $c2 -o ./ -t $params.thread
    """
}

// rrna removal of trimmed meta-transcriptomics fastq file
// tools: bowtie2

process DeeprRNA {

    cache params.cache
    publishDir "$rrnaRemove_outdir", mode: params.publish_mode

    input:
    path c1 from trim_fq1_02
    path c2 from trim_fq2_02

    output:
    path "*"
    path "${sample}_filter_1.fastq.gz" into filter_fq1
    path "${sample}_filter_2.fastq.gz" into filter_fq2

    script:
    """
    bowtie2 -x $params.rrna_db -1 $c1 -2 $c2 --no-head --no-unal --no-mixed -p $params.thread --un-conc-gz ${sample}_filter.fastq.gz --al-conc-gz ${sample}_rRNA.fastq.gz -S ${sample}_rRNA.sam

    rm ${sample}_rRNA.sam

    mv ${sample}_filter.fastq.1.gz ${sample}_filter_1.fastq.gz
    mv ${sample}_filter.fastq.2.gz ${sample}_filter_2.fastq.gz
    """
}

// copy trimmed and rrna-removed fastq file into multiple channels
filter_fq1.into{filter_fq1_01; filter_fq1_02}
filter_fq2.into{filter_fq2_01; filter_fq2_02}

// reads alignment
// tools: script - bbmap.sh & samtools

process readAlignment {

    cache params.cache
    publishDir "$readsMapping_outdir", mode: params.publish_mode

    input:
    path f1 from filter_fq1_01
    path f2 from filter_fq2_01

    output:
    path "*"
    path "${sample}_04_bbmap.sorted.bam" into mapped_bam

    script:
    """
    bbmap.sh in=$f1 in2=$f2 ref=$params.ref nodisk threads=$params.thread covstats=${sample}_depth_bbmap.depth out=${sample}_04_bbmap_temp.sam

    samtools view -bS -S ${sample}_04_bbmap_temp.sam -o ${sample}_04_bbmap.bam
    samtools sort -@ $params.thread -O bam ${sample}_04_bbmap.bam -o ${sample}_04_bbmap.sorted.bam

    rm ${sample}_04_bbmap_temp.sam
    """
}

// calculate and obtain gene count matrix
// tools: featureCounts

process readSummarization {
    
    cache params.cache
    publishDir "$featureCounts_outdir", mode: params.publish_mode

    input:
    path bam from mapped_bam

    output:
    path "*"
    
    project_dir = projectDir

    script:
    """
    featureCounts -a $params.gtf -p $bam -t transcript -g gene_id -T $params.thread -o ${sample}_bp_trim.txt
    python ${workflow.projectDir}/bin/06_tpm.py ${sample}_bp_trim.txt ${sample}_bp_trim_tpm.txt
    """
}

process taxonomyAnnotation {

    cache params.cache
    publishDir "$kraken2_outdir", mode: params.publish_mode

    input:
    path f1 from filter_fq1_02
    path f2 from filter_fq2_02
    
    output:
    path "*"

    script:
    """
    kraken2 --db $params.kraken2_db --output ${sample}_06_kraken2_nr.out --report ${sample}_06_kraken2_nr.report --paired --gzip-compressed $f1 $f2 --threads $params.thread
    """
}

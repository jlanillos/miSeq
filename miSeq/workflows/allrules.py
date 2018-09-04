#!python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#configfile: "/Users/javierlanillos/Desktop/Neoantigens/Project/Project_RNA/src/config_files/config_template_javier.json"


#source activate py36
import json
# Directories:
#fastqc_dir = config["analysis"]["analysis_dir"] + config["analysis"]["sample_id"]  + "/" + config["analysis"]["result"]  + "fastqc/"
#rule_dir = "/Users/javierlanillos/Desktop/Neoantigens/Project/Project_RNA/"

cutadapt_dir = config["analysis"]["analysis_dir"]  + config["analysis"]["result"]  + "cutadapt/"
fastqc_dir = config["analysis"]["analysis_dir"]  + config["analysis"]["result"]  + "fastqc/"
bam_dir = config["analysis"]["analysis_dir"]  + config["analysis"]["result"]  + "bam/"
result_dir = config["analysis"]["analysis_dir"] + config["analysis"]["result"]
cufflinks_dir = config["analysis"]["analysis_dir"]  + config["analysis"]["result"] + "cuffout/"
rule_dir = "/Users/javierlanillos/Desktop/Neoantigens/Project/Project_RNA/"
vcf_dir = config["analysis"]["analysis_dir"] + config["analysis"]["result"]  + "vcf/"
vep_dir =  config["analysis"]["analysis_dir"] + config["analysis"]["result"]  + "vep/"




######## RULE_ALL #######

rule all:
  input:
#    expand(fastqc_dir + "{sample}_1_fastqc.zip", sample = config["samples"]),
#    expand(fastqc_dir + "{sample}_2_fastqc.zip", sample = config["samples"]),
    expand(cutadapt_dir + "{sample}_1.ca.fastq.gz", sample = config["samples"]),
    expand(cutadapt_dir + "{sample}_2.ca.fastq.gz", sample = config["samples"]),
    expand(cutadapt_dir + "{sample}_1.ca.short.fastq.gz", sample = config["samples"]),
    expand(cutadapt_dir + "{sample}_2.ca.short.fastq.gz", sample = config["samples"]),
    expand(bam_dir + "{sample}.bam", sample=config["samples"]),
    expand(bam_dir + "{sample}.sorted.bam", sample=config["samples"]),
    expand(bam_dir + "{sample}.sorted.mrkdup.bam", sample = config["samples"]),
#### Add RealignerTargetCreator and IndelRealignment (so, add GATK3.8) when using other variant callers. For Mutect2 is not longer necessary with GATK4
####expand(bam_dir + "{sample}" + ".sorted.mrkdup.bam" + ".intervals", sample = config["samples"]),
####expand(bam_dir + "{sample}" + ".sorted.mrkdup.ralgn.bam", sample = config["samples"]),
####expand(bam_dir + "{sample}" + ".sorted.mrkdup.ralgn.bsrcl.bam", sample = config["samples"]),
####expand(bam_dir + "{sample}" + ".sorted.mrkdup.ralgn.bsrcl.bam.list", sample = config["samples"]),
    expand(result_dir + "{sample}" + ".sorted.mrkdup.bam.picard.bedintervals", sample = config["samples"]),
    expand(result_dir + "{sample}.sorted.mrkdup.hsmetric", sample = config["samples"]),
    expand(result_dir + "{sample}.sorted.alignmetric", sample = config["samples"]),
    bam_dir + "tumor.bam",
    expand(vcf_dir + config["vcf"]["mutect"]["default"]),
    vcf_dir + config["vcf"]["mutect"]["default"],
    vcf_dir + "filt_split_" + config["vcf"]["mutect"]["default"],
    vep_dir + "vep_" + config["vcf"]["mutect"]["default"],
###    expand(bam_dir + "normal_tumor.contest", sample=config["samples"]),
###    expand(bam_dir + "tumor_normal.contest", sample=config["samples"]),

######### FASTQC ###############
# Get summary statistics from sequencing data
rule fastq:
  input:
    read1= config["path"]["fastq"] + "{sample}" + "1_001.fastq.gz",
    read2= config["path"]["fastq"] + "{sample}" + "2_001.fastq.gz",
  output:
    read1 = fastqc_dir + "{sample}_1_fastqc.zip",
    read2 = fastqc_dir + "{sample}_2_fastqc.zip",
  params:
    fastqc_dir=fastqc_dir
  shell:
    "fastqc {input.read1} --outdir {params.fastqc_dir};"
    "fastqc {input.read2} --outdir {params.fastqc_dir};"


########## CUTADAPT ############
# remove adapter from paired end reads
rule cutadapt:
  input:
    read1= config["path"]["fastq"] + "{sample}" + "1_001.fastq.gz",
    read2= config["path"]["fastq"] + "{sample}" + "2_001.fastq.gz",
  output:
    read1 = cutadapt_dir + "{sample}_1.ca.fastq.gz",
    read2 = cutadapt_dir + "{sample}_2.ca.fastq.gz",
    read1_short = cutadapt_dir + "{sample}_1.ca.short.fastq.gz",
    read2_short = cutadapt_dir + "{sample}_2.ca.short.fastq.gz"
  params:
    adapter = config["QC"]["adapter"],
    min_seq_length = config["QC"]["min_seq_length"],
  threads: 4
  shell:
    "cutadapt "
            "-b {params.adapter} -B {params.adapter} "
            "--minimum-length {params.min_seq_length} "
            "--too-short-output {output.read1_short} "
            "--too-short-paired-output {output.read2_short} "
            "-o {output.read1} "
            "-p {output.read2} "
            "{input.read1} {input.read2}; "

########## BWA MEM ############

# Take cutadapt (trimmed) output fastq files and align with bwa mem. Also, convert the output to sam format
rule bwa_mem:
  input:
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
    #fa = config["path"]["genomefa"] + config["references"]["genomefa"],
    read1 = cutadapt_dir + "{sample}_1.ca.fastq.gz",
    read2 = cutadapt_dir + "{sample}_2.ca.fastq.gz",
    #refidx = expand(config["path"]["genomefa"] + config["references"]["genomefa"] + ".{prefix}", prefix=["amb","ann","bwt","pac","sa"])
    refidx = expand(config["path"]["genomefa"] + config["references"]["genomefa_GATK"] + ".{prefix}", prefix=["amb","ann","bwt","pac","sa"]),
  output:
    bamout = bam_dir + "{sample}.bam",
  params:
    header_1 = "'@RG\\tID:" +  "{sample}" + "\\tSM:" + "{sample}" + "\\tPL:ILLUMINAi'",
  threads: 4
  shell:
    "bwa mem "
        "-t {threads} "
        "-R  {params.header_1} "
        "-M "
        "-v 1 "
        "{input.fa} {input.read1} {input.read2} "
        "| samtools view -Sb - > {output.bamout}; "


######## SAMTOOLS_SORT #########
rule samtools_sort_index:
    input:
        bam_dir + "{sample}.bam",
    output:
        bam_dir + "{sample}.sorted.bam",
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output} {input}"


########## PICARD TOOLS MARK DUPLICATES ############

rule picard_mark_duplicate:
  input:
    bam_dir + "{sample}" + ".sorted.bam",
  output:
    mrkdup = bam_dir + "{sample}" + ".sorted.mrkdup.bam",
  log:
    stats = bam_dir + "{sample}.sorted.mrkdup.txt",
  shell:
    "picard MarkDuplicates "
        "INPUT={input} "
        "OUTPUT={output.mrkdup} "
        "VALIDATION_STRINGENCY=SILENT "
        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 "
        "REMOVE_DUPLICATES=FALSE "
        "METRICS_FILE='{log.stats}'; "
    "samtools index {output.mrkdup}; "


######## GATK REALIGNER TARGET CREATOR #########

rule GATK_RealignerTargetCreator:
  input:
    bam = bam_dir + "{sample}" + ".sorted.mrkdup.bam",
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
    #fa = config["path"]["genomefa"] + config["references"]["genomefa"],
    dbindels = config["path"]["genomefa"] + config["references"]["dbsnp"],
  output:
    bam_dir + "{sample}" + ".sorted.mrkdup.bam" + ".intervals",
  shell:
    "gatk "
        "-T RealignerTargetCreator "
        "-I {input.bam} "
        "-R {input.fa} "
        "--known {input.dbindels} "
        "--out {output}; "

######## GATK INDEL REALIGNMENT #########

rule GATK_IndelRealigner:
  input:
    bam = bam_dir + "{sample}" + ".sorted.mrkdup.bam",
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
    knownsites = config["path"]["genomefa"] + config["references"]["dbsnp"],
    intervals = bam_dir + "{sample}" + ".sorted.mrkdup.bam" + ".intervals",
  output:
    bam_dir + "{sample}" + ".sorted.mrkdup.ralgn.bam"
  shell:
    "gatk "
        "-T IndelRealigner "
        "-I {input.bam} "
        "-R {input.fa} "
        "-known {input.knownsites} "
        "-targetIntervals {input.intervals} "
        "--out {output}"

######### GATK BASE RECALIBRATOR ##########

rule GATK_baserecalibrator:
  input:
    bam = bam_dir + "{sample}" + ".sorted.mrkdup.ralgn.bam",
    #fa = config["path"]["genomefa"] + config["references"]["genomefa"],
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
    knownsites = config["path"]["genomefa"] + config["references"]["dbsnp"],
  output:
    bam_bsrcl = bam_dir + "{sample}" + ".sorted.mrkdup.ralgn.bsrcl.bam",
    bam_bsrcl_list = bam_dir + "{sample}" + ".sorted.mrkdup.ralgn.bsrcl.bam.list",
  shell:
    "gatk "
        "-T BaseRecalibrator "
        "-R {input.fa} "
        "--knownSites {input.knownsites} "
        "-I {input.bam} "
        "--out {output.bam_bsrcl_list}; "
    "gatk "
        "-T PrintReads "
        "-R {input.fa} "
        "-I {input.bam} "
        "--out {output.bam_bsrcl}; "
    "samtools index {input.bam};"

rule cp_tumor_bam:
  input:
    #bam_bsrcl = bam_dir + "{sample}" + ".sorted.mrkdup.ralgn.bsrcl.bam",
    #bam_bsrcl = expand(bam_dir + "{sample}" + ".sorted.mrkdup.ralgn.bsrcl.bam", sample = config["samples"])
    bam_mrkdup = expand(bam_dir + "{sample}" + ".sorted.mrkdup.bam", sample = config["samples"]),
  params:
    bam2rm = bam_dir,
  output:
    tumorBam = bam_dir + "tumor.bam",
  shell:
    "cp {input.bam_mrkdup} {output.tumorBam}; "
    "samtools index {output.tumorBam};"
#rm {params.bam2rm}*.sorted.bam; rm {params.bam2rm}*.sorted.mrkdup*;"


########## PICARD COLLECTHSMETRICS ###########

rule picard_collectHSmetric:
  input:
    fadict = (config["path"]["genomefa"] + config["references"]["genomefa_GATK"]).replace(".fasta",".dict"),
    bed = config["path"]["panel"] + config["bed"]["variant_panel"],
    bam = bam_dir + "tumor.bam",
  output:
    hsmetric = result_dir + "{sample}.sorted.mrkdup.hsmetric",
    bed2interval = result_dir + "{sample}" + ".sorted.mrkdup.bam.picard.bedintervals",
  shell:
    "picard "
      "BedToIntervalList "
      "I={input.bed} "
      "O={output.bed2interval} "
      "SD={input.fadict}; "
    "picard "
      "CollectHsMetrics "
      "BI={output.bed2interval} "
      "TI={output.bed2interval} "
      "I={input.bam} "
      "O={output.hsmetric} "
      "COVERAGE_CAP=500 "
      "METRIC_ACCUMULATION_LEVEL=ALL_READS "
      "METRIC_ACCUMULATION_LEVEL=LIBRARY"


######### PICARD COLLECTALIGNMENTSUMMARYMETRICS ##########

rule picard_CollectAlignmentSummaryMetrics:
  input:
    bam = bam_dir + "{sample}.sorted.bam",
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
  output:
    metric = result_dir + "{sample}.sorted.alignmetric",
  params:
    adapter = config["QC"]["adapter"],
  shell:
    "picard "
      "CollectAlignmentSummaryMetrics "
      "R={input.fa} "
      "I={input.bam} "
      "O={output.metric} "
      "ADAPTER_SEQUENCE={params.adapter} "
      "METRIC_ACCUMULATION_LEVEL=ALL_READS "
      "METRIC_ACCUMULATION_LEVEL=LIBRARY;"



##############################################################
#####################  VARIANT CALLING #######################
##############################################################

#####################  Mutect 2  #############################

rule mutect2_somatic_tumor_only:
  input:
    fa = config["path"]["genomefa"] + config["references"]["genomefa_GATK"],
    dbsnpALL = config["path"]["genomefa"] + config["references"]["dbsnpALL"],
    COSMICcodMut = config["path"]["genomefa"] + config["references"]["cosmic"],
    bamTumor =  bam_dir + "tumor.bam",
    bed = config["path"]["panel"] + config["bed"]["variant_panel"],
    germline = config["path"]["genomefa"] + config["references"]["germline"]
  output:
    m2_vcf = vcf_dir + config["vcf"]["mutect"]["default"],
  params:
    tumor_sample_name = list(config["samples"].keys())[0],
    gatk4 = config["analysis"]["gatk4"]
  threads: 6
  shell:
    "java -jar {params.gatk4} Mutect2 "
        "-R {input.fa} "
        "-I {input.bamTumor} "
        "-tumor {params.tumor_sample_name} "
        "-L {input.bed} "
        "-O {output.m2_vcf};"
#        "--germline-resource {input.germline} "  # Add this only if totally sure that some variants of interest will not be filtered out

rule splitmultiallelic:
  input:
    m2_vcf = vcf_dir + config["vcf"]["mutect"]["default"]
  output:
    m2_vcf_gz = vcf_dir + config["vcf"]["mutect"]["merged"],  # Do not take much into account the name "merged". Nothing was merged, but only a VCF compressed
    m2_vcf_splitalleles = vcf_dir + "split_" + config["vcf"]["mutect"]["default"]
  shell:
    "bgzip -c {input.m2_vcf} > {output.m2_vcf_gz}; tabix -p vcf {output.m2_vcf_gz};" # Be aware of the use of bcftools norm to split multiallelic locations
    " bcftools norm -m -both -o {output.m2_vcf_splitalleles} -Ov {output.m2_vcf_gz}" # splitting multiallelic sites in the VCF



rule mutect2_filterMutectCalls:
  input:
   m2_vcf_splitalleles = vcf_dir + "split_" + config["vcf"]["mutect"]["default"]
  output:
    filt_m2_vcf = vcf_dir + "filt_split_" + config["vcf"]["mutect"]["default"]
  params:
    gatk4 = config["analysis"]["gatk4"]
  threads: 6
  shell:
    "java -jar {params.gatk4} FilterMutectCalls "
        "-V {input.m2_vcf_splitalleles} "
        "-O {output.filt_m2_vcf};"

#############################################################
#####################  VCF ANNOTATION #######################
#############################################################

##############  VARIANT-EFFECT-PREDICTOR  ###################


rule vep_annotation_Mutect2:
  input:
    filt_m2_vcf = vcf_dir + "filt_split_" + config["vcf"]["mutect"]["default"],
  output:
    vep_vcf = vep_dir + "vep_" + config["vcf"]["mutect"]["default"],
  params:
    vep_exec = config["path"]["vep"],
    processVEP = config["path"]["code"] + "VEPoutputparser.py",
    searchmutations = config["path"]["code"] + "searchmutations.py",
    outsearchmut = '/home/jlanillos/Disco4tb/Projects/miSeq/results/VCF/',
    outVEPprocess = vep_dir # Where is 'FILTERED_annotations.csv' (created in VEPoutputparser.py
    #vep_refseq = config["references"]["vep_cache_refseq"]
  threads: 4
  shell:
    "{params.vep_exec}vep "
        "-i {input.filt_m2_vcf} "
        "-o {output.vep_vcf} "
        "--fork {threads} "
        "--vcf "
        "--poly p "
        "--sift p "
        "--variant_class "
        "-format vcf "
        "--offline "
        "-v "
        "--force_overwrite "
        "--use_given_ref "
        "--af_gnomad; "
        "python {params.processVEP} "
        "--file {output.vep_vcf}; "
        "python {params.searchmutations} --file {params.outVEPprocess} --out {params.outsearchmut}"
        #       "--refseq {params.vep_refseq}   "


  #!python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#__author__ = "Javier Lanillos"

# This file contains the whole set of rules to run the RNAseq pipeline necessary for neoantigen prediction using pvacseq.
# The goal is to start from input tumor paired end RNAseq fastq files and obtain three different types of outputs:
    # VCF files
    # Cufflinks output files: genes.fpkm, isoforms.fpmk
    # Bam-readcount: indels.bam-read-count and snvs.bam-readcount
# Also, this pipeline will keep intermediate files such as BAM (sorted, filtered...)
# These output files will be used as input files for pvacseq


rule_dir = "../"


include:
  rule_dir + "snakemake_rules/quality_control/cutadapt.rule"
include:
  rule_dir + "snakemake_rules/align/hisat2.rule"
include:
  rule_dir + "snakemake_rules/align/samtools.rule"
include:
  rule_dir + "snakemake_rules/quality_control/fastqc.rule"
include:
  rule_dir + "snakemake_rules/align/samtools.rule"
include:
  rule_dir + "snakemake_rules/quality_control/picard.rule"
include:
  rule_dir + "snakemake_rules/quality_control/sambamba_depth.rule"
include:
  rule_dir + "snakemake_rules/quality_control/collectqc.rule"
include:
  rule_dir + "snakemake_rules/quality_control/GATK.rule"
include:
  rule_dir + "snakemake_rules/transcript_analysis/cufflinks.rule
include:
  rule_dir + "snakemake_rules/transcript_analysis/bamreadcount.rule"


cutadapt_dir = config["analysis"]["analysis_dir"] + config["analysis"]["sample_id"]  + "/" + config["analysis"]["result"]  + "cutadapt/"
fastqc_dir = config["analysis"]["analysis_dir"] + config["analysis"]["sample_id"]  + "/" + config["analysis"]["result"]  + "fastqc/"
bam_dir = config["analysis"]["analysis_dir"]  + config["analysis"]["sample_id"]  + "/" + config["analysis"]["result"]  + "bam/"
vcf_dir = config["analysis"]["analysis_dir"]  + config["analysis"]["sample_id"]  + "/" + config["analysis"]["result"]  + "vcf/"
bamreadcount_out_dir = config["analysis"]["analysis_dir"] + config["analysis"]["sample_id"] + "/" + config["analysis"]["result"] + "bamreadcount/"
cufflinks_out = config["analysis"]["analysis_dir"] + config["analysis"]["sample_id"] + "/" + config["analysis"]["result"] + "cuffout/"
result_dir = config["analysis"]["analysis_dir"]  + config["analysis"]["sample_id"]  + "/" + config["analysis"]["result"]

manta_dir = "manta/"
strelka_dir = "strelka/"
vardict_dir = "vardict/"

rule all:
  input:
    #expand(bam_dir + "{mysample}.bam", mysample=config["samples"]),
    #expand(bam_dir + "{mysample}.sorted.bam", mysample=config["samples"]),
    #expand(bam_dir + "{mysample}.sorted.rmdup.bam", mysample=config["samples"]),
    #expand(bam_dir + "{mysample}.sorted.rmdup.hsmetric", mysample=config["samples"]),
    #expand(bam_dir + "{mysample}.sorted.rmdup.cov.bed", mysample=config["samples"]),
    #expand(bam_dir + "{mysample}.sorted.rmdup.intervals", mysample=config["samples"]),
    #expand(bam_dir + "{mysample}.sorted.rmdup.ralgn.bam", mysample=config["samples"]),
    #expand(bam_dir + "{mysample}.sorted.rmdup.ralgn.bsrcl.bam", mysample=config["samples"]),
    result_dir + "qc/" + "qc_report.pdf",
    vcf_dir + "combined_vcf.vcf",
    #expand(vcf_dir + manta_dir + "{vcf_file}", vcf_file=config["vcf"]["manta"]),
    #expand(vcf_dir + strelka_dir + "{vcf_file}", vcf_file=config["vcf"]["strelka"]),
    #vcf_dir + vardict_dir + config["vcf"]["vardict"],

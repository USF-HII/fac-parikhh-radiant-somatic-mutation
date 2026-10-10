from os.path import abspath, basename, dirname, join, splitext

shell.prefix("set -euo pipefail; set -x;")

TOP = os.getcwd() # BASE
BIO = os.environ["BIO"]
RUN = os.environ.get("RUN", "main")
WRK = join(TOP, "tmp", "wrk", RUN)
MTR = join(TOP, "tmp", "mtr", RUN) # metrics

def top(*path):
    return os.path.join(TOP, *path)

def bio(*path):
    return os.path.join(BIO, *path)

def mtr(*path):
    return os.path.join(MTR, *path)

def wrk(*path):
    return os.path.join(WRK, *path)

#----------------------------------------------------------------------------------------------------
# helper funcs
#----------------------------------------------------------------------------------------------------

def build_db():
    db = {}
    with open(top("mta", "manifest.tsv")) as f:
        for line in f:
            sample, sex, path = line.rstrip("\n").split("\t")
            db[sample] = { "sex": sex, "path": bio("lab", path) }
    return db

#----------------------------------------------------------------------------------------------------
# variables
#----------------------------------------------------------------------------------------------------

DB = build_db()

MOSAICHUNTER_BLAT = bio("ref", "prj", "somatic-mutation", "mosaichunter", "blat-v369")
MOSAICHUNTER_JAR = bio("ref", "prj", "somatic-mutation", "mosaichunter", "mosaichunter-4bdadaa7.jar")

HG38_REF = bio("ref", "broad", "hg38", "v0", "Homo_sapiens_assembly38.fasta")
HG38_DBSNP_VCF = bio("ref", "broad", "hg38", "v0", "Homo_sapiens_assembly38.dbsnp138.vcf.gz")

HG38_REP_REG_BED = top("tmp", "data", "ucsc_repetitive_region_hg38.bed")
HG38_INDEL_BED = top("tmp", "data", "ucsc_indel_hg38.bed")
HG38_COMMON_BED = top("tmp", "data", "ucsc_common_hg38.bed")

UCSC_HG38_SUPDUP = bio("ref", "ucsc", "hg38", "database", "genomicSuperDups.txt.gz")
UCSC_HG38_RMSK = bio("ref", "ucsc", "hg38", "database", "rmsk.txt.gz")
UCSC_HG38_COMM = bio("ref", "ucsc", "hg38", "database", "snp151Common.txt.gz")

if os.environ.get("CHROMOSOMES"):
    CHROMOSOMES = os.environ["CHROMOSOMES"].split(",")
else:
    CHROMOSOMES = [f"chr{n}" for n in range(1, 23)] + ["chrX"]

if "." in RUN:
    INP_PCT = RUN.split(".")[1]
else:
    INP_PCT = ""

#----------------------------------------------------------------------------------------------------
# containers
#----------------------------------------------------------------------------------------------------

CONTAINER_CMD = "bio-img"

JAVA = " ".join([CONTAINER_CMD, bio("img", "prd", "gen", "openjdk", "8.simg"), "java"])
SAMTOOLS = " ".join([CONTAINER_CMD, bio("img", "prd", "gen", "htslib", "1.23.simg"), "samtools"])

#----------------------------------------------------------------------------------------------------
# targets
#----------------------------------------------------------------------------------------------------

rule target_main:
    input:
        expand(wrk("out", "mosaichunter", "{sample}.tsv"), sample=DB),

#-----------------------------------------------------------------------------------------------------
# pragmas
#-----------------------------------------------------------------------------------------------------

wildcard_constraints: chr="[^/]+"
wildcard_constraints: sample="[^/]+"

localrules: mosaichunter_agg

#----------------------------------------------------------------------------------------------------
# rules
#----------------------------------------------------------------------------------------------------

rule inp_bam:
    input:
        ref=HG38_REF,
        cram=lambda wc: DB[wc.sample]["path"],
    output:
        bam=wrk("inp", "bam", "{sample}.bam"),
    resources:
        mem="5G",
        runtime="8h",
    params:
        pct=INP_PCT,
    threads:
        1
    benchmark:
        mtr("inp", "bam", "{sample}.txt"),
    shell:
        """
        if [[ -n "{params.pct}" ]]; then
          _opts="-s 42.{params.pct}"
        else
          _opts=""
        fi

        {SAMTOOLS} view $_opts --no-PG --reference={input.ref} --bam --output={output.bam} {input.cram}

        {SAMTOOLS} index --bai {output.bam}

        /bin/touch --date='+1 hour' {output.bam}.bai
        """

rule inp_ref:
    input:
        ref=HG38_REF,
        split_ref_py=top("scripts", "split-ref.py"),
    output:
        ref=wrk("inp", "ref", "{chr}.fa")
    resources:
       mem="5G",
       runtime="1h",
    threads:
        1
    benchmark:
        mtr("inp", "ref", "{chr}.txt"),
    shell:
        """
        python3 {input.split_ref_py} {input.ref} {wildcards.chr} > {output.ref}

        {SAMTOOLS} faidx {output.ref}
        """

rule inp_ucsc:
    input:
        supdup=UCSC_HG38_SUPDUP,
        rmsk=UCSC_HG38_RMSK,
        comm=UCSC_HG38_COMM,
    output:
        comm=wrk("inp", "ucsc", "common.bed"),
        indel=wrk("inp", "ucsc", "indel.bed"),
        repreg=wrk("inp", "ucsc", "repreg.bed"),
        segdup=wrk("inp", "ucsc", "sedgup.bed"),
    shell:
        """
        gunzip -c {input.comm}   | cut -f2-4 | sort --version-sort > {output.comm}

        gunzip -c {input.comm} \
          | awk '$12 == "in-del" || $12 == "insertion" || $12 == "deletion" {{ print }}' \
          | cut -f2-4 \
          | sort --version-sort \
          > {output.indel}

        gunzip -c {input.supdup} | cut -f2-4 | sort --version-sort > {output.segdup}

        gunzip -c {input.rmsk}   | cut -f6-8 | sort --version-sort > {output.repreg}
        """

rule mosaichunter:
    input:
        jar=MOSAICHUNTER_JAR,
        blat=MOSAICHUNTER_BLAT,
        ref=wrk("inp", "ref", "{chr}.fa"),
        dbsnp=HG38_DBSNP_VCF,
        rep_reg_bed=wrk("inp", "ucsc", "repreg.bed"),
        indel_bed=wrk("inp", "ucsc", "indel.bed"),
        common_bed=wrk("inp", "ucsc", "common.bed"),
        bam=wrk("inp", "bam", "{sample}.bam"),
    output:
        tsv=wrk("out", "mosaichunter", "{sample}", "{chr}.tsv"),
    resources:
       mem="20G",
       runtime="21d",
    threads:
        1
    params:
        sex=lambda wc: DB[wc.sample]["sex"],
        java_opts="-Xmx20g",
    shadow:
        "shallow"
    benchmark:
        mtr("out", "mosaichunter", "{sample}", "{chr}.txt"),
    shell:
        """
        {JAVA} {params.java_opts} -jar {input.jar} \
          genome \
          -P misaligned_reads_filter.blat_path={input.blat} \
          -P input_file={input.bam} \
          -P reference_file={input.ref} \
          -P mosaic_filter.dbsnp_file={input.dbsnp} \
          -P repetitive_region_filter.bed_file={input.rep_reg_bed} \
          -P indel_region_filter.bed_file={input.indel_bed} \
          -P common_site_filter.bed_file={input.common_bed} \
          -P mosaic_filter.sex={params.sex} \
          -P valid_references={wildcards.chr} \
          -P chr={wildcards.chr} \
          -P output_dir=$(pwd)

        /bin/mv final.passed.tsv {output.tsv}
        """

rule mosaichunter_agg:
    input:
        tsvs=expand(wrk("out", "mosaichunter", "{{sample}}", "{chr}.tsv"), chr=CHROMOSOMES),
    output:
        tsv=wrk("out", "mosaichunter", "{sample}.tsv")
    shell:
        """
        header="
        chr
        pos
        refbase
        depth
        bases
        baseqs
        majorbase
        majordepth
        minorbase
        minordepth
        major_minor_allele_with_dbsnp_freq
        log10_prior_major_homozygous
        log10_prior_heterozygous
        log10_prior_minor_homozygous
        log10_prior_mosaic
        log10_likelihood_major_homozygous
        log10_likelihood_heterozygous
        log10_likelihood_minor_homozygous
        log10_likelihood_mosaic
        log10_posterior_major_homozygous
        log10_posterior_heterozygous
        log10_posterior_minor_homozygous
        log10_posterior_mosaic
        mosaic_posterior_probability
        "

        echo $header | tr ' ' '\t' > {output.tsv}

        cat {input.tsvs} >> {output.tsv}
        """

# vim: ft=snakemake tabstop=4 shiftwidth=4 softtabstop=0 expandtab

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
# variables
#----------------------------------------------------------------------------------------------------

PFX = bio("lab", "radiant", "broad")

DB = {
  "SM-L5BY4": join(PFX, "Human-WGS-BAM-1/RADIANT_Set001-002_41samples/RADIANT_Set001-002_41samples_BAM/134542-0359493956_PM21-00438-A_SM-L5BY4_v1_WGS_GCP.bam"),
}

MOSAICHUNTER_BLAT = bio("ref", "prj", "somatic-mutation", "mosaichunter", "blat-v369")
MOSAICHUNTER_JAR = bio("ref", "prj", "somatic-mutation", "mosaichunter", "mosaichunter-4bdadaa7.jar")

HG38_REF = bio("ref", "broad", "hg38", "v0", "Homo_sapiens_assembly38.fasta")
HG38_DBSNP_VCF = bio("ref", "broad", "hg38", "v0", "Homo_sapiens_assembly38.dbsnp138.vcf.gz")

HG38_REP_REG_BED = top("tmp", "data", "ucsc_repetitive_region_hg38.bed")
HG38_INDEL_BED = top("tmp", "data", "ucsc_indel_hg38.bed")
HG38_COMMON_BED = top("tmp", "data", "ucsc_common_hg38.bed")

CHROMOSOMES = ["chr21", "chr22"]

#----------------------------------------------------------------------------------------------------
# containers
#----------------------------------------------------------------------------------------------------

CONTAINER_CMD = "bio-img"

JAVA8 = " ".join(CONTAINER_CMD, bio("img", "prd", "gen", "openjdk", "8.simg"), "java8")

#----------------------------------------------------------------------------------------------------
# targets
#----------------------------------------------------------------------------------------------------

rule target_main:
    input:
        expand(wrk("mosaichunter", "{sample}", "{chr}.tsv"), sample=DB, chr=CHROMOSOMES),

#-----------------------------------------------------------------------------------------------------
# pragmas
#-----------------------------------------------------------------------------------------------------

wildcard_constraints: chr="[^/]+"
wildcard_constraints: sample="[^/]+"

#----------------------------------------------------------------------------------------------------
# rules
#----------------------------------------------------------------------------------------------------

rule mosaichunter:
    input:
        jar=MOSAICHUNTER_JAR,
        blat=MOSAICHUNTER_BLAT,
        ref=HG38_REF,
        dbsnp=HG38_DBSNP_VCF,
        rep_reg_bed=HG38_REP_REG_BED,
        indel_bed=HG38_INDEL_BED,
        common_bed=HG38_COMMON_BED,
        bam=lambda wc: DB[wc.sample],
    output:
        tsv=wrk("mosaichunter", "{sample}", "{chr}.tsv"),
    resources:
       mem="96G",
       runtime="7d",
    threads:
        1
    params:
        work=wrk("mosaichunter", "{sample}", "{chr}.d"),
        sex="M",
    benchmark:
        mtr("mosaichunter", "{sample}", "{chr}.txt"),
    shell:
        """
        /bin/mkdir -p {params.work} && cd {params.work}

        {JAVA8} -Xmx90G -jar {input.jar} \
          genome \
          -P misaligned_reads_filter.blat_path={input.blat} \
          -P input_file={input.bam}
          -P reference_file={input.ref} \
          -P mosaic_filter.dbsnp_file={input.dbsnp} \
          -P repetitive_region_filter.bed_file={input.rep_reg_bed} \
          -P indel_region_filter.bed_file={input.indel_bed} \
          -P common_site_filter.bed_file={input.common_bed} \
          -P mosaic_filter.sex={params.sex} \
          -P valid_references={wildcards.chr} \
          -P output_dir=$(pwd)

        /bin/mv final.passed.tsv {output.tsv}
        """

# vim: ft=snakemake tabstop=4 shiftwidth=4 softtabstop=0 expandtab

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
    "SM-L5BY4": { "id": "134542", "sex": "F", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set001-002_41samples/134542-0359493956_PM21-00438-A_SM-L5BY4_v1_WGS_GCP.bam"), },
    "SM-N8UGS": { "id": "134198", "sex": "M", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/134198-0359493977_PM24-00181-A_SM-N8UGS.bam"), },
    "SM-N8UGU": { "id": "135718", "sex": "M", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135718-0052920102_PM24-00306-A_SM-N8UGU.bam"), },
    "SM-N8UGV": { "id": "135634", "sex": "F", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135634-0052670702_PM24-00307-A_SM-N8UGV.bam"), },
    "SM-N8UGW": { "id": "135685", "sex": "M", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135685-0052729302_PM24-00308-A_SM-N8UGW.bam"), },
    "SM-N8URS": { "id": "135681", "sex": "F", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135681-0052614402_PM24-00312-A_SM-N8URS.bam"), },
    "SM-N8URT": { "id": "135347", "sex": "M", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135347-0052714802_PM24-00311-A_SM-N8URT.bam"), },
    "SM-N8URU": { "id": "135753", "sex": "M", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135753-0052786102_PM24-00310-A_SM-N8URU.bam"), },
    "SM-N8URV": { "id": "135700", "sex": "M", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135700-0052614602_PM24-00309-A_SM-N8URV.bam"), },
    "SM-NL7S4": { "id": "135635", "sex": "M", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135635-0052902102_PM24-00313-A_SM-NL7S4.bam"), },
    "SM-NL7S5": { "id": "135577", "sex": "F", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135577-0052786502_PM24-00315-A_SM-NL7S5.bam"), },
    "SM-NL7S6": { "id": "135264", "sex": "F", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135264-0052831502_PM24-00316-A_SM-NL7S6.bam"), },
    "SM-NL7S7": { "id": "135508", "sex": "F", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135508-0052786302_PM24-00317-A_SM-NL7S7.bam"), },
    "SM-NL7S8": { "id": "135169", "sex": "F", "path": join(PFX, "Human-WGS-BAM-1/RADIANT_Set035_14samples/RADIANT_Set035_14samples_BAM/135169-0052902302_PM24-00318-A_SM-NL7S8.bam"), },
}

MOSAICHUNTER_BLAT = bio("ref", "prj", "somatic-mutation", "mosaichunter", "blat-v369")
MOSAICHUNTER_JAR = bio("ref", "prj", "somatic-mutation", "mosaichunter", "mosaichunter-4bdadaa7.jar")

HG38_REF = bio("ref", "broad", "hg38", "v0", "Homo_sapiens_assembly38.fasta")
HG38_DBSNP_VCF = bio("ref", "broad", "hg38", "v0", "Homo_sapiens_assembly38.dbsnp138.vcf.gz")

HG38_REP_REG_BED = top("tmp", "data", "ucsc_repetitive_region_hg38.bed")
HG38_INDEL_BED = top("tmp", "data", "ucsc_indel_hg38.bed")
HG38_COMMON_BED = top("tmp", "data", "ucsc_common_hg38.bed")

CHROMOSOMES = ["chr1", "chr2", "chr21", "chr22", "chrX"]

#----------------------------------------------------------------------------------------------------
# containers
#----------------------------------------------------------------------------------------------------

CONTAINER_CMD = "bio-img"

JAVA = " ".join([CONTAINER_CMD, bio("img", "prd", "gen", "openjdk", "8.simg"), "java"])

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

        {JAVA} -Xmx90G -jar {input.jar} \
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
          -P output_dir=$(pwd)

        /bin/mv final.passed.tsv {output.tsv}
        """

# vim: ft=snakemake tabstop=4 shiftwidth=4 softtabstop=0 expandtab

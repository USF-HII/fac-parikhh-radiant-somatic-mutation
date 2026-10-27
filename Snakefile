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

MOSAICHUNTER_HEADER = """
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
""".strip().split()

MOSAICHUNTER_INPUT_INDICES = [MOSAICHUNTER_HEADER.index(f) for f in ["chr", "pos", "pos", "majorbase", "minorbase"]]

ANNOVAR_SCRIPT = bio("ref", "annovar.2503", "table_annovar.pl")

ANNOVAR_DB = bio("ref", "annovar", "humandb")

ANNOVAR_SPEC = """
    f:avsnp151
    g:refGene
    f:clinvar_20250721
    f:dbnsfp47a
    f:dbnsfp47a_interpro
    f:gnomad41_exome
    f:gnomad41_genome
    f:GTEx_v8_eQTL
    f:GTEx_v8_sQTL
    f:intervar_20250721
    f:logofunc
""".strip().split()

ANNOVAR_OPERATION = ",".join([i.split(":")[0] for i in ANNOVAR_SPEC])
ANNOVAR_PROTOCOL  = ",".join([i.split(":")[1] for i in ANNOVAR_SPEC])

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

ANNOVAR = " ".join(["/usr/bin/env", "LC_ALL=C", CONTAINER_CMD, bio("img", "prd", "gen", "perl", "5.40.1.simg"), ANNOVAR_SCRIPT])
JAVA = " ".join([CONTAINER_CMD, bio("img", "prd", "gen", "openjdk", "8.simg"), "java"])
SAMTOOLS = " ".join([CONTAINER_CMD, bio("img", "prd", "gen", "htslib", "1.23.simg"), "samtools"])

#----------------------------------------------------------------------------------------------------
# targets
#----------------------------------------------------------------------------------------------------

if os.environ.get("SAMPLES"):
    SAMPLES = os.environ["SAMPLES"].split(",")
    DB = { k: v for k, v in DB.items() if k in SAMPLES }

rule target_main:
    input:
        expand(wrk("out", "mosaichunter_annotate_final_pivot", "{sample}.txt"), sample=DB),
        expand(wrk("out", "mosaichunter_annotate_final", "{sample}.tsv"), sample=DB),
        expand(wrk("out", "mosaichunter", "{sample}.tsv"), sample=DB),

#-----------------------------------------------------------------------------------------------------
# pragmas
#-----------------------------------------------------------------------------------------------------

wildcard_constraints: chr="[^/]+"
wildcard_constraints: sample="[^/]+"

localrules: mosaichunter_agg
localrules: mosaichunter_annotate_input
localrules: mosaichunter_annotate_final
localrules: mosaichunter_annotate_final_pivot

#----------------------------------------------------------------------------------------------------
# rules
#----------------------------------------------------------------------------------------------------

rule inp_bam:
    input:
        ref=HG38_REF,
        cram=lambda wc: DB[wc.sample]["path"],
    output:
        bam=temp(wrk("inp", "bam", "{sample}.bam")),
    resources:
        mem="5G",
        runtime="8h",
        io=1,
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
       attempt=lambda wc, attempt: attempt, # params then uses resources.attempt
       mem=lambda wc, attempt: f"{20 * attempt}G",
       runtime=lambda wc, attempt: f"{attempt}d",
    threads:
        1
    params:
        sex=lambda wc: DB[wc.sample]["sex"],
        java_opts=lambda wc, resources: f"-Xmx{20 * resources.attempt}g",
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
        cat {input.tsvs} >> {output.tsv}
        """

rule mosaichunter_annotate_input:
    input:
        tsv=wrk("out", "mosaichunter", "{sample}.tsv"),
    output:
        tsv=wrk("out", "mosaichunter_annotate_input", "{sample}.tsv"),
    resources:
        mem="10G",
        runtime="4h",
    threads:
        1
    benchmark:
        mtr("out", "mosaichunter_annotate_input", "{sample}.txt"),
    run:
        with open(input.tsv) as f, open(output.tsv, "w") as f_out:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                fields = [fields[i] for i in MOSAICHUNTER_INPUT_INDICES] + fields # add fields for annovar lookup at beginning

                print("\t".join(fields), file=f_out)

rule mosaichunter_annotate:
    input:
        annovar_db=ANNOVAR_DB,
        tsv=wrk("out", "mosaichunter_annotate_input", "{sample}.tsv"),
    params:
        annovar_protocol=ANNOVAR_PROTOCOL,
        annovar_operation=ANNOVAR_OPERATION,
    output:
        tsv=wrk("out", "mosaichunter_annotate", "{sample}.tsv"),
    shadow:
        "shallow"
    resources:
        mem="10G",
        runtime="4h",
    threads:
        1
    benchmark:
        mtr("out", "moscaichunter_annotate", "{sample}.txt"),
    shell:
        """
        {ANNOVAR} \
          --outfile=out \
          --otherinfo \
          --nastring=NA \
          --polish \
          --buildver=hg38 \
          --protocol={params.annovar_protocol} \
          --operation={params.annovar_operation} \
          {input.tsv} \
          {input.annovar_db}

        /bin/cp out.hg38_multianno.txt {output.tsv}
        """

rule mosaichunter_annotate_final:
    input:
        tsv=wrk("out", "mosaichunter_annotate", "{sample}.tsv"),
    output:
        tsv=wrk("out", "mosaichunter_annotate_final", "{sample}.tsv"),
    resources:
        mem="10G",
        runtime="4h",
    threads:
        1
    benchmark:
        mtr("out", "mosaichunter_annotate_final", "{sample}.txt"),
    run:
        with open(input.tsv) as f, open(output.tsv, "w") as f_out:
            header = f.readline().rstrip("\n").split("\t")
            header = MOSAICHUNTER_HEADER + header[len(MOSAICHUNTER_INPUT_INDICES):-len(MOSAICHUNTER_HEADER)]
            print("\t".join(header), file=f_out)

            for line in f:
                fields = line.rstrip("\n").split("\t")
                fields = fields[-len(MOSAICHUNTER_HEADER):] + fields[len(MOSAICHUNTER_INPUT_INDICES):-len(MOSAICHUNTER_HEADER)]
                print("\t".join(fields), file=f_out)

rule mosaichunter_annotate_final_pivot:
    input:
        tsv=wrk("out", "mosaichunter_annotate_final", "{sample}.tsv"),
    output:
        txt=wrk("out", "mosaichunter_annotate_final_pivot", "{sample}.txt")
    resources:
        mem="10G",
        runtime="4h",
    threads:
        1
    benchmark:
        mtr("out", "mosaichunter_annotate_final_pivot", "{sample}.txt"),
    run:
        with open(input.tsv) as f, open(output.txt, "w") as f_out:
            header = f.readline().rstrip("\n").split("\t")
            length = max(len(f) for f in header)

            for line in f:
                fields = line.rstrip("\n").split("\t")

                for col, val in zip(header, fields):
                    print(f"{col:<{length}}\t{val}", file=f_out)

                print(file=f_out)

# vim: ft=snakemake tabstop=4 shiftwidth=4 softtabstop=0 expandtab

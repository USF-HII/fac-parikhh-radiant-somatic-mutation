#!/usr/bin/env bash

# ----------------------------------------------------------------------------------------------------
# util funcs
# ----------------------------------------------------------------------------------------------------

_tmp() { # tmp() [path] - makes dir under ${TMP}/<task>[/path] and echoes path for capture
  local task=$(echo ${FUNCNAME[1]} | sed 's/^task_//')
  local dir="${TMP}/${task}"
  if [[ $# -eq 1 ]]; then dir=${dir}/$1; fi
  if ! [[ -d ${dir} ]]; then /bin/mkdir -p ${dir}; fi
  echo ${dir}
}

_pivot() {
  local tsv=$1
  local n=${2:-2}

  paste \
    <(sed -n    1p ${tsv} | tr '\t' '\n') \
    <(sed -n ${n}p ${tsv} | tr '\t' '\n') \
  |  column --separator=$'\t' --table
}

# ----------------------------------------------------------------------------------------------------
# rad
# ----------------------------------------------------------------------------------------------------

rad_init() {
  echo tsv=/home/bioinfo/outbox/parikhh/for_kevin/radiant_wgs/radiant_participant_sex_2026_09_28.txt
  echo bam_idxs=\"$(ls ${BIO}/lab/radiant/broad/Human-WGS-BAM-[0-9]*.index)\"
  echo cram_idxs=\"$(ls ${BIO}/lab/radiant/broad/Human-WGS-CRAM-[0-9]*.index)\"
}

rad_tsv_ids() {
  source <(rad_init)
  tail -n+2 ${tsv} | cut -f1
}

rad_tsv_ids_sex() {
  # note we remove the few entries with NA for sex
  source <(rad_init)
  tail -n+2 ${tsv} | grep -v '[[:space:]]NA[[:space:]]' \
    | awk '{ print $1 "\t" ($3==1 ? "M" : "F")}' \
    | sort --version-sort
}

rad_cram_ls() {
  source <(rad_init)

  for idx in ${cram_idxs}; do
    pfx=$(echo radiant/broad/$(basename $idx .index))
    cut -d: -f1 $idx | grep '\.cram$' | sed "s|^|${pfx}/|"
  done | sort --version-sort
}

rad_cram_man() {
  tmp=$(_tmp) && /bin/mkdir -p ${tmp} && cd ${tmp}
  rad_cram_ls > cram.txt

  while read id sex; do
    if [[ $id  == "SM-L5KUV" ]]; then # corrupt
      continue
    fi

    path=$(grep ${id} cram.txt)

    if [[ -n ${path} ]]; then
      echo $id $sex $path | tr ' ' '\t'
    fi
  done < <(rad_tsv_ids_sex)
}

# ----------------------------------------------------------------------------------------------------
# global vars
# ----------------------------------------------------------------------------------------------------

TOP=$(readlink -f $(dirname "$0"))
TMP=${TOP}/tmp/task

# ----------------------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------------------

if [[ $# -lt 1 ]]; then
  sed -n 's/^\([a-zA-Z_][a-zA-Z0-9_]*\) *().*/\1/p' "$0" | grep -v '^_' | sort --version-sort
  exit 1
else
  "$@"
fi

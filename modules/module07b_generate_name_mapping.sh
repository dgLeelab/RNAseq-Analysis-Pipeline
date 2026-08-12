#!/usr/bin/env bash
# =========================================================
# module07b_generate_name_mapping.sh
# LOC -> gene name 재매핑표 생성 (해석/annotation 용, 분석 자체에는 불필요)
#
# 사용: ./module07b_generate_name_mapping.sh <nuclear.gff> <출력 tsv>
# =========================================================
set -euo pipefail

GFF="${1:?사용법: $0 <nuclear.gff> <output loc2name.tsv>}"
OUTFILE="${2:-loc2name.tsv}"

awk -F'\t' '$3=="gene"{
  gname=""; gnum="";
  if(match($9,/gene=[^;]+/))    gname=substr($9,RSTART+5,RLENGTH-5);
  if(match($9,/GeneID:[0-9]+/)) gnum=substr($9,RSTART+7,RLENGTH-7);
  if(gname!="" && gnum!="") print "LOC"gnum"\t"gname
}' "${GFF}" > "${OUTFILE}"

echo "[module07b] 완료: ${OUTFILE} ($(wc -l < "${OUTFILE}")개)"

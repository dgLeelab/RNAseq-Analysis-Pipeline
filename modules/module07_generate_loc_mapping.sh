#!/usr/bin/env bash
# =========================================================
# module07_generate_loc_mapping.sh
# gene_id -> LOC(NCBI GeneID) 매핑표 생성
#
# annotation 파일의 attribute 포맷을 자동 감지합니다:
#   - GTF 스타일:  gene_id "값";
#   - GFF3 스타일: gene_id=값;   또는  ID=gene-값;...Dbxref=GeneID:숫자
# 어느 쪽인지 module00 출력을 보고 자동으로 정규식을 선택하므로,
# GTF/GFF3 어느 것을 넣어도 동작합니다.
#
# 사용: ./module07_generate_loc_mapping.sh <module00 전처리된 nuclear annotation> <출력 tsv>
# 예:   ./module07_generate_loc_mapping.sh Gmax_v4.0_nuclear.gff gene_id2loc.tsv
# =========================================================
set -euo pipefail

GFF="${1:?사용법: $0 <nuclear.gff|.gtf> <output gene_id2loc.tsv>}"
OUTFILE="${2:-gene_id2loc.tsv}"

# 첫 gene/mRNA/transcript feature 라인 하나로 포맷 판별
SAMPLE_LINE=$(awk -F'\t' '$3=="gene" || $3=="mRNA" || $3=="transcript" {print; exit}' "${GFF}")

if echo "${SAMPLE_LINE}" | grep -q 'gene_id "'; then
    FORMAT="GTF"
elif echo "${SAMPLE_LINE}" | grep -q 'gene_id='; then
    FORMAT="GFF3_gene_id_eq"
elif echo "${SAMPLE_LINE}" | grep -q 'ID=gene-'; then
    FORMAT="GFF3_ID_prefix"
else
    echo "[module07] WARNING: attribute 포맷을 자동 감지하지 못했습니다."
    echo "  샘플 라인: ${SAMPLE_LINE}"
    echo "  기본값(GFF3_gene_id_eq)으로 시도합니다. 결과(${OUTFILE})가 비어있으면"
    echo "  이 스크립트의 정규식을 실제 포맷에 맞게 직접 수정하세요."
    FORMAT="GFF3_gene_id_eq"
fi

echo "[module07] 감지된 annotation 포맷: ${FORMAT}"

case "${FORMAT}" in
    GTF)
        awk -F'\t' '$3=="gene"{
          gid=""; gnum="";
          if(match($9,/gene_id "[^"]+"/)) { gid=substr($9,RSTART+9,RLENGTH-10) }
          if(match($9,/GeneID:[0-9]+/)) gnum=substr($9,RSTART+7,RLENGTH-7);
          if(gid!="" && gnum!="") print gid"\tLOC"gnum
        }' "${GFF}" > "${OUTFILE}"
        ;;
    GFF3_gene_id_eq)
        awk -F'\t' '$3=="gene"{
          gid=""; gnum="";
          if(match($9,/gene_id=[^;]+/)) gid=substr($9,RSTART+8,RLENGTH-8);
          if(match($9,/GeneID:[0-9]+/)) gnum=substr($9,RSTART+7,RLENGTH-7);
          if(gid!="" && gnum!="") print gid"\tLOC"gnum
        }' "${GFF}" > "${OUTFILE}"
        ;;
    GFF3_ID_prefix)
        awk -F'\t' '$3=="gene"{
          gid=""; gnum="";
          if(match($9,/ID=gene-[^;]+/)) gid=substr($9,RSTART+8,RLENGTH-8);
          if(match($9,/GeneID:[0-9]+/)) gnum=substr($9,RSTART+7,RLENGTH-7);
          if(gid!="" && gnum!="") print gid"\tLOC"gnum
        }' "${GFF}" > "${OUTFILE}"
        ;;
esac

echo "[module07] 완료: ${OUTFILE}"
echo "[module07] 매핑 개수: $(wc -l < "${OUTFILE}")"

if [[ ! -s "${OUTFILE}" ]]; then
    echo "[module07] ERROR: 매핑 결과가 비어있습니다! annotation 포맷 자동감지가 실패했을 수 있습니다."
    echo "  샘플 라인을 직접 확인하세요: awk -F'\t' '\$3==\"gene\"{print; exit}' ${GFF}"
    exit 1
fi

echo "[module07] 중복 gene_id 확인 (비어야 정상):"
cut -f1 "${OUTFILE}" | sort | uniq -d

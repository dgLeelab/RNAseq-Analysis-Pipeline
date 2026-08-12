#!/usr/bin/env bash
# =========================================================
# list_raw_prefixes.sh
# raw fastq 디렉토리에서 R1 파일들의 prefix(샘플 구분자) 목록을 뽑아줍니다.
# 이 목록을 보고 sample_map.tsv (sample_id <TAB> raw_prefix)를 채우면 됩니다.
#
# 사용: ./list_raw_prefixes.sh
# =========================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=========================================="
echo "RAW_DIR: ${RAW_DIR}"
echo "=========================================="
echo ""
echo "R1 파일 목록 (prefix 추정 포함):"
echo ""

for f in "${RAW_DIR}"/*_R1_*.fastq.gz "${RAW_DIR}"/*_R1.fastq.gz; do
    [[ -f "${f}" ]] || continue
    fname=$(basename "${f}")
    # _S숫자_L숫자_R1_001.fastq.gz 형태에서 prefix만 추출
    prefix=$(echo "${fname}" | sed -E 's/_S[0-9]+_L[0-9]+_R1_[0-9]+\.fastq\.gz$//; s/_R1\.fastq\.gz$//')
    printf "  %-40s -> prefix: %s\n" "${fname}" "${prefix}"
done

echo ""
echo "=========================================="
echo "sample_map.tsv 뼈대 (${PROJECT_DIR}/sample_map.tsv 로 저장 후 sample_id 칼럼을 채우세요)"
echo "=========================================="
echo ""
echo -e "# sample_id\traw_prefix"
for f in "${RAW_DIR}"/*_R1_*.fastq.gz "${RAW_DIR}"/*_R1.fastq.gz; do
    [[ -f "${f}" ]] || continue
    fname=$(basename "${f}")
    prefix=$(echo "${fname}" | sed -E 's/_S[0-9]+_L[0-9]+_R1_[0-9]+\.fastq\.gz$//; s/_R1\.fastq\.gz$//')
    echo -e "<CHANGE_ME>\t${prefix}"
done

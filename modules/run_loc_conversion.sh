#!/usr/bin/env bash
# =========================================================
# run_loc_conversion.sh
# module07(gene_id2loc.tsv) 이후: count matrix, TPM matrix를 각각
# LOC 식별자로 치환하고, 이후 단계가 쓸 최종 파일명으로 정리합니다.
#
#   gene_count_matrix.csv        -> genes.counts.matrix (LOC, TSV)
#   expression.TPM.matrix        -> expression.TPM.matrix (LOC로 덮어씀, 원본은 .orig 보관)
#
# 사용: ./run_loc_conversion.sh [gene_id2loc.tsv] [count_csv] [tpm_matrix]
# 실행 위치: module05, module06, module07 결과가 모두 있는 디렉토리
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAP_FILE="${1:-gene_id2loc.tsv}"
COUNT_CSV="${2:-gene_count_matrix.csv}"
TPM_MATRIX="${3:-expression.TPM.matrix}"

for f in "${MAP_FILE}" "${COUNT_CSV}" "${TPM_MATRIX}"; do
    [[ -f "${f}" ]] || { echo "[run_loc_conversion] ERROR: ${f} 없음"; exit 1; }
done

echo "[1/2] count matrix LOC 치환"
R_ENVIRON_USER=/dev/null Rscript "${SCRIPT_DIR}/loc_convert.R" \
    "${MAP_FILE}" "${COUNT_CSV}" "gene_count_matrix.loc.csv"

# csv -> tsv 변환 (헤더의 R 타입 표기 등 불필요한 라인 방어적으로 제외)
grep -v "<class" gene_count_matrix.loc.csv | sed 's/,/\t/g' > genes.counts.matrix
echo "  -> genes.counts.matrix 생성"

echo "[2/2] TPM matrix LOC 치환"
R_ENVIRON_USER=/dev/null Rscript "${SCRIPT_DIR}/loc_convert_tpm.R" \
    "${MAP_FILE}" "${TPM_MATRIX}" "${TPM_MATRIX}.loc"

cp "${TPM_MATRIX}" "${TPM_MATRIX}.orig"
cp "${TPM_MATRIX}.loc" "${TPM_MATRIX}"
echo "  -> ${TPM_MATRIX} (LOC로 치환됨, 원본은 ${TPM_MATRIX}.orig 보관)"

echo "[run_loc_conversion] 완료"

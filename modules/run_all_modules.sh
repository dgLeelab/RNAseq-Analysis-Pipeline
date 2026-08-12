#!/usr/bin/env bash
# =========================================================
# run_all_modules.sh
# module00 ~ module09 전체 순서 실행 (참고용 오케스트레이터)
#
# 프로젝트마다 organelle scaffold ID, threshold(P/C) 등이 다르므로
# 실제로는 각 module을 상황 보며 개별 실행하는 것을 권장합니다.
# 이 스크립트는 "정상 케이스"의 전체 흐름을 보여주는 참고용입니다.
#
# 사용 (예시, 대두 프로젝트 기준):
#   ./run_all_modules.sh \
#       <원본genome.gff> \
#       NC_007942.1 NC_020455.1
# =========================================================
set -euo pipefail

RAW_GFF="${1:?사용법: $0 <원본 GFF3> [organelle_scaffold_id ...]}"
shift
ORGANELLE_IDS=("$@")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NUCLEAR_GFF="nuclear_annotation.gff"

echo ">>> module00: annotation 전처리"
"${SCRIPT_DIR}/module00_prep_annotation.sh" "${RAW_GFF}" "${NUCLEAR_GFF}" "${ORGANELLE_IDS[@]}"

echo ">>> module04: StringTie"
"${SCRIPT_DIR}/module04_stringtie.sh" "${NUCLEAR_GFF}"

echo ">>> module00b: 폴더명 검증 (StringTie 결과 생성된 뒤에만 의미 있음)"
"${SCRIPT_DIR}/module00b_check_folder_names.sh"

echo ">>> module05: prepDE (count matrix)"
"${SCRIPT_DIR}/run_module05.sh"

echo ">>> module06: TPM 추출"
python3 "${SCRIPT_DIR}/module06_TPM_Extraction.py" expression.TPM.matrix

echo ">>> module07: gene_id -> LOC 매핑표 생성"
"${SCRIPT_DIR}/module07_generate_loc_mapping.sh" "${NUCLEAR_GFF}" gene_id2loc.tsv

echo ">>> LOC 치환 (count matrix + TPM matrix)"
"${SCRIPT_DIR}/run_loc_conversion.sh" gene_id2loc.tsv gene_count_matrix.csv expression.TPM.matrix

echo ">>> module08: TMM normalization (raw count matrix 기준)"
"${SCRIPT_DIR}/module08_TMM_normalization.sh" genes.counts.matrix

echo ">>> samples.file 생성"
"${SCRIPT_DIR}/generate_samples_file.sh" genes.counts.matrix samples.file

echo ">>> module09: PtR + edgeR DE + clustering (기본 P=1e-3, C=2)"
"${SCRIPT_DIR}/module09_matrix_process.sh" -P 1e-3 -C 2

echo ">>> 검증"
"${SCRIPT_DIR}/module_validate.sh" genes.counts.matrix expression.TMM.matrix samples.file

echo ""
echo "전체 파이프라인 완료."

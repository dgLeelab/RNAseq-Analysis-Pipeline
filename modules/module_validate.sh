#!/usr/bin/env bash
# =========================================================
# module_validate.sh
# 최종 matrix들에 대한 4가지 핵심 검증 (모두 "출력 없음/0" 이 정상)
#
# 사용: ./module_validate.sh [genes.counts.matrix] [expression.TMM.matrix] [samples.file]
# =========================================================
set -uo pipefail

COUNT_MATRIX="${1:-genes.counts.matrix}"
TMM_MATRIX="${2:-expression.TMM.matrix}"
SAMPLES_FILE="${3:-samples.file}"

echo "===================================================="
echo "1) 식별자 중복 확인 (비어야 정상)"
echo "===================================================="
tail -n +2 "${COUNT_MATRIX}" | cut -f1 | sort | uniq -d

echo ""
echo "===================================================="
echo "2) 컬럼 수 일치 확인 (출력 없어야 정상)"
echo "===================================================="
awk -F'\t' 'NR==1{n=NF} NF!=n{print NR,NF,n}' "${COUNT_MATRIX}"

echo ""
echo "===================================================="
echo "3) counts ⊆ TMM 확인 (반드시 0, 가장 중요)"
echo "===================================================="
tail -n +2 "${COUNT_MATRIX}" | cut -f1 | LC_ALL=C sort > /tmp/module_validate_c.txt
tail -n +2 "${TMM_MATRIX}"   | cut -f1 | LC_ALL=C sort > /tmp/module_validate_t.txt
diff_count=$(comm -23 /tmp/module_validate_c.txt /tmp/module_validate_t.txt | wc -l)
echo "counts에만 있고 TMM에 없는 유전자 수: ${diff_count} (0 이어야 정상)"
rm -f /tmp/module_validate_c.txt /tmp/module_validate_t.txt

echo ""
echo "===================================================="
echo "4) samples.file 샘플명 == matrix 헤더 컬럼명 일치 확인"
echo "===================================================="
matrix_cols=$(head -1 "${COUNT_MATRIX}" | tr '\t' '\n' | tail -n +2 | sort)
samples_col2=$(cut -f2 "${SAMPLES_FILE}" | sort)
if [[ "${matrix_cols}" == "${samples_col2}" ]]; then
    echo "OK: 일치함"
else
    echo "MISMATCH 발견:"
    diff <(echo "${matrix_cols}") <(echo "${samples_col2}")
fi

echo ""
echo "검증 완료."

#!/usr/bin/env bash
# =========================================================
# generate_samples_file.sh
# genes.counts.matrix의 헤더(컬럼명)로부터 samples.file(group<TAB>sample)을
# 직접 생성합니다. matrix 컬럼명과 100% 일치를 보장하기 위해, samples.txt가
# 아니라 실제 matrix 헤더에서 파싱합니다 (실수로 이름이 어긋나는 것 방지).
#
# 컬럼명이 "{group}_{replicate}" 형식이라고 가정하고, 마지막 "_숫자"를
# replicate로, 그 앞부분 전체를 group으로 사용합니다.
#
# 사용: ./generate_samples_file.sh [genes.counts.matrix] [출력 samples.file]
# =========================================================
set -uo pipefail

MATRIX="${1:-genes.counts.matrix}"
OUT_FILE="${2:-samples.file}"

if [[ ! -f "${MATRIX}" ]]; then
    echo "[generate_samples_file] ERROR: ${MATRIX} 없음"
    exit 1
fi

head -1 "${MATRIX}" | tr '\t' '\n' | tail -n +2 | \
    awk -F'_' '{grp=$1; for(i=2;i<NF;i++) grp=grp"_"$i; print grp"\t"$0}' > "${OUT_FILE}"

echo "[generate_samples_file] 생성 완료: ${OUT_FILE}"
cat "${OUT_FILE}"
echo ""
echo "[generate_samples_file] 그룹핑 확인:"
cut -f1 "${OUT_FILE}" | sort | uniq -c

echo ""
echo "[generate_samples_file] 명명 규칙 검증 (컬럼명이 {group}_{replicate숫자} 형식인지):"
bad=0
while IFS=$'\t' read -r grp sample; do
    suffix="${sample##*_}"
    if ! [[ "${suffix}" =~ ^[0-9]+$ ]]; then
        echo "  [WARNING] ${sample} : replicate 접미사가 숫자가 아닙니다 (그룹핑이 잘못됐을 수 있음)"
        bad=1
    fi
done < "${OUT_FILE}"
if [[ "${bad}" -eq 0 ]]; then
    echo "  OK: 모든 샘플명이 {group}_{숫자} 형식입니다."
else
    echo "  -> samples.file을 열어 그룹핑이 의도한 대로 됐는지 직접 확인하세요."
fi

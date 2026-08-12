#!/usr/bin/env bash
# =========================================================
# module08_TMM_normalization.sh
#
# genes.counts.matrix (run_loc_conversion.sh 결과, raw count) -> expression.TMM.matrix
# Trinity의 run_TMM_scale_matrix.pl 사용 (PATH에 있어야 함, trinity_env 활성화 필요)
#
# run_TMM_scale_matrix.pl --matrix 는 raw read count matrix를 요구합니다
# (TPM/FPKM 등 이미 정규화된 값을 넣으면 안 됨). TPM.matrix를 넣으면
# lib.size가 모든 샘플에서 동일(1e6)해져 TMM 계수가 왜곡됩니다.
#
# 실행 위치: run_loc_conversion.sh 결과(genes.counts.matrix)가 있는 디렉토리
# 사용: ./module08_TMM_normalization.sh [입력 raw count matrix, 기본: genes.counts.matrix]
# =========================================================
set -uo pipefail

COUNT_MATRIX="${1:-genes.counts.matrix}"
OUT_MATRIX="expression.TMM.matrix"

if [[ ! -f "${COUNT_MATRIX}" ]]; then
    echo "[module08] ERROR: ${COUNT_MATRIX} 없음. run_loc_conversion.sh를 먼저 실행하세요."
    exit 1
fi

if ! command -v run_TMM_scale_matrix.pl &> /dev/null; then
    echo "[module08] ERROR: run_TMM_scale_matrix.pl 을 PATH에서 찾을 수 없습니다."
    echo "           trinity_env를 활성화했는지 확인하세요:"
    echo "           conda activate trinity_env"
    exit 1
fi

echo "[module08] TMM normalization 시작: ${COUNT_MATRIX} -> ${OUT_MATRIX}"

# run_TMM_scale_matrix.pl은 정규화된 매트릭스를 STDOUT으로 출력하고,
# <matrix>.TMM_info.txt(정규화 계수)와 <matrix>.runTMM.R을 부산물로 남긴다.
run_TMM_scale_matrix.pl --matrix "${COUNT_MATRIX}" > "${OUT_MATRIX}" 2> module08_tmm.log

if [[ -s "${OUT_MATRIX}" ]]; then
    echo "[module08] 완료: ${OUT_MATRIX} (정규화 계수: ${COUNT_MATRIX}.TMM_info.txt)"
else
    echo "[module08] ERROR: ${OUT_MATRIX}가 생성되지 않았습니다."
    echo "           module08_tmm.log 확인 필요"
    exit 1
fi

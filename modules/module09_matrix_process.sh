#!/usr/bin/env bash

# =========================================================
# module09_matrix_process.sh
#
# PtR(샘플 상관/PCA QC) + run_DE_analysis.pl(edgeR pairwise DE)
# + analyze_diff_expr.pl(TMM matrix 기준 클러스터링용 union DE 추출)
# + define_clusters_by_cutting_tree.pl
#
# 필요: Trinity PtR, run_DE_analysis.pl, analyze_diff_expr.pl,
# define_clusters_by_cutting_tree.pl (PATH, trinity_env)
#
# 필요 입력 파일 (같은 디렉토리):
# - genes.counts.matrix
# - expression.TMM.matrix
# - samples.file
#
# 사용:
# ./module09_matrix_process.sh [-P Pvalue] [-C FoldChange]
#
# 예:
# ./module09_matrix_process.sh -P 1e-3 -C 2
#
# 결과:
# P0.001_FC2/
# =========================================================

set -uo pipefail


# ---------------------------------------------------------
# 입력 파일
# ---------------------------------------------------------

COUNT_MATRIX="genes.counts.matrix"
TMM_MATRIX="expression.TMM.matrix"
SAMPLES_FILE="samples.file"


# ---------------------------------------------------------
# threshold
# ---------------------------------------------------------

PVALUE="1e-3"
FOLDCHANGE="2"

while getopts "P:C:" opt; do
    case "${opt}" in
        P) PVALUE="${OPTARG}" ;;
        C) FOLDCHANGE="${OPTARG}" ;;
        *) ;;
    esac
done


# ---------------------------------------------------------
# 입력 파일 확인
# ---------------------------------------------------------

for f in "${COUNT_MATRIX}" "${TMM_MATRIX}" "${SAMPLES_FILE}"; do
    if [[ ! -f "${f}" ]]; then
        echo "[module09] ERROR: ${f} 없음."

        if [[ "${f}" == "${SAMPLES_FILE}" ]]; then
            echo "  generate_samples_file.sh 를 먼저 실행하세요."
        fi

        exit 1
    fi
done


# =========================================================
# ★ 결과 폴더 설정
# =========================================================

WORKDIR="$(pwd)"

# 1e-3 -> 0.001
# 5e-2 -> 0.05
PVALUE_LABEL=$(awk -v p="${PVALUE}" 'BEGIN { printf "%.10g", p+0 }')

OUTDIR="P${PVALUE_LABEL}_FC${FOLDCHANGE}"

echo "[module09] 결과 폴더: ${OUTDIR}"

mkdir -p "${OUTDIR}"


# cd 이후에도 입력 파일을 읽을 수 있도록 절대경로 지정
COUNT_MATRIX="${WORKDIR}/${COUNT_MATRIX}"
TMM_MATRIX="${WORKDIR}/${TMM_MATRIX}"
SAMPLES_FILE="${WORKDIR}/${SAMPLES_FILE}"

cd "${OUTDIR}" || exit 1


# ---------------------------------------------------------
# R 환경
# ---------------------------------------------------------

export R_ENVIRON_USER=/dev/null


# =========================================================
# 1. PtR
# =========================================================

echo "[module09-1] PtR: 샘플 상관관계 + PCA QC (raw count matrix 기준)"

PtR \
    --matrix "${COUNT_MATRIX}" \
    --samples "${SAMPLES_FILE}" \
    --log2 \
    --min_rowSums 10 \
    --compare_replicates \
    --CPM \
    --sample_cor_matrix \
    --center_rows \
    --prin_comp 3 \
    > module09_ptr.log 2>&1

echo "  -> module09_ptr.log 확인 (PCA≈PC1/PC2 분산비, replicate 응집도 등)"


# =========================================================
# 2. edgeR pairwise DE
# =========================================================

echo "[module09-2] run_DE_analysis.pl (edgeR, 품종 pairwise DE)"

run_DE_analysis.pl \
    --matrix "${COUNT_MATRIX}" \
    --method edgeR \
    --samples_file "${SAMPLES_FILE}" \
    > module09_de_pairwise.log 2>&1


DE_DIR=$(ls -d edgeR* 2>/dev/null | head -1)

if [[ -z "${DE_DIR}" ]]; then
    echo "[module09] ERROR: run_DE_analysis.pl 결과 디렉토리(edgeR*)를 찾지 못했습니다."
    echo "           module09_de_pairwise.log 확인 필요"
    exit 1
fi


cp "${DE_DIR}"/*.DE_results . 2>/dev/null || true

echo "  -> ${DE_DIR}/*.DE_results 를 현재 디렉토리로 복사함 (66쌍 등 pairwise 결과)"


# =========================================================
# 3. union DE
# =========================================================

echo "[module09-3] analyze_diff_expr.pl (TMM matrix 기준 union DE 추출, -P ${PVALUE} -C ${FOLDCHANGE})"

analyze_diff_expr.pl \
    --max_genes_clust 15000 \
    --matrix "${TMM_MATRIX}" \
    --samples "${SAMPLES_FILE}" \
    -P "${PVALUE}" \
    -C "${FOLDCHANGE}" \
    > module09_analyze.log 2>&1


# =========================================================
# 4. RData 자동 감지
# =========================================================

echo "[module09-4] RData 자동 감지 (파일명이 과학적 표기법일 수 있음, 하드코딩 금지)"

RDATA=$(ls diffExpr.*.matrix.RData 2>/dev/null | head -1)

if [[ -z "${RDATA}" ]]; then
    echo "[module09] ERROR: diffExpr.*.matrix.RData 파일이 생성되지 않았습니다."
    echo "           -P/-C 임계값이 너무 엄격해서 유의 유전자가 0개일 수 있습니다."
    echo "           (본 프로젝트 기준 권장값: P<0.001, C2)"
    echo "           -> module09_analyze.log 확인 필요"
    echo "           마지막 줄이 'Reading matrix file.'에서 멈춰있으면 threshold 완화 필요"
    exit 1
fi

echo "  [감지됨] ${RDATA}"


# =========================================================
# 5. cluster
# =========================================================

echo "[module09-5] define_clusters_by_cutting_tree.pl"

define_clusters_by_cutting_tree.pl \
    -R "${RDATA}" \
    --Ptree 30 \
    --no_column_reordering \
    > module09_cluster.log 2>&1


# =========================================================
# 완료
# =========================================================

echo
echo "[module09] 완료."
echo
echo "  결과 디렉토리:"
echo "    ${WORKDIR}/${OUTDIR}"
echo
echo "  주요 산출물:"
echo "    - module09_ptr.log 관련 PDF (샘플 상관/PCA)"
echo "    - *.edgeR.DE_results (품종쌍 DE 결과)"
echo "    - ${RDATA}"
echo "    - diffExpr.*.genes_vs_samples_heatmap.pdf"
echo "    - clusters_fixed_P_30.*"

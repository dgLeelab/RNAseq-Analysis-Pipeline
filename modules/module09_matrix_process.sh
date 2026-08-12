#!/usr/bin/env bash
# =========================================================
# module09_matrix_process.sh
# PtR(샘플 상관/PCA QC) + run_DE_analysis.pl(edgeR pairwise DE)
# + analyze_diff_expr.pl(TMM matrix 기준 클러스터링용 union DE 추출)
# + define_clusters_by_cutting_tree.pl
#
# 필요: Trinity PtR, run_DE_analysis.pl, analyze_diff_expr.pl,
#       define_clusters_by_cutting_tree.pl (PATH, trinity_env)
# 필요 입력 파일 (같은 디렉토리):
#   - genes.counts.matrix     (raw count, LOC 식별자, run_loc_conversion.sh 결과)
#   - expression.TPM.matrix   (TMM 정규화 전 TPM, LOC 식별자로 이미 치환됨)
#   - samples.file            (generate_samples_file.sh 결과)
#
# ★★★ 가장 중요한 함정 (실제로 겪은 것) ★★★
#   1) analyze_diff_expr.pl 에 -P/-C 를 명시하지 않으면 기본값
#      (P0.001_C2)으로 실행되어, 폴더명/의도한 threshold와 실제 생성
#      파일명이 어긋날 수 있음.
#   2) 이 Trinity 버전은 RData 파일명에 과학적 표기법을 사용합니다:
#      -P 1e-3 -C 2  ->  diffExpr.P1e-3_C2.matrix.RData  (P0.001_C2 아님!)
#      따라서 define_clusters -R 값은 절대 하드코딩하지 말고 자동 감지.
#   3) threshold가 너무 엄격하면(P<0.0001 & C3 등) 유의 유전자가 0개가
#      되어 RData 자체가 생성되지 않음. 로그에 "Reading matrix file."
#      까지만 찍히고 조용히 끝나면 이 경우이니 threshold를 완화할 것.
#      본 프로젝트(같은 종 품종 간 비교) 기준 P<0.001, C2(FC>=2)가 적정.
#
# .Renviron의 R_LIBS가 R의 S7 패키지와 충돌하므로 R_ENVIRON_USER=/dev/null 필수.
#
# 사용: ./module09_matrix_process.sh [-P Pvalue] [-C FoldChange]
#   예) ./module09_matrix_process.sh -P 1e-3 -C 2
# =========================================================
set -uo pipefail

COUNT_MATRIX="genes.counts.matrix"
TMM_MATRIX="expression.TMM.matrix"
SAMPLES_FILE="samples.file"

PVALUE="1e-3"
FOLDCHANGE="2"
while getopts "P:C:" opt; do
    case "${opt}" in
        P) PVALUE="${OPTARG}" ;;
        C) FOLDCHANGE="${OPTARG}" ;;
        *) ;;
    esac
done

for f in "${COUNT_MATRIX}" "${TMM_MATRIX}" "${SAMPLES_FILE}"; do
    if [[ ! -f "${f}" ]]; then
        echo "[module09] ERROR: ${f} 없음."
        [[ "${f}" == "${SAMPLES_FILE}" ]] && \
            echo "  generate_samples_file.sh 를 먼저 실행하세요."
        exit 1
    fi
done

export R_ENVIRON_USER=/dev/null

echo "[module09-1] PtR: 샘플 상관관계 + PCA QC (raw count matrix 기준)"
PtR --matrix "${COUNT_MATRIX}" \
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

echo "[module09-3] analyze_diff_expr.pl (TMM matrix 기준 union DE 추출, -P ${PVALUE} -C ${FOLDCHANGE})"
analyze_diff_expr.pl \
    --max_genes_clust 15000 \
    --matrix "${TMM_MATRIX}" \
    --samples "${SAMPLES_FILE}" \
    -P "${PVALUE}" -C "${FOLDCHANGE}" \
    > module09_analyze.log 2>&1

echo "[module09-4] RData 자동 감지 (파일명이 과학적 표기법일 수 있음, 하드코딩 금지)"
RDATA=$(ls diffExpr.*.matrix.RData 2>/dev/null | head -1)

if [[ -z "${RDATA}" ]]; then
    echo "[module09] ERROR: diffExpr.*.matrix.RData 파일이 생성되지 않았습니다."
    echo "           -P/-C 임계값이 너무 엄격해서 유의 유전자가 0개일 수 있습니다."
    echo "           (본 프로젝트 기준 권장값: P<0.001, C2)"
    echo "           -> module09_analyze.log 확인 필요 (마지막 줄이 'Reading matrix file.'"
    echo "              에서 멈춰있으면 threshold 완화 필요)"
    exit 1
fi
echo "  [감지됨] ${RDATA}"

echo "[module09-5] define_clusters_by_cutting_tree.pl"
define_clusters_by_cutting_tree.pl \
    -R "${RDATA}" \
    --Ptree 30 \
    --no_column_reordering \
    > module09_cluster.log 2>&1

echo "[module09] 완료."
echo "  주요 산출물:"
echo "    - module09_ptr.log 관련 PDF (샘플 상관/PCA)"
echo "    - *.edgeR.DE_results (품종쌍 DE 결과)"
echo "    - ${RDATA}"
echo "    - diffExpr.*.genes_vs_samples_heatmap.pdf"
echo "    - clusters_fixed_P_30.* (유전자 클러스터)"

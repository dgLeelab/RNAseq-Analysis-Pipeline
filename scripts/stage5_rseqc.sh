#!/usr/bin/env bash
# =========================================================
# Stage 5: Alignment QC (RSeQC) + Strand 자동판별
# 사용 env: qc_tools_env (rseqc)
#
# 처리 순서:
#   1) strand_info.sh가 없으면, 첫 번째 유효 샘플만 먼저 "순차" 처리해서
#      strand를 판별하고 strand_info.sh 저장
#   2) 나머지 샘플들은 PARALLEL_SAMPLES_LIGHT개씩 동시 병렬 처리
#   3) geneBody_coverage.py는 전체 샘플 합산이라 병렬화 대상 아님 (마지막에 1회 실행)
#
# 이미 결과 파일이 있는 샘플은 자동 SKIP.
# =========================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env "${ENV_QC_TOOLS}"

log "Stage 5: RSeQC 시작 (동시 ${PARALLEL_SAMPLES_LIGHT}개 샘플)"

: > "${DIR_RSEQC}/bam_list_for_genebody.txt"

# ---- 샘플 1개 처리 함수 (infer_experiment + read_distribution) ----
rseqc_one_sample () {
    local sample="$1"
    local bam="${DIR_ALIGN}/${sample}.sorted.bam"
    local infer_out="${DIR_RSEQC}/${sample}.infer_experiment.txt"
    local dist_out="${DIR_RSEQC}/${sample}.read_distribution.txt"

    if [[ -f "${infer_out}" && -f "${dist_out}" ]]; then
        log "  [ALREADY DONE] ${sample} -> skip"
        [[ -f "${bam}" ]] && echo "${bam}" >> "${DIR_RSEQC}/bam_list_for_genebody.txt"
        return 0
    fi

    if [[ ! -f "${bam}" ]]; then
        log "  [SKIP] ${sample}: sorted bam 없음"
        return 0
    fi

    log "  [infer_experiment] ${sample}"
    infer_experiment.py -i "${bam}" -r "${RSEQC_BED}" \
        > "${infer_out}" 2>> "${LOG_DIR}/stage5_rseqc.log"

    log "  [read_distribution] ${sample}"
    read_distribution.py -i "${bam}" -r "${RSEQC_BED}" \
        > "${dist_out}" 2>> "${LOG_DIR}/stage5_rseqc.log"

    echo "${bam}" >> "${DIR_RSEQC}/bam_list_for_genebody.txt"
    log "  [완료] ${sample}"
}
export -f rseqc_one_sample
export DIR_ALIGN DIR_RSEQC RSEQC_BED LOG_DIR

determine_strand () {
    local sample="$1"
    local infer_file="${DIR_RSEQC}/${sample}.infer_experiment.txt"

    local fwd_frac rev_frac
    fwd_frac=$(grep -oP '(?<=\+\+,1--,2\+-,2-\+"\: )[0-9.]+' "${infer_file}" || echo "0")
    rev_frac=$(grep -oP '(?<=\+-,1-\+,2\+\+,2--"\: )[0-9.]+' "${infer_file}" || echo "0")

    log "  [Strand 판별] ${sample} 기준 forward=${fwd_frac}, reverse=${rev_frac}"

    local is_fwd is_rev
    is_fwd=$(awk -v f="${fwd_frac}" 'BEGIN{print (f>0.6)?1:0}')
    is_rev=$(awk -v r="${rev_frac}" 'BEGIN{print (r>0.6)?1:0}')

    local STRANDEDNESS HISAT2_STRAND_OPT STRINGTIE_STRAND_OPT DUPRADAR_STRANDED
    if [[ "${is_fwd}" == "1" ]]; then
        STRANDEDNESS="forward"; HISAT2_STRAND_OPT="--rna-strandness FR"
        STRINGTIE_STRAND_OPT="--fr"; DUPRADAR_STRANDED=1
    elif [[ "${is_rev}" == "1" ]]; then
        STRANDEDNESS="reverse"; HISAT2_STRAND_OPT="--rna-strandness RF"
        STRINGTIE_STRAND_OPT="--rf"; DUPRADAR_STRANDED=2
    else
        STRANDEDNESS="unstranded"; HISAT2_STRAND_OPT=""
        STRINGTIE_STRAND_OPT=""; DUPRADAR_STRANDED=0
    fi

    cat > "${STRAND_INFO_FILE}" << STRANDEOF
# 자동 생성됨 (Stage 5, 기준 샘플: ${sample})
STRANDEDNESS="${STRANDEDNESS}"
HISAT2_STRAND_OPT="${HISAT2_STRAND_OPT}"
STRINGTIE_STRAND_OPT="${STRINGTIE_STRAND_OPT}"
DUPRADAR_STRANDED=${DUPRADAR_STRANDED}
STRANDEOF

    log "  [Strand 판별 결과] ${STRANDEDNESS} -> ${STRAND_INFO_FILE} 저장됨"
}

# ---- 1단계: strand_info.sh 없으면 첫 유효 샘플만 순차 처리해서 판별 ----
FIRST_SAMPLE=""
if [[ ! -f "${STRAND_INFO_FILE}" ]]; then
    while read -r sample; do
        [[ -z "${sample}" ]] && continue
        bam="${DIR_ALIGN}/${sample}.sorted.bam"
        [[ -f "${bam}" ]] || continue
        FIRST_SAMPLE="${sample}"
        break
    done < "${SAMPLE_LIST}"

    if [[ -z "${FIRST_SAMPLE}" ]]; then
        log "  [ERROR] Strand 판별할 샘플이 하나도 없습니다. Stage 4(align) 결과를 확인하세요."
        deactivate_env
        exit 1
    fi

    log "Stage 5-1: strand 판별용 기준 샘플 순차 처리: ${FIRST_SAMPLE}"
    rseqc_one_sample "${FIRST_SAMPLE}"
    determine_strand "${FIRST_SAMPLE}"
else
    log "  [Strand 이미 판별됨] $(grep STRANDEDNESS "${STRAND_INFO_FILE}") -> 재판별 skip"
fi

# ---- 2단계: 나머지 샘플 병렬 처리 ----
log "Stage 5-2: 나머지 샘플 병렬 처리 시작"
job_count=0
while read -r sample; do
    [[ -z "${sample}" ]] && continue
    [[ "${sample}" == "${FIRST_SAMPLE}" ]] && continue   # 이미 처리한 기준 샘플 제외

    rseqc_one_sample "${sample}" &
    job_count=$((job_count + 1))

    if (( job_count >= PARALLEL_SAMPLES_LIGHT )); then
        wait -n
        job_count=$((job_count - 1))
    fi

done < "${SAMPLE_LIST}"

wait   # 남은 job 전부 대기

# ---- geneBody_coverage (전체 샘플 합산, 병렬화 대상 아님) ----
log "  [geneBody_coverage] 전체 샘플 합산 실행"
geneBody_coverage.py \
    -i "$(paste -sd, "${DIR_RSEQC}/bam_list_for_genebody.txt")" \
    -r "${RSEQC_BED}" \
    -o "${DIR_RSEQC}/all_samples" \
    >> "${LOG_DIR}/stage5_rseqc.log" 2>&1

log "Stage 5: MultiQC 통합"
multiqc "${DIR_RSEQC}" -o "${DIR_RSEQC}" -n multiqc_rseqc --force \
    >> "${LOG_DIR}/stage5_multiqc.log" 2>&1

deactivate_env

log "Stage 5 완료 (strand: $(grep STRANDEDNESS "${STRAND_INFO_FILE}"))"

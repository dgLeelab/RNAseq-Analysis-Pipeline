#!/usr/bin/env bash
# =========================================================
# Stage 4: Alignment (HISAT2 -> sorted BAM via samtools)
# 사용 env: 시스템 PATH (hisat2), base (samtools)
#
# PARALLEL_SAMPLES 개씩 동시에 HISAT2를 돌립니다 (샘플당 THREADS_PER_SAMPLE 스레드).
# 순차 처리보다 서버 코어를 훨씬 효율적으로 씁니다.
#
# 이미 결과(.sorted.bam)가 있는 샘플은 자동 skip.
# =========================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env "base"

HISAT2_STRAND_OPT=""
if [[ -f "${STRAND_INFO_FILE}" ]]; then
    # shellcheck disable=SC1090
    source "${STRAND_INFO_FILE}"
    log "Stage 4: 기존 strand 정보 감지 (${STRANDEDNESS}) -> HISAT2에 반영"
fi

log "Stage 4: HISAT2 alignment 시작 (동시 ${PARALLEL_SAMPLES}개 샘플, 샘플당 ${THREADS_PER_SAMPLE}스레드)"

# ---- 샘플 1개 처리 함수 (백그라운드 job으로 실행됨) ----
align_one_sample () {
    local sample="$1"
    local r1="${DIR_TRIM}/${sample}_R1.trim.fastq.gz"
    local r2="${DIR_TRIM}/${sample}_R2.trim.fastq.gz"
    local sam_out="${DIR_ALIGN}/${sample}.sam"
    local bam_sorted="${DIR_ALIGN}/${sample}.sorted.bam"
    local summary="${DIR_ALIGN}/${sample}.hisat2_summary.txt"

    if [[ -f "${bam_sorted}" ]]; then
        log "  [ALREADY DONE] ${sample} -> skip"
        return 0
    fi

    if [[ ! -f "${r1}" || ! -f "${r2}" ]]; then
        log "  [SKIP] ${sample}: trimmed fastq 없음"
        return 0
    fi

    log "  [HISAT2 시작] ${sample}"
    hisat2 \
        -p "${THREADS_PER_SAMPLE}" \
        -x "${HISAT2_INDEX_PREFIX}" \
        -1 "${r1}" -2 "${r2}" \
        ${HISAT2_STRAND_OPT} \
        --summary-file "${summary}" \
        -S "${sam_out}" \
        >> "${LOG_DIR}/stage4_hisat2.log" 2>&1

    log "  [samtools sort/index 시작] ${sample}"
    samtools sort -@ "${THREADS_PER_SAMPLE}" -o "${bam_sorted}" "${sam_out}" \
        >> "${LOG_DIR}/stage4_samtools.log" 2>&1
    samtools index "${bam_sorted}"

    rm -f "${sam_out}"
    log "  [완료] ${sample}"
}
export -f align_one_sample
export DIR_TRIM DIR_ALIGN HISAT2_INDEX_PREFIX THREADS_PER_SAMPLE HISAT2_STRAND_OPT LOG_DIR

# ---- 병렬 실행 (동시 실행 개수를 PARALLEL_SAMPLES로 제한) ----
job_count=0
while read -r sample; do
    [[ -z "${sample}" ]] && continue

    align_one_sample "${sample}" &
    job_count=$((job_count + 1))

    if (( job_count >= PARALLEL_SAMPLES )); then
        wait -n   # 실행 중인 job 중 하나가 끝날 때까지 대기
        job_count=$((job_count - 1))
    fi

done < "${SAMPLE_LIST}"

wait   # 남은 job들 전부 끝날 때까지 대기

log "Stage 4: MultiQC 통합"
multiqc "${DIR_ALIGN}" -o "${DIR_ALIGN}" -n multiqc_align --force \
    >> "${LOG_DIR}/stage4_multiqc.log" 2>&1

deactivate_env

log "Stage 4 완료"

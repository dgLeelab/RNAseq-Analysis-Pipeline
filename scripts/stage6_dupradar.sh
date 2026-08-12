#!/usr/bin/env bash
# =========================================================
# Stage 6: Duplication check
#   6-1) samtools markdup (duplicate 마킹, base env)
#   6-2) dupRadar (R, dupradar_env) - 발현량 대비 dup 해석
#        strand 값은 Stage 5에서 자동판별된 STRAND_INFO_FILE을 사용
# =========================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

if [[ ! -f "${STRAND_INFO_FILE}" ]]; then
    log "  [ERROR] ${STRAND_INFO_FILE} 없음. Stage 5(RSeQC)를 먼저 실행하세요."
    exit 1
fi
# shellcheck disable=SC1090
source "${STRAND_INFO_FILE}"
log "Stage 6: strand=${STRANDEDNESS} (dupRadar stranded=${DUPRADAR_STRANDED}) 적용"

# ---- 6-1) samtools markdup (base) ----
activate_env "base"
log "Stage 6-1: samtools markdup 시작"

while read -r sample; do
    [[ -z "${sample}" ]] && continue
    bam_in="${DIR_ALIGN}/${sample}.sorted.bam"
    bam_namesorted="${DIR_DUP}/${sample}.namesorted.bam"
    bam_fixmate="${DIR_DUP}/${sample}.fixmate.bam"
    bam_possorted="${DIR_DUP}/${sample}.possorted.bam"
    bam_markdup="${DIR_DUP}/${sample}.markdup.bam"

    if [[ ! -f "${bam_in}" ]]; then
        log "  [SKIP] ${sample}: sorted bam 없음"
        continue
    fi

    log "  [markdup] ${sample}"
    samtools sort -@ "${THREADS}" -n -o "${bam_namesorted}" "${bam_in}"
    samtools fixmate -@ "${THREADS}" -m "${bam_namesorted}" "${bam_fixmate}"
    samtools sort -@ "${THREADS}" -o "${bam_possorted}" "${bam_fixmate}"
    samtools markdup -@ "${THREADS}" -s "${bam_possorted}" "${bam_markdup}" \
        2>> "${LOG_DIR}/stage6_markdup_stats.log"
    samtools index "${bam_markdup}"

    rm -f "${bam_namesorted}" "${bam_fixmate}" "${bam_possorted}"

done < "${SAMPLE_LIST}"

deactivate_env
log "Stage 6-1 완료"

# ---- 6-2) dupRadar (dupradar_env) ----
activate_env "${ENV_DUPRADAR}"
log "Stage 6-2: dupRadar 분석 시작"

while read -r sample; do
    [[ -z "${sample}" ]] && continue
    bam_markdup="${DIR_DUP}/${sample}.markdup.bam"
    if [[ ! -f "${bam_markdup}" ]]; then
        log "  [SKIP] ${sample}: markdup bam 없음"
        continue
    fi

    log "  [dupRadar] ${sample}"
    Rscript "${SCRIPT_DIR}/run_dupradar.R" \
        "${bam_markdup}" \
        "${REF_GTF}" \
        "${DIR_DUP}/${sample}_dupradar" \
        "${DUPRADAR_STRANDED}" \
        >> "${LOG_DIR}/stage6_dupradar.log" 2>&1

done < "${SAMPLE_LIST}"

deactivate_env
log "Stage 6 완료"

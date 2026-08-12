#!/usr/bin/env bash
# =========================================================
# Stage 3: Contamination screening (FastQ Screen)
# 사용 env: qc_tools_env (fastq-screen, bowtie2)
# 이미 결과 파일이 있는 샘플은 자동으로 SKIP (재실행 시 시간 절약)
# =========================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env "${ENV_QC_TOOLS}"

log "Stage 3: FastQ Screen 시작"

while read -r sample; do
    [[ -z "${sample}" ]] && continue

    r1="${DIR_TRIM}/${sample}_R1.trim.fastq.gz"

    # 이미 결과가 있으면 skip
    existing=$(find "${DIR_SCREEN}" -maxdepth 1 -iname "${sample}_R1.trim_screen.txt" 2>/dev/null | head -1)
    if [[ -n "${existing}" ]]; then
        log "  [ALREADY DONE] ${sample} -> skip"
        continue
    fi

    if [[ ! -f "${r1}" ]]; then
        log "  [SKIP] ${sample}: trimmed fastq 없음 (${r1})"
        continue
    fi

    log "  [FastQ Screen] ${sample}"
    fastq_screen \
        --conf "${FASTQ_SCREEN_CONF}" \
        --aligner bowtie2 \
        --threads "${THREADS}" \
        --outdir "${DIR_SCREEN}" \
        "${r1}" \
        >> "${LOG_DIR}/stage3_fastqscreen.log" 2>&1

done < "${SAMPLE_LIST}"

log "Stage 3: MultiQC 통합"
multiqc "${DIR_SCREEN}" -o "${DIR_SCREEN}" -n multiqc_screen --force \
    >> "${LOG_DIR}/stage3_multiqc.log" 2>&1

deactivate_env

log "Stage 3 완료"
log "  -> DIR_SCREEN 결과 확인 후 오염 의심 샘플 있으면 다음 단계 진행 전 판단 필요"

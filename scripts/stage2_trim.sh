#!/usr/bin/env bash
# =========================================================
# Stage 2: Adapter/Quality trimming (fastp)
# 사용 env: 시스템 PATH (fastp) - conda 불필요
# 이미 결과 파일이 있는 샘플은 자동으로 SKIP (재실행 시 시간 절약)
# =========================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

log "Stage 2: fastp trimming 시작"

while read -r sample; do
    [[ -z "${sample}" ]] && continue

    out1="${DIR_TRIM}/${sample}_R1.trim.fastq.gz"
    out2="${DIR_TRIM}/${sample}_R2.trim.fastq.gz"

    # 이미 결과가 있으면 skip
    if [[ -f "${out1}" && -f "${out2}" ]]; then
        log "  [ALREADY DONE] ${sample} -> skip"
        continue
    fi

    r1=$(find_fastq "${sample}" 1) || { log "  [SKIP] ${sample}: R1 fastq 못 찾음"; continue; }
    r2=$(find_fastq "${sample}" 2) || { log "  [SKIP] ${sample}: R2 fastq 못 찾음"; continue; }

    log "  [fastp] ${sample}"
    fastp \
        -i "${r1}" -I "${r2}" \
        -o "${out1}" -O "${out2}" \
        --thread "${THREADS}" \
        --detect_adapter_for_pe \
        --json "${DIR_TRIM}/${sample}.fastp.json" \
        --html "${DIR_TRIM}/${sample}.fastp.html" \
        >> "${LOG_DIR}/stage2_fastp.log" 2>&1

done < "${SAMPLE_LIST}"

log "Stage 2: MultiQC 통합 (fastp report 포함)"
multiqc "${DIR_TRIM}" -o "${DIR_TRIM}" -n multiqc_trim --force \
    >> "${LOG_DIR}/stage2_multiqc.log" 2>&1

log "Stage 2 완료"

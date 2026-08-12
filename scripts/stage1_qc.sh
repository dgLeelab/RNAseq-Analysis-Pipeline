#!/usr/bin/env bash
# =========================================================
# Stage 1: Raw read QC (FastQC + MultiQC)
# 사용 env: 시스템 PATH (fastqc, multiqc) - conda 불필요
# 이미 결과 파일이 있는 샘플은 자동으로 SKIP (재실행 시 시간 절약)
#
# 주의: FastQC 결과 파일명은 sample_id가 아니라 실제 입력 fastq
#       파일명(raw_prefix 기반)을 그대로 사용하므로, skip 체크는
#       find_fastq()로 알아낸 실제 파일명 기준으로 해야 함.
# =========================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

log "Stage 1: FastQC (raw) 시작"

while read -r sample; do
    [[ -z "${sample}" ]] && continue

    r1=$(find_fastq "${sample}" 1) || { log "  [SKIP] ${sample}: R1 fastq 못 찾음"; continue; }
    r2=$(find_fastq "${sample}" 2) || { log "  [SKIP] ${sample}: R2 fastq 못 찾음"; continue; }

    # FastQC 출력 파일명은 입력 파일 basename 기준 (예: 01-1_S41_L002_R1_001_fastqc.zip)
    r1_base=$(basename "${r1}" .fastq.gz)
    r2_base=$(basename "${r2}" .fastq.gz)

    if [[ -f "${DIR_QC_RAW}/${r1_base}_fastqc.zip" && -f "${DIR_QC_RAW}/${r2_base}_fastqc.zip" ]]; then
        log "  [ALREADY DONE] ${sample} -> skip"
        continue
    fi

    log "  [FastQC] ${sample} (${r1}, ${r2})"
    fastqc -t "${THREADS}" -o "${DIR_QC_RAW}" "${r1}" "${r2}" \
        >> "${LOG_DIR}/stage1_fastqc.log" 2>&1

done < "${SAMPLE_LIST}"

log "Stage 1: MultiQC 통합"
multiqc "${DIR_QC_RAW}" -o "${DIR_QC_RAW}" -n multiqc_raw --force \
    >> "${LOG_DIR}/stage1_multiqc.log" 2>&1

log "Stage 1 완료"

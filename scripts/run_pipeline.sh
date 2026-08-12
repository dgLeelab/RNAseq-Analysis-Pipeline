#!/usr/bin/env bash
# =========================================================
# run_pipeline.sh
# Stage 1~6 순차 실행 (QC -> Trim -> Screen -> Align -> RSeQC/Strand -> dupRadar)
# 정량(StringTie 이후)은 ../modules/ 이하 스크립트로 이어서 진행합니다.
#
# 사용: ./run_pipeline.sh [start_stage] [end_stage]
#   예) ./run_pipeline.sh 1 6   (전체)
#       ./run_pipeline.sh 4 5   (align + strand 판별만 먼저 테스트)
# =========================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

START="${1:-1}"
END="${2:-6}"

STAGES=(
    "1:${SCRIPT_DIR}/stage1_qc.sh:FastQC/MultiQC (raw QC)"
    "2:${SCRIPT_DIR}/stage2_trim.sh:fastp (trimming + dup rate 1차)"
    "3:${SCRIPT_DIR}/stage3_screen.sh:FastQ Screen (contamination)"
    "4:${SCRIPT_DIR}/stage4_align.sh:HISAT2 (alignment)"
    "5:${SCRIPT_DIR}/stage5_rseqc.sh:RSeQC + strand 자동판별"
    "6:${SCRIPT_DIR}/stage6_dupradar.sh:samtools markdup + dupRadar"
)

log "=========================================="
log "RNA-seq 파이프라인 시작 (stage ${START} ~ ${END})"
log "=========================================="

for entry in "${STAGES[@]}"; do
    stage_num="${entry%%:*}"
    rest="${entry#*:}"
    script_path="${rest%%:*}"
    desc="${rest#*:}"

    if (( stage_num < START || stage_num > END )); then
        continue
    fi

    log "------------------------------------------"
    log ">>> Stage ${stage_num}: ${desc}"
    log "------------------------------------------"

    bash "${script_path}"

    if [[ $? -ne 0 ]]; then
        log "[ERROR] Stage ${stage_num} 실패. 파이프라인 중단."
        exit 1
    fi
done

log "=========================================="
log "파이프라인 종료 (stage ${START} ~ ${END})"
log "=========================================="

log "최종 통합 MultiQC 리포트 생성"
multiqc "${OUT_ROOT}" -o "${DIR_MULTIQC}" -n multiqc_final \
    --ignore "${DIR_MULTIQC}/*" \
    >> "${LOG_DIR}/final_multiqc.log" 2>&1

log "완료. 통합 리포트: ${DIR_MULTIQC}/multiqc_final.html"
log "다음 단계: ../modules/ 에서 정량(StringTie) 이후 진행"

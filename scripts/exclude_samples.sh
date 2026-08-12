#!/usr/bin/env bash
# =========================================================
# exclude_samples.sh
# QC 결과 검토 후 특정 샘플을 samples.txt에서 제외.
# 제외된 샘플은 samples_excluded.log에 사유와 함께 기록되고,
# 원본 samples.txt는 samples.txt.bak_{timestamp}로 백업됨.
#
# 사용: ./exclude_samples.sh <sample_id> "<제외 사유>"
# 예:   ./exclude_samples.sh CMJ076_1 "RNA degradation, rRNA screen 이상"
#
# 여러 개 한 번에 빼려면 반복 호출하면 됨.
# =========================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SAMPLE_ID="${1:?사용법: ./exclude_samples.sh <sample_id> \"<사유>\"}"
REASON="${2:-사유 미기재}"

if ! grep -qx "${SAMPLE_ID}" "${SAMPLE_LIST}"; then
    log "[WARNING] ${SAMPLE_ID} 가 현재 samples.txt에 없습니다 (이미 제외됐거나 오타 가능성)"
    exit 1
fi

# 백업
cp "${SAMPLE_LIST}" "${SAMPLE_LIST}.bak_$(date '+%Y%m%d_%H%M%S')"

# 제외 처리
grep -vx "${SAMPLE_ID}" "${SAMPLE_LIST}" > "${SAMPLE_LIST}.tmp"
mv "${SAMPLE_LIST}.tmp" "${SAMPLE_LIST}"

# 로그 기록
EXCLUDE_LOG="${PROJECT_DIR}/samples_excluded.log"
echo -e "$(date '+%Y-%m-%d %H:%M:%S')\t${SAMPLE_ID}\t${REASON}" >> "${EXCLUDE_LOG}"

log "[제외 완료] ${SAMPLE_ID} (${REASON})"
log "  -> 남은 샘플 수: $(wc -l < "${SAMPLE_LIST}")"
log "  -> 제외 이력: ${EXCLUDE_LOG}"

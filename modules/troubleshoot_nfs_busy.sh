#!/usr/bin/env bash
# =========================================================
# troubleshoot_nfs_busy.sh
# NFS 파일시스템에서 "rm: cannot remove ... .nfsXXXX: Device or resource busy"
# 에러가 날 때, 어떤 프로세스가 해당 파일을 붙잡고 있는지 확인/정리합니다.
#
# 사용: ./troubleshoot_nfs_busy.sh <문제되는 디렉토리>
# =========================================================
set -uo pipefail

TARGET_DIR="${1:?사용법: $0 <디렉토리>}"

echo "[troubleshoot] ${TARGET_DIR} 를 붙잡고 있는 프로세스 확인:"
lsof +D "${TARGET_DIR}" 2>/dev/null

echo ""
echo "[troubleshoot] R/Trinity 관련 프로세스 확인:"
ps aux | grep -E "R |Rscript|analyze_diff|PtR|define_clusters" | grep -v grep

echo ""
echo "위 목록에서 관련 PID를 찾아 종료하세요:"
echo "  kill <PID>      # 정상 종료 시도"
echo "  kill -9 <PID>   # 안 죽으면 강제 종료"
echo ""
echo "주의: 다른 터미널 세션에서 해당 디렉토리에 cd 되어 있는 것만으로도"
echo "      NFS busy가 발생할 수 있습니다. 모든 세션에서 해당 디렉토리를"
echo "      벗어난 뒤 다시 시도하세요."

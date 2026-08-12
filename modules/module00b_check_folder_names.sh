#!/usr/bin/env bash
# =========================================================
# module00b_check_folder_names.sh
# StringTie 결과 폴더명은 이후 matrix 컬럼명 -> samples.file 그룹 라벨로
# 그대로 이어지므로, 오타가 있으면 여기서 미리 잡아야 합니다.
#
# 기대 샘플 수는 하드코딩하지 않고, config.sh의 samples.txt에서 자동으로
# 셉니다 (config.sh를 못 찾으면 기대값 비교 없이 현황만 보여줍니다).
# 인자로 직접 기대값을 주면 그 값이 우선합니다.
#
# 사용: ./module00b_check_folder_names.sh [기대 샘플 수(선택), 생략 시 samples.txt에서 자동 계산]
# 실행 위치: StringTie 결과 상위 디렉토리 (각 {sample}/ 폴더가 있는 곳)
#            예: results/07_stringtie/ 에서 ../../../modules/module00b_...sh 로 실행
# =========================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_CANDIDATE="${SCRIPT_DIR}/../scripts/config.sh"

EXPECTED_COUNT="${1:-}"

if [[ -z "${EXPECTED_COUNT}" ]]; then
    if [[ -f "${CONFIG_CANDIDATE}" ]]; then
        # shellcheck disable=SC1090
        SAMPLE_LIST_PATH=$(grep -oP '^SAMPLE_LIST="\K[^"]+' "${CONFIG_CANDIDATE}" 2>/dev/null || true)
        # config.sh 안의 변수 치환(${PROJECT_DIR} 등)까지 반영해서 실제 경로 얻기
        if [[ -n "${SAMPLE_LIST_PATH}" ]]; then
            RESOLVED_SAMPLE_LIST=$(bash -c "source '${CONFIG_CANDIDATE}' &>/dev/null; echo \"\${SAMPLE_LIST}\"")
            if [[ -f "${RESOLVED_SAMPLE_LIST}" ]]; then
                EXPECTED_COUNT=$(grep -c . "${RESOLVED_SAMPLE_LIST}")
                echo "[module00b] config.sh의 samples.txt(${RESOLVED_SAMPLE_LIST})에서 기대 샘플 수 자동 계산: ${EXPECTED_COUNT}"
            fi
        fi
    fi
fi

if [[ -z "${EXPECTED_COUNT}" ]]; then
    echo "[module00b] WARNING: 기대 샘플 수를 알 수 없습니다 (config.sh/samples.txt를 못 찾음)."
    echo "  현재 폴더 개수만 보여드리고, 기대값 비교는 생략합니다."
    echo "  기대값을 직접 지정하려면: $0 <기대_샘플수>"
fi

echo ""
echo "[module00b] 현재 샘플 폴더 목록:"
ls -d */ 2>/dev/null | sed 's#/$##' | sort

actual_count=$(ls -d */ 2>/dev/null | wc -l)
echo ""
if [[ -n "${EXPECTED_COUNT}" ]]; then
    echo "[module00b] 폴더 개수: ${actual_count} (기대값: ${EXPECTED_COUNT})"
    if [[ "${actual_count}" -ne "${EXPECTED_COUNT}" ]]; then
        echo "[module00b] WARNING: 개수가 기대값과 다릅니다. 누락/중복 폴더 확인 필요."
    fi
else
    echo "[module00b] 폴더 개수: ${actual_count} (기대값 없음, 비교 생략)"
fi

echo ""
echo "[module00b] 그룹(품종 등)별 replicate 개수 ({group}_{숫자} 형식 가정):"
ls -d */ 2>/dev/null | sed 's#/$##' | sed -E 's/_[0-9]+$//' | sort | uniq -c

echo ""
echo "[module00b] 오타 의심 패턴 체크 (마지막 접미사가 양의 정수가 아닌 경우):"
found_bad=0
ls -d */ 2>/dev/null | sed 's#/$##' | while read -r name; do
    suffix="${name##*_}"
    if ! [[ "${suffix}" =~ ^[0-9]+$ ]]; then
        echo "  [의심] ${name} (replicate 번호가 숫자가 아님 -> {group}_{숫자} 형식이 아닐 수 있음)"
        found_bad=1
    fi
done

echo ""
echo "오타 발견 시 예시: mv CMJ03_3 CMJ032_3"

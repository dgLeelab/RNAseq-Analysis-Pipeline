#!/usr/bin/env bash
# module05 실행 wrapper
# module04_stringtie.sh와 같은 디렉토리(각 sample 하위폴더가 있는 곳)에서 실행
# 출력: gene_count_matrix.csv, transcript_count_matrix.csv (prepDE.py 기본 파일명)
# 사용: ./run_module05.sh
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[module05] prepDE.py (python3) 실행 -> gene/transcript count matrix"
echo "  주의: StringTie가 -e 옵션으로 실행된 GTF여야 coverage가 인식됩니다."
python3 "${SCRIPT_DIR}/module05_prepDE_after_Stringtie.py3" -i .

if [[ ! -f gene_count_matrix.csv ]]; then
    echo "[module05] ERROR: gene_count_matrix.csv 생성 실패."
    echo "  'Error: no GTF files found' 가 떴다면 -i 경로 문제이거나,"
    echo "  각 샘플 폴더의 GTF가 -e 옵션 없이 생성됐을 수 있습니다."
    exit 1
fi

echo "[module05] 완료: gene_count_matrix.csv, transcript_count_matrix.csv 생성됨"

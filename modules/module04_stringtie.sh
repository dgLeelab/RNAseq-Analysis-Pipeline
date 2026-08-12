#!/usr/bin/env bash
# =========================================================
# module04_stringtie.sh
# HISAT2 정렬 BAM -> 샘플별 StringTie 정량 (-e -B)
#
# 주의: -G에는 반드시 module00_prep_annotation.sh로 전처리된 annotation을
#       사용해야 합니다 (organelle trans-splicing 유전자, 빈 transcript_id,
#       gene_id 내 세미콜론 문제 해결된 버전). HISAT2 정렬 자체는 organelle
#       포함 genome으로 했어도 무방합니다 (annotation에서만 제외하면 됨).
#
# strand 옵션은 Stage 5(RSeQC)가 생성한 strand_info.sh를 자동으로 읽어
# --fr/--rf/미지정을 자동 적용합니다 (하드코딩 금지).
#
# 실행 위치: results/ 디렉토리 기준 상대경로 사용 (04_hisat2와 같은 레벨)
# =========================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_CANDIDATE="${SCRIPT_DIR}/../scripts/config.sh"

if [[ -f "${CONFIG_CANDIDATE}" ]]; then
    # shellcheck disable=SC1090
    source "${CONFIG_CANDIDATE}"
    BAM_DIR="${DIR_ALIGN}"
    STRAND_FILE="${STRAND_INFO_FILE}"
else
    BAM_DIR="../04_hisat2"
    STRAND_FILE="../strand_info.sh"
fi

# -G에 쓸 annotation (module00 결과) - 반드시 명시적으로 지정
GFF="${1:?사용법: $0 <module00 전처리된 nuclear annotation .gff/.gtf>}"

thread_per_sample="${THREADS_PER_SAMPLE:-24}"
parallel_samples="${PARALLEL_SAMPLES:-8}"

STRINGTIE_STRAND_OPT=""
if [[ -f "${STRAND_FILE}" ]]; then
    # shellcheck disable=SC1090
    source "${STRAND_FILE}"
    echo "[module04] strand 자동판별 결과 사용: ${STRANDEDNESS} (${STRINGTIE_STRAND_OPT:-unstranded})"
else
    echo "[module04] WARNING: strand_info.sh 없음 -> unstranded로 진행 (Stage 5 먼저 실행 권장)"
fi

date
echo "StringTie start! (동시 ${parallel_samples}개 샘플, 샘플당 ${thread_per_sample}스레드, GFF=${GFF})"

run_stringtie_one () {
    local bam="$1"
    local sample
    sample=$(basename "$bam" .sorted.bam)

    mkdir -p "$sample"

    if [[ -f "$sample/transcripts.gtf" ]]; then
        echo "[ALREADY DONE] $sample -> skip"
        return 0
    fi

    echo "[StringTie 시작] $sample"
    stringtie -eB \
        -G "$GFF" \
        -A "$sample/gene_abundances.tsv" \
        -o "$sample/transcripts.gtf" \
        -p "$thread_per_sample" \
        ${STRINGTIE_STRAND_OPT} \
        "$bam"
    echo "[완료] $sample"
}
export -f run_stringtie_one
export GFF thread_per_sample STRINGTIE_STRAND_OPT

job_count=0
for bam in "${BAM_DIR}"/*.sorted.bam
do
    run_stringtie_one "$bam" &
    job_count=$((job_count + 1))

    if (( job_count >= parallel_samples )); then
        wait -n
        job_count=$((job_count - 1))
    fi
done

wait

date
echo "done!"

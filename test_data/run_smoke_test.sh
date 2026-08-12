#!/usr/bin/env bash
# =========================================================
# run_smoke_test.sh
# 합성 테스트 데이터로 전체 파이프라인(Stage 1~6 + modules 일부)을
# 자동 실행해서, clone 직후 환경이 제대로 세팅됐는지 수 분 내에 검증합니다.
#
# 사전 조건: setup.sh 실행 완료, conda env(qc_tools_env/dupradar_env/
#            trinity_env) 및 fastqc/fastp/hisat2/stringtie/samtools 설치됨
#
# 사용: ./run_smoke_test.sh
# =========================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="${REPO_DIR}/test_data/_smoke_run"

echo "=========================================="
echo " Smoke Test: 합성 데이터로 전체 파이프라인 검증"
echo "=========================================="

rm -rf "${TEST_DIR}"
mkdir -p "${TEST_DIR}"

echo "[1/6] 테스트 데이터 생성"
python3 "${REPO_DIR}/test_data/generate_test_data.py" "${TEST_DIR}"

echo "[2/6] 테스트용 config.sh 생성"
CONDA_SH_DETECTED=$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh

cat > "${TEST_DIR}/config.sh" << CONFEOF
#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${TEST_DIR}"
RAW_DIR="\${PROJECT_DIR}/00_raw"
SAMPLE_LIST="\${PROJECT_DIR}/samples.txt"
SAMPLE_MAP="\${PROJECT_DIR}/sample_map.tsv"

REF_DIR="\${PROJECT_DIR}/ref"
REF_GENOME_FA="\${REF_DIR}/genome.fasta"
REF_GTF="\${REF_DIR}/genome.gtf"
HISAT2_INDEX_PREFIX="\${REF_DIR}/hisat2_index/genome"
BOWTIE2_TARGET_INDEX="\${REF_DIR}/bowtie2_index/genome"
RSEQC_BED="\${REF_DIR}/genome.bed"
FASTQ_SCREEN_CONF="\${REF_DIR}/fastq_screen.conf"
SCREEN_REFS_DIR="\${REF_DIR}/screen_refs"

OUT_ROOT="\${PROJECT_DIR}/results"
DIR_QC_RAW="\${OUT_ROOT}/01_fastqc_raw"
DIR_TRIM="\${OUT_ROOT}/02_fastp"
DIR_SCREEN="\${OUT_ROOT}/03_fastq_screen"
DIR_ALIGN="\${OUT_ROOT}/04_hisat2"
DIR_RSEQC="\${OUT_ROOT}/05_rseqc"
DIR_DUP="\${OUT_ROOT}/06_dupradar"
DIR_MULTIQC="\${OUT_ROOT}/multiqc"
LOG_DIR="\${OUT_ROOT}/logs"
STRAND_INFO_FILE="\${OUT_ROOT}/strand_info.sh"

THREADS=2
THREADS_PER_SAMPLE=2
PARALLEL_SAMPLES=2
PARALLEL_SAMPLES_LIGHT=4

ENV_DUPRADAR="dupradar_env"
ENV_QC_TOOLS="qc_tools_env"
ENV_TRINITY="trinity_env"

CONDA_SH="${CONDA_SH_DETECTED}"

activate_env () {
    local env_name="\$1"
    set +u
    source "\${CONDA_SH}"
    conda activate "\${env_name}"
    set -u
}

deactivate_env () {
    set +u
    conda deactivate
    set -u
}

find_fastq () {
    local sample="\$1"
    local read_num="\$2"
    local prefix="\${sample}"
    if [[ -f "\${SAMPLE_MAP}" ]]; then
        local mapped
        mapped=\$(awk -F'\t' -v s="\${sample}" '\$1==s {print \$2; exit}' "\${SAMPLE_MAP}")
        [[ -n "\${mapped}" ]] && prefix="\${mapped}"
    fi
    local candidates=(
        "\${RAW_DIR}/\${prefix}_R\${read_num}.fastq.gz"
        "\${RAW_DIR}/\${prefix}_S"*"_L"*"_R\${read_num}_"*".fastq.gz"
        "\${RAW_DIR}/\${prefix}"*"_R\${read_num}"*".fastq.gz"
    )
    for pattern in "\${candidates[@]}"; do
        for f in \${pattern}; do
            [[ -f "\${f}" ]] && { echo "\${f}"; return 0; }
        done
    done
    return 1
}

mkdir -p "\${DIR_QC_RAW}" "\${DIR_TRIM}" "\${DIR_SCREEN}" "\${DIR_ALIGN}" \\
         "\${DIR_RSEQC}" "\${DIR_DUP}" "\${DIR_MULTIQC}" "\${LOG_DIR}"

log () {
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$*"
}
CONFEOF

echo "[3/6] HISAT2 index 구축 (작은 genome이라 수 초 내 완료)"
mkdir -p "${TEST_DIR}/ref/hisat2_index"
hisat2-build -p 2 "${TEST_DIR}/ref/genome.fasta" "${TEST_DIR}/ref/hisat2_index/genome" > /dev/null

echo "[4/6] RSeQC용 genome.bed 생성"
if command -v gff3ToGenePred &> /dev/null; then
    gff3ToGenePred -useName "${TEST_DIR}/ref/genome.gff3" "${TEST_DIR}/ref/genome.genePred"
    genePredToBed "${TEST_DIR}/ref/genome.genePred" "${TEST_DIR}/ref/genome.bed"
else
    echo "  WARNING: gff3ToGenePred 없음 -> Stage 5는 스킵될 수 있습니다"
    echo "  (conda activate qc_tools_env 후 ucsc-gff3togenepred, ucsc-genepredtobed 설치 필요)"
    touch "${TEST_DIR}/ref/genome.bed"
fi

echo "[5/6] FastQ Screen 참조 생략 (테스트용 최소 데이터라 스킵, Stage 3는 target genome만 대상)"
mkdir -p "${TEST_DIR}/ref/screen_refs"/{rRNA,ecoli,phix,human}
# 최소한의 fastq_screen.conf만 생성 (target genome만 포함, 스모크 테스트 목적)
if command -v bowtie2-build &> /dev/null; then
    mkdir -p "${TEST_DIR}/ref/bowtie2_index"
    bowtie2-build "${TEST_DIR}/ref/genome.fasta" "${TEST_DIR}/ref/bowtie2_index/genome" > /dev/null 2>&1 || true
fi
cat > "${TEST_DIR}/ref/fastq_screen.conf" << SCREENEOF
BOWTIE2 $(which bowtie2 2>/dev/null || echo bowtie2)
THREADS 2
DATABASE  Target  ${TEST_DIR}/ref/bowtie2_index/genome
SCREENEOF

echo "[6/6] Stage 1~6 실행"
cd "${REPO_DIR}/scripts"
CONFIG_BACKUP=""
if [[ -f config.sh ]]; then
    CONFIG_BACKUP="config.sh.smoketest_backup"
    mv config.sh "${CONFIG_BACKUP}"
fi
cp "${TEST_DIR}/config.sh" config.sh

set +e
bash run_pipeline.sh 1 6
RESULT=$?
set -e

rm -f config.sh
if [[ -n "${CONFIG_BACKUP}" ]]; then
    mv "${CONFIG_BACKUP}" config.sh
fi

echo ""
echo "=========================================="
if [[ "${RESULT}" -eq 0 ]]; then
    echo " Smoke Test 성공: Stage 1~6 전부 정상 실행됨"
else
    echo " Smoke Test 실패: 위 로그와 ${TEST_DIR}/results/logs/ 를 확인하세요"
fi
echo "=========================================="
echo "테스트 결과 위치: ${TEST_DIR}/results/"
echo "정리하려면: rm -rf ${TEST_DIR}"
exit "${RESULT}"

#!/usr/bin/env bash
# =========================================================
# setup_screen_references.sh
# FastQ Screen용 오염원 reference 다운로드 + bowtie2 index 구축
#   - rRNA (SILVA), E.coli, PhiX, Human  (범용, species 안 바뀜 -> 재사용 가능)
#   - Target genome (프로젝트마다 다름, config.sh의 REF_GENOME_FA 기준으로 자동 구축)
#
# config.sh를 source해서 모든 경로를 자동으로 가져옵니다 (하드코딩 없음).
# 실행 전 반드시 config.sh의 REF_DIR, REF_GENOME_FA 등을 프로젝트에 맞게
# 설정해두세요.
#
# 사용: ./setup_screen_references.sh
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

THREADS_BUILD="${THREADS:-16}"

mkdir -p "${SCREEN_REFS_DIR}"/{rRNA,ecoli,phix,human}
mkdir -p "$(dirname "${BOWTIE2_TARGET_INDEX}")"

activate_env "${ENV_QC_TOOLS}"

# ---------------------------------------------------------
# 0) Target genome bowtie2 index (config.sh의 REF_GENOME_FA 기준)
#    주의: HISAT2 index(.ht2)와 bowtie2 index(.bt2)는 포맷이 달라 호환 안 됨.
#          FastQ Screen용으로 별도 구축 필요.
# ---------------------------------------------------------
if [[ -f "${BOWTIE2_TARGET_INDEX}.1.bt2" || -f "${BOWTIE2_TARGET_INDEX}.1.bt2l" ]]; then
    echo "[0/4] Target genome bowtie2 index 이미 존재 -> skip"
else
    echo "[0/4] Target genome bowtie2 index 구축: ${REF_GENOME_FA}"
    bowtie2-build --threads "${THREADS_BUILD}" "${REF_GENOME_FA}" "${BOWTIE2_TARGET_INDEX}"
fi

# ---------------------------------------------------------
# 1) rRNA (SILVA SSU + LSU, Non-redundant)
# ---------------------------------------------------------
cd "${SCREEN_REFS_DIR}/rRNA"
if [[ -f "rRNA_index.1.bt2" || -f "rRNA_index.1.bt2l" ]]; then
    echo "[1/4] rRNA index 이미 존재 -> skip"
else
    echo "[1/4] rRNA (SILVA) 다운로드"
    # SILVA release 번호는 주기적으로 바뀌므로 https://www.arb-silva.de/download/arb-files/ 에서 최신 버전 확인 권장
    wget -c "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
    wget -c "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"
    gunzip -k SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
    gunzip -k SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz
    cat SILVA_138.2_SSURef_NR99_tax_silva.fasta SILVA_138.2_LSURef_NR99_tax_silva.fasta > rRNA_combined.fasta
    bowtie2-build --threads "${THREADS_BUILD}" rRNA_combined.fasta rRNA_index
fi

# ---------------------------------------------------------
# 2) E. coli (NCBI RefSeq K-12 MG1655)
# ---------------------------------------------------------
cd "${SCREEN_REFS_DIR}/ecoli"
if [[ -f "ecoli_index.1.bt2" || -f "ecoli_index.1.bt2l" ]]; then
    echo "[2/4] E. coli index 이미 존재 -> skip"
else
    echo "[2/4] E. coli 다운로드"
    wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    gunzip -k GCF_000005845.2_ASM584v2_genomic.fna.gz
    bowtie2-build --threads "${THREADS_BUILD}" GCF_000005845.2_ASM584v2_genomic.fna ecoli_index
fi

# ---------------------------------------------------------
# 3) PhiX (Illumina control)
# ---------------------------------------------------------
cd "${SCREEN_REFS_DIR}/phix"
if [[ -f "phix_index.1.bt2" || -f "phix_index.1.bt2l" ]]; then
    echo "[3/4] PhiX index 이미 존재 -> skip"
else
    echo "[3/4] PhiX 다운로드"
    wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
    gunzip -k GCF_000819615.1_ViralProj14015_genomic.fna.gz
    bowtie2-build --threads "${THREADS_BUILD}" GCF_000819615.1_ViralProj14015_genomic.fna phix_index
fi

# ---------------------------------------------------------
# 4) Human (GRCh38, primary assembly) - 용량 매우 큼 (fasta ~3GB)
# ---------------------------------------------------------
cd "${SCREEN_REFS_DIR}/human"
if [[ -f "human_index.1.bt2" || -f "human_index.1.bt2l" ]]; then
    echo "[4/4] Human index 이미 존재 -> skip"
else
    echo "[4/4] Human GRCh38 다운로드 (용량 큼, 시간 소요)"
    wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    gunzip -k GCF_000001405.40_GRCh38.p14_genomic.fna.gz
    bowtie2-build --threads "${THREADS_BUILD}" GCF_000001405.40_GRCh38.p14_genomic.fna human_index
fi

deactivate_env

# ---------------------------------------------------------
# fastq_screen.conf 자동 생성 (config.sh 경로 기준)
# ---------------------------------------------------------
BOWTIE2_BIN=$(source "${CONDA_SH}" && conda activate "${ENV_QC_TOOLS}" && which bowtie2)

cat > "${FASTQ_SCREEN_CONF}" << CONFEOF
BOWTIE2 ${BOWTIE2_BIN}
THREADS ${THREADS_BUILD}

DATABASE  Target        ${BOWTIE2_TARGET_INDEX}
DATABASE  rRNA          ${SCREEN_REFS_DIR}/rRNA/rRNA_index
DATABASE  EColi         ${SCREEN_REFS_DIR}/ecoli/ecoli_index
DATABASE  PhiX          ${SCREEN_REFS_DIR}/phix/phix_index
DATABASE  Human         ${SCREEN_REFS_DIR}/human/human_index
CONFEOF

echo ""
echo "=========================================="
echo "모든 reference index 구축 완료"
echo "fastq_screen.conf 생성됨: ${FASTQ_SCREEN_CONF}"
echo "=========================================="
cat "${FASTQ_SCREEN_CONF}"

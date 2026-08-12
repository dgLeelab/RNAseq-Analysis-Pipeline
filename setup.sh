#!/usr/bin/env bash
# =========================================================
# setup.sh
# git clone 직후 딱 한 번 실행하는 초기 설정 스크립트.
#
# 하는 일:
#   1) scripts/config.sh 가 없으면 config.sh.example 을 복사해서 생성
#   2) 필수 프로그램(conda, hisat2, stringtie, samtools 등) 존재 여부 확인
#   3) config.sh를 대화형으로 채워넣을지, 나중에 직접 편집할지 안내
#
# 사용: ./setup.sh
# =========================================================
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${REPO_DIR}/scripts/config.sh"
CONFIG_EXAMPLE="${REPO_DIR}/scripts/config.sh.example"

echo "=========================================="
echo " RNA-seq Pipeline 초기 설정"
echo "=========================================="
echo ""

# ---- 1) config.sh 생성 ----
if [[ -f "${CONFIG_FILE}" ]]; then
    echo "[1/3] scripts/config.sh 이미 존재 -> 그대로 둡니다."
else
    cp "${CONFIG_EXAMPLE}" "${CONFIG_FILE}"
    echo "[1/3] scripts/config.sh.example -> scripts/config.sh 로 복사했습니다."

    # CONDA_SH는 자동 감지 가능하므로 여기서 바로 채워줌 (sed 치환 실수 방지)
    if command -v conda &> /dev/null; then
        CONDA_BASE=$(conda info --base 2>/dev/null || true)
        if [[ -n "${CONDA_BASE}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
            DETECTED_CONDA_SH="${CONDA_BASE}/etc/profile.d/conda.sh"
            # macOS/Linux sed 차이 없이 안전하게 python으로 치환
            python3 - "${CONFIG_FILE}" "${DETECTED_CONDA_SH}" << 'PYEOF'
import sys, re
config_file, conda_sh = sys.argv[1], sys.argv[2]
with open(config_file) as f:
    content = f.read()
content = re.sub(r'^CONDA_SH=.*$', f'CONDA_SH="{conda_sh}"', content, count=1, flags=re.MULTILINE)
with open(config_file, "w") as f:
    f.write(content)
PYEOF
            echo "  -> CONDA_SH 자동 감지 및 반영: ${DETECTED_CONDA_SH}"
        else
            echo "  -> CONDA_SH 자동 감지 실패. config.sh에서 직접 확인/수정하세요 (conda info --base 참고)."
        fi
    fi
fi

# ---- 2) 필수 프로그램 확인 ----
echo ""
echo "[2/3] 필수 프로그램 확인:"
check_cmd () {
    local cmd="$1"
    if command -v "${cmd}" &> /dev/null; then
        echo "  OK   : ${cmd} ($(command -v "${cmd}"))"
    else
        echo "  MISSING : ${cmd}  (PATH에 없음, conda env 활성화 또는 설치 필요)"
    fi
}
check_cmd conda
check_cmd fastqc
check_cmd fastp
check_cmd hisat2
check_cmd samtools
check_cmd stringtie
check_cmd multiqc
echo "  (fastq_screen, bowtie2, rseqc, gff3ToGenePred 는 qc_tools_env conda env에"
echo "   있는지 확인하세요: conda activate qc_tools_env 후 다시 이 스크립트를 실행하면"
echo "   추가로 체크됩니다.)"

if command -v conda &> /dev/null; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    if conda env list 2>/dev/null | grep -q "qc_tools_env"; then
        conda activate qc_tools_env 2>/dev/null || true
        check_cmd fastq_screen
        check_cmd bowtie2
        check_cmd read_distribution.py
        check_cmd gff3ToGenePred
        conda deactivate 2>/dev/null || true
    else
        echo "  MISSING : qc_tools_env conda env 자체가 없습니다."
        echo "    생성 명령: conda create -n qc_tools_env -c conda-forge -c bioconda fastq-screen bowtie2 rseqc ucsc-gff3togenepred ucsc-genepredtobed"
    fi
    if conda env list 2>/dev/null | grep -q "dupradar_env"; then
        echo "  OK   : dupradar_env 존재"
    else
        echo "  MISSING : dupradar_env conda env가 없습니다. (Stage 6 dupRadar용)"
    fi
    if conda env list 2>/dev/null | grep -q "trinity_env"; then
        echo "  OK   : trinity_env 존재"
    else
        echo "  MISSING : trinity_env conda env가 없습니다. (modules/module08, module09용)"
    fi
fi

# ---- 3) config.sh 수정 안내 ----
echo ""
echo "[3/3] 다음 단계:"
echo ""
echo "  scripts/config.sh 를 열어 아래 항목을 프로젝트에 맞게 수정하세요:"
echo "    - PROJECT_DIR          : 프로젝트 루트 (00_raw/, samples.txt 가 있는 곳)"
echo "    - REF_DIR               : reference 파일들 루트"
echo "    - REF_GENOME_FA, REF_GTF"
echo "    - HISAT2_INDEX_PREFIX   : hisat2-build 로 만든 index prefix"
echo "    - CONDA_SH              : conda 초기화 스크립트 경로 (보통 <conda_base>/etc/profile.d/conda.sh)"
echo "    - THREADS, THREADS_PER_SAMPLE, PARALLEL_SAMPLES : 서버 코어 수(nproc)에 맞게 조정"
echo ""
echo "  수정 후 확인:"
echo "    bash -n scripts/config.sh && echo OK"
echo ""
echo "  Reference index가 아직 없다면:"
echo "    cd scripts"
echo "    hisat2-build -p \$(nproc) \${REF_GENOME_FA} \${HISAT2_INDEX_PREFIX}   # (config.sh 값 참고해서 수동 실행)"
echo "    ./setup_screen_references.sh          # FastQ Screen용 reference 자동 구축"
echo "    (RSeQC용 genome.bed는 README의 'RSeQC BED 생성' 절 참고 -"
echo "     gff3ToGenePred + genePredToBed 사용 권장, organelle scaffold 사전 제외 필요)"
echo ""
echo "  samples.txt, sample_map.tsv 준비:"
echo "    ./list_raw_prefixes.sh   # raw fastq 파일명 확인 후 sample_map.tsv 뼈대 생성"
echo ""
echo "  모든 준비가 끝나면:"
echo "    cd scripts && nohup ./run_pipeline.sh > run_pipeline.log 2>&1 &"
echo ""
echo "설정 안내 끝."

# RNA-seq Analysis Pipeline (HISAT2 · StringTie · edgeR)

Bulk RNA-seq(paired-end short read) 데이터를 대상으로 한
QC → 트리밍 → 오염 스크리닝 → 정렬 → strand 판별 → 중복률 QC → 정량 →
차등발현/클러스터링까지 이어지는 파이프라인입니다. 특정 종/샘플에
종속되지 않도록 작성되어, reference genome/GTF와 `samples.txt`만
바꾸면 다른 프로젝트에도 재사용할 수 있습니다. bash 기반 stage/module
스크립트를 순서대로 실행하는 구조이며, Snakemake/Nextflow 같은
워크플로 매니저는 사용하지 않습니다.

> **현재 상태**: Stage 1~6(QC~정렬~strand 판별~dupRadar)과
> `modules/`(StringTie 정량~edgeR 차등발현~클러스터링)이 콩(Glycine max)
> 12품종 × 3 반복(36샘플) 프로젝트에서 처음부터 끝까지 실행 검증되었습니다.

## 파이프라인 구조

```
[사전준비, 1회성]
  scripts/list_raw_prefixes.sh          raw fastq 파일명 -> sample_map.tsv 뼈대 생성
  scripts/setup_screen_references.sh    FastQ Screen용 오염 참조(rRNA/E.coli/PhiX/Human) index 구축

[stage 1] FastQC + MultiQC              (raw read QC)
     |  00_raw/*.fastq.gz
     v
[stage 2] fastp                         (adapter/quality trimming)
     |  02_fastp/{sample}_R{1,2}.trim.fastq.gz
     v
[stage 3] FastQ Screen                  (오염원 스크리닝)
     v
[stage 4] HISAT2 + samtools sort/index  (정렬, 병렬 처리)
     |  04_hisat2/{sample}.sorted.bam
     v
[stage 5] RSeQC + strand 자동판별       (infer_experiment / read_distribution / geneBody_coverage)
     |  results/strand_info.sh 생성 (이후 단계가 자동으로 읽어감)
     v
[stage 6] samtools markdup + dupRadar   (중복률 QC, PCR bias vs 발현량 기반 dup 구분)
     v
[modules/] StringTie -> prepDE -> TPM -> LOC 매핑 -> TMM -> PtR/edgeR DE·클러스터링
```

각 stage/module 스크립트는 `config.sh`를 source(또는 상대경로로 참조)해서
공통 경로/자원 설정을 가져오며, 이미 처리된 샘플은 결과 파일 존재 여부로
자동 skip합니다(재실행 시 이어서 처리 가능). 샘플 수나 조건 구성에 대한
가정은 코드 안에 없으며, `samples.txt`에 나열한 sample_id 목록을 그대로
순회합니다.

## 사용 툴 / 버전

| 툴 | 버전 | 비고 |
|---|---|---|
| FastQC | v0.12.1 | |
| fastp | v1.0.1 | |
| FastQ Screen | 0.16.0 | conda env `qc_tools_env` |
| bowtie2 | 2.5.5 | FastQ Screen 및 오염 참조 index 빌드용 |
| HISAT2 | 2.2.1 | splice-aware(`_tran`) index 권장 |
| samtools | 1.21 | |
| RSeQC | 5.0.5 | conda env `qc_tools_env` |
| StringTie | 3.0.0 | |
| MultiQC | v1.35 | |
| dupRadar / Rsubread | (conda env `dupradar_env`) | |
| R / edgeR | (conda base) | TMM normalization, LOC 매핑 후처리 |
| Trinity (util scripts) | (conda env `trinity_env`) | `run_TMM_scale_matrix.pl`, `PtR`, `run_DE_analysis.pl`, `analyze_diff_expr.pl`, `define_clusters_by_cutting_tree.pl` |
| conda | miniconda3 | 환경: `qc_tools_env`, `dupradar_env`, `trinity_env` |

모든 stage/module은 `nohup bash 스크립트명 &` 형태로 백그라운드 실행하는
것을 전제로 작성되었습니다. 스케줄러(SLURM/PBS)는 사용하지 않습니다.

## 디렉토리 구조

```
.
├── README.md
├── .gitignore
├── scripts/
│   ├── config.sh.example              # 복사해서 config.sh로 만든 뒤 경로 수정
│   ├── list_raw_prefixes.sh           # 사전준비: sample_map.tsv 뼈대 생성
│   ├── setup_screen_references.sh     # 사전준비: FastQ Screen 참조 index 구축
│   ├── exclude_samples.sh             # QC 후 샘플 제외 (samples.txt에서 자동 제거 + 이력 기록)
│   ├── run_pipeline.sh                # stage 1~6 순차 실행 오케스트레이터
│   ├── stage1_qc.sh
│   ├── stage2_trim.sh
│   ├── stage3_screen.sh
│   ├── stage4_align.sh
│   ├── stage5_rseqc.sh
│   ├── stage6_dupradar.sh
│   └── run_dupradar.R
└── modules/                            # StringTie 이후 정량/DE/클러스터링
    ├── module00_prep_annotation.sh          # organelle 제거 + gene라인 제거 + 세미콜론 치환
    ├── module00b_check_folder_names.sh      # StringTie 결과 폴더명(=matrix 컬럼명) 검증
    ├── module04_stringtie.sh                # strand 자동판별 반영, 병렬 처리
    ├── module05_prepDE_after_Stringtie.py3  # python3 호환 prepDE (StringTie 배포본은 python2)
    ├── run_module05.sh
    ├── module06_TPM_Extraction.py
    ├── module07_generate_loc_mapping.sh     # gene_id -> LOC 매핑표
    ├── module07b_generate_name_mapping.sh   # LOC -> gene name (해석용, 선택)
    ├── loc_convert.R                        # count matrix LOC 치환
    ├── loc_convert_tpm.R                    # TPM matrix LOC 치환
    ├── run_loc_conversion.sh                # 위 두 R스크립트 실행 wrapper
    ├── module08_TMM_normalization.sh
    ├── generate_samples_file.sh             # matrix 헤더에서 samples.file 자동 생성
    ├── module09_matrix_process.sh           # PtR + run_DE_analysis.pl + analyze_diff_expr.pl + clustering
    ├── module_validate.sh                   # 4가지 핵심 검증
    ├── troubleshoot_nfs_busy.sh             # NFS ".nfsXXXX busy" 문제 대응
    └── run_all_modules.sh                   # 전체 흐름 참고용 오케스트레이터
```

## 사용 방법

### 0. 초기 설정 (clone 후 최초 1회)

```bash
git clone <this-repo-url>
cd <repo-name>
./setup.sh
```

`setup.sh`가 하는 일:
- `scripts/config.sh.example` → `scripts/config.sh` 복사 (이미 있으면 건드리지 않음)
- **최초 생성 시** `conda info --base`로 `CONDA_SH` 값을 자동 감지해서 채움
- conda/fastqc/hisat2/stringtie 등 필수 프로그램 및 `qc_tools_env`/`dupradar_env`/
  `trinity_env` conda env 존재 여부 확인
- 다음 단계(경로 수정, reference 준비) 안내 출력

이후 `scripts/config.sh`를 열어 `PROJECT_DIR`, `REF_DIR`,
`REF_GENOME_FA`, `REF_GTF`, `HISAT2_INDEX_PREFIX`,
`THREADS`/`THREADS_PER_SAMPLE`/`PARALLEL_SAMPLES`(서버 코어 수에 맞게)를
프로젝트에 맞게 수정하면, 다른 스크립트들은 전부 이 `config.sh` 하나만
보고 동작합니다 (하드코딩된 경로 없음). **`CONDA_SH`는 자동으로
채워지지만, 값이 제대로 들어갔는지 직접 한 번 확인하는 것을
권장합니다** (`grep "CONDA_SH=" scripts/config.sh`) — 관련 문제는
[트러블슈팅](#configsh의-conda_sh-관련-자주-겪는-실수) 참고.

### 1. 사전 준비 (reference 파일 실체 생성 + 샘플 매핑)

`config.sh`에 `REF_GENOME_FA`, `REF_GTF`, `HISAT2_INDEX_PREFIX`,
`RSEQC_BED`, `FASTQ_SCREEN_CONF` 등의 **경로를 적어두는 것과, 그 경로에
실제 파일이 존재하는 것은 별개**입니다. `REF_GENOME_FA`, `REF_GTF`는
직접 준비(다운로드 등)해야 하고, 나머지(index, bed, conf 파일)는 아래
순서대로 명령어를 실행해야 생성됩니다. 이 단계를 건너뛰고 바로
`run_pipeline.sh`를 실행하면 실패합니다.

```bash
cd scripts
```

**1) HISAT2 index 생성** (`HISAT2_INDEX_PREFIX` 경로에 생성됨)
```bash
mkdir -p "$(dirname <config.sh의 HISAT2_INDEX_PREFIX 값>)"
hisat2-build -p 4 <config.sh의 REF_GENOME_FA 값> <config.sh의 HISAT2_INDEX_PREFIX 값>
```

**2) RSeQC용 `genome.bed` 생성** (`RSEQC_BED` 경로에 생성됨, `qc_tools_env` 필요)

AGAT보다 아래 방법을 권장합니다(자세한 이유는
[RSeQC BED 생성](#rseqc-bed-생성) 참고). **원본 GFF3**가 필요합니다
(GTF가 아님):
```bash
conda activate qc_tools_env
gff3ToGenePred -useName <원본genome.gff3> genome.genePred
genePredToBed genome.genePred <config.sh의 RSEQC_BED 값>
conda deactivate
```

**3) FastQ Screen reference + `fastq_screen.conf` 생성** (`SCREEN_REFS_DIR`, `FASTQ_SCREEN_CONF`에 생성됨)

`config.sh`를 자동으로 읽어서 target genome bowtie2 index와
rRNA/EColi/PhiX/Human index를 모두 구축하고 `fastq_screen.conf`까지
자동 생성합니다 (1회성, Human 포함 시 용량/시간이 꽤 큽니다):
```bash
./setup_screen_references.sh
```

**4) 샘플 목록 준비**

```bash
./list_raw_prefixes.sh          # raw fastq 파일명을 보고 sample_map.tsv 뼈대 생성
```
`samples.txt`에 분석할 sample_id를 한 줄에 하나씩 작성하고(**파일 끝
개행 필수** — 없으면 마지막 샘플이 누락됩니다),
`sample_map.tsv`에는 `sample_id<TAB>raw_prefix` 형식으로 raw 파일명과의
매핑을 채워 넣습니다. 샘플 개수나 조건/반복 구성에는 제약이 없습니다.

**작은 테스트 데이터로 위 과정을 한 번에 자동 검증하고 싶다면**
`test_data/run_smoke_test.sh`를 사용하세요 — 1~4번을 전부 자동으로
처리하고 Stage 1~6까지 실행합니다 ([테스트 데이터](#테스트-데이터) 참고).

### 2. Stage 1~6 실행

```bash
cd scripts
nohup ./run_pipeline.sh > run_pipeline.log 2>&1 &     # 전체 실행
./run_pipeline.sh 4 5                                   # stage 4~5만 실행 (align + strand 판별 테스트)
```

각 stage는 개별 실행도 가능합니다(`bash stage4_align.sh` 등). QC
결과를 보고 이상 샘플을 제외하려면:
```bash
./exclude_samples.sh <sample_id> "제외 사유"
```

### 3. Strand 자동판별

Stage 5(RSeQC)가 첫 번째 유효 샘플의 `infer_experiment.py` 결과를 보고
strand(unstranded/forward/reverse)를 자동 판별해서
`results/strand_info.sh`에 저장합니다. Stage 4(HISAT2)와
`modules/module04_stringtie.sh`는 이 파일이 있으면 자동으로 읽어서
`--rna-strandness`, `--fr`/`--rf` 옵션을 적용합니다 (하드코딩 불필요).

### 4. modules/ — 정량 및 DE 분석

**주의**: 아래 module 스크립트들은 `modules/` 디렉토리 안이 아니라,
**정량 작업용 디렉토리**(StringTie가 샘플별 하위폴더를 만들 위치, 예:
`results/07_stringtie/`)에서 실행해야 합니다. `modules/`로 `cd`해서
실행하면 결과 폴더가 `modules/` 안에 잘못 생성됩니다. 아래 예시는
`results/07_stringtie/`에서 실행하는 기준이며, `modules/`까지의
상대경로는 `../../../modules/`입니다 (repo 루트 아래
`scripts/`, `modules/`와 나란히 `test_data/` 등이 있고, 그 아래
`results/07_stringtie/`가 있는 구조 기준 — `PROJECT_DIR`을 repo
바깥에 두는 등 구조가 다르면 경로 깊이도 달라지니, 실제로 실행하기
전에 아래로 파일이 보이는지 먼저 확인하세요):
```bash
ls ../../../modules/module00_prep_annotation.sh
```

```bash
mkdir -p <PROJECT_DIR>/results/07_stringtie   # 정량 작업 디렉토리 (이름은 자유, 이 예시는 07_stringtie)
cd <PROJECT_DIR>/results/07_stringtie

# 0) annotation 전처리 (organelle trans-splicing 유전자 등 문제 해결)
"../../../modules/module00_prep_annotation.sh" \
    /path/to/genomic.gff nuclear_annotation.gff \
    NC_007942.1 NC_020455.1   # chloroplast, mitochondria (예시, 프로젝트마다 확인)

# 1) StringTie 정량 (현재 디렉토리에 샘플별 하위폴더 생성, config.sh 자동 인식)
"../../../modules/module04_stringtie.sh" nuclear_annotation.gff

# 폴더명(= matrix 컬럼명 = samples.file 그룹 라벨) 오타 확인 — 가장 먼저 잡을 것
"../../../modules/module00b_check_folder_names.sh"

# 2) count matrix
"../../../modules/run_module05.sh"

# 3) TPM matrix
python3 "../../../modules/module06_TPM_Extraction.py" expression.TPM.matrix

# 4) gene_id -> LOC 매핑표
"../../../modules/module07_generate_loc_mapping.sh" nuclear_annotation.gff gene_id2loc.tsv

# 5) count/TPM matrix를 LOC 식별자로 치환
"../../../modules/run_loc_conversion.sh" gene_id2loc.tsv gene_count_matrix.csv expression.TPM.matrix

# 6) TMM normalization (Trinity 유틸리티, trinity_env 활성화 필요)
#    주의: 반드시 raw count matrix(genes.counts.matrix)를 입력으로 써야 합니다.
#    이미 정규화된 TPM matrix를 넣으면 lib.size가 모든 샘플에서 동일(1e6)해져
#    TMM 계수가 왜곡됩니다 (자세한 내용은 트러블슈팅 절 참고).
conda activate trinity_env
"../../../modules/module08_TMM_normalization.sh" genes.counts.matrix

# 7) samples.file 생성 (matrix 헤더 기준, 정확한 컬럼명 매칭 보장)
"../../../modules/generate_samples_file.sh" genes.counts.matrix samples.file

# 8) module09_matrix_process.sh: PtR + edgeR DE + 클러스터링
#    P/C(pvalue, fold change) 값을 바꿔가며 여러 번 재실행하는 경우가 많음.
#    재실행마다 이전 결과(edgeR*, diffExpr.*, *.DE_results 등)가 섞일 수 있으니
#    threshold별로 결과를 구분해서 보관하는 걸 권장:
#      mkdir -p runs/P1e-3_C2 && bash ../../../modules/module09_matrix_process.sh -P 1e-3 -C 2 \
#        && mv module09_*.log diffExpr.* clusters_fixed* *.DE_results edgeR* runs/P1e-3_C2/ 2>/dev/null
export R_ENVIRON_USER=/dev/null
"../../../modules/module09_matrix_process.sh" -P 1e-3 -C 2

# 9) module_validate.sh: 최종 검증 (4가지 핵심 체크)
"../../../modules/module_validate.sh" genes.counts.matrix expression.TMM.matrix samples.file
```

전체 흐름을 한 번에 보고 싶다면 `run_all_modules.sh`를 참고하세요
(실제 실행은 프로젝트 상황에 따라 조정 필요).

## 테스트 데이터

외부 다운로드 없이 전체 파이프라인을 수 분 내로 검증할 수 있는 작은
합성(synthetic) 데이터셋을 `test_data/`에서 생성할 수 있습니다.

```bash
cd test_data
./run_smoke_test.sh
```

이 한 줄이 하는 일:
1. `generate_test_data.py`로 작은 genome(7kb, nuclear 2개 + organelle
   유사 scaffold 1개), GTF/GFF3(7개 유전자), 4개 샘플(2그룹×2반복)의
   합성 fastq를 생성
2. 테스트 전용 `config.sh`를 임시로 만들어서 HISAT2 index, `genome.bed`
   자동 구축
3. `run_pipeline.sh 1 6` 실행

성공하면 콘솔에 `Smoke Test 성공` 메시지가 뜨고, `test_data/_smoke_run/`
아래에서 각 stage 결과를 직접 확인할 수 있습니다. 실패하면
`test_data/_smoke_run/results/logs/`를 확인하세요.

데이터만 별도로 생성해서 직접 단계별로 실행해보고 싶다면:
```bash
python3 test_data/generate_test_data.py <원하는 출력 디렉토리>
```
생성되는 것: `ref/genome.fasta`, `ref/genome.gtf`, `ref/genome.gff3`,
`00_raw/*.fastq.gz`, `samples.txt`. 이 경로들을 `config.sh`의
`PROJECT_DIR`, `REF_GENOME_FA`, `REF_GTF` 등에 넣고 [사용
방법](#사용-방법)의 1번(사전 준비)부터 그대로 따라 하면 됩니다.

## 알려진 함정 / 트러블슈팅 (실제로 겪은 것들)

### Annotation 파싱 에러 (module00에서 해결)
- **organelle trans-splicing 유전자**: chloroplast/mitochondria의
  `nad1`, `rps12` 등은 exon이 genome 상 여러 조각·역순·양쪽 strand로
  흩어져 있어(`strand="?"`, `start>end`) StringTie가 파싱 중 죽습니다.
  → `-G` annotation에서 organelle scaffold를 통째로 제외 (HISAT2 정렬
  자체는 organelle 포함 genome으로 했어도 무방).
- **`transcript_id ""`(빈 값)**: GTF의 `gene` feature 라인이 이 값을
  가지고 있어 `"no valid ID found for GFF record"` 에러 발생 →
  `gene` feature 라인 자체를 제거 (exon/transcript/CDS는 그대로 유지되어
  분석 결과에 영향 없음).
- **gene_id 내 세미콜론**: `CYCB1;1` 같은 유전자명이 GTF 속성 구분자(`;`)와
  충돌 → 속성 값 내부의 세미콜론만 언더스코어로 치환.

### RSeQC BED 생성
- **AGAT 변환 비권장**: `agat_convert_sp_gff2bed.pl` 결과는 RSeQC가
  기대하는 표준 BED12와 미묘하게 달라(`thickStart`가 `.`, `itemRgb`가
  콤마 포함 문자열 등) 연쇄적인 파싱 에러가 남.
- **권장 방법**: 원본 GFF3 + `gff3ToGenePred -useName` + `genePredToBed`
  (GTF 대신 GFF3를 직접 쓰면 `transcript_id` 충돌 문제도 회피됨). 여기서도
  organelle trans-splicing 유전자는 사전에 제외해야 함(`invalid strand
  for mRNA`).

### prepDE.py
- StringTie 배포본에 포함된 `prepDE.py`는 python2 문법(`print` statement)
  이라 python2가 없는 서버에서는 실행 불가 → 이 저장소의
  `module05_prepDE_after_Stringtie.py3`(python3 호환 버전) 사용.
- `Error: no GTF files found` → 샘플 폴더/GTF 경로 문제, StringTie가
  **반드시 `-e` 옵션**으로 실행된 GTF여야 coverage가 인식됩니다.

### TPM matrix를 TMM 정규화 입력으로 잘못 쓴 것 (module08, 개념적 버그)
`run_TMM_scale_matrix.pl --matrix`는 **raw read count matrix**를
요구합니다(공식 옵션 설명: `matrix of raw read counts (not
normalized!)`). 그런데 초기 버전은 이미 정규화된 `expression.TPM.matrix`
(module06 결과)를 입력으로 쓰고 있었습니다.
- TPM은 정의상 샘플마다 합이 정확히 1,000,000이 되도록 정규화된 값이라,
  이걸 raw count처럼 넣으면 edgeR `DGEList`가 계산하는 `lib.size`가
  모든 샘플에서 사실상 동일(~1e6)해집니다. TMM이 원래 보정해야 할
  "실제 라이브러리 크기 차이"라는 정보 자체가 입력 단계에서 사라져,
  TPM 위에 다시 정규화를 얹는 이중 정규화가 됩니다.
- 별개로, `run_TMM_scale_matrix.pl`은 정규화 결과를 **STDOUT**으로
  출력하는데, 초기 버전은 이 STDOUT을 버려두고 존재하지도 않는
  `<입력>.TMM_normalized.txt` 파일을 찾다가 항상 에러로 종료했습니다
  (즉 `expression.TMM.matrix`가 한 번도 실제로 생성되지 못했음).
- **수정**: `module08_TMM_normalization.sh`의 기본 입력을
  `genes.counts.matrix`(raw count, LOC 매핑됨)로 바꾸고, STDOUT을
  `expression.TMM.matrix`로 직접 리다이렉트하도록 고쳤습니다. `run_all_modules.sh`도 같이 수정됨.
- module09의 PtR/`run_DE_analysis.pl`은 원래도 `genes.counts.matrix`
  (raw count)를 정확히 쓰고 있었으므로 영향 없었고, `analyze_diff_expr.pl`
  이 쓰는 `expression.TMM.matrix`만 이번 수정으로 정상화됩니다.

### define_clusters_by_cutting_tree.pl의 RData 파일명
가장 많이 겪은 함정입니다.
- `analyze_diff_expr.pl`에 `-P`/`-C`를 명시하지 않으면 기본값(P0.001_C2)
  으로 실행되어 의도한 threshold와 실제 파일명이 어긋남.
- 이 Trinity 버전은 파일명에 **과학적 표기법**을 씁니다:
  `-P 1e-3 -C 2` → `diffExpr.P1e-3_C2.matrix.RData` (`P0.001_C2`
  아님!). `-R` 값을 하드코딩하지 말고 항상 자동 감지할 것:
  ```bash
  RDATA=$(ls diffExpr.*.matrix.RData)
  define_clusters_by_cutting_tree.pl -R "$RDATA" --Ptree 30 --no_column_reordering
  ```
- threshold가 너무 엄격하면(P<0.0001 & C3 등) 유의 유전자가 0개가 되어
  RData 자체가 생성되지 않습니다. 로그가 `"Reading matrix file."`에서
  멈추면 이 경우이니 threshold를 완화하세요 (같은 종 품종 간 비교라면
  P<0.001, C2 정도가 적정인 경우가 많았습니다).

### 병렬 실행 관련
- `conda activate` 내부 스크립트(예: 다른 언어 런타임의 activate hook)가
  `set -u`(bash 엄격 모드)와 충돌해 `unbound variable` 에러를 낼 수 있음
  → conda activate/deactivate 실행 구간만 일시적으로 `set +u` 처리.
- 같은 로그 파일에 여러 프로세스를 동시에 실행하면(중복 실행, 터미널
  꼬임 등) 로그가 섞이고 공유 파일(예: BAM 리스트)에 경합이 생겨
  중복/손상된 산출물이 만들어질 수 있음 → 병렬 처리 시 각 샘플 결과를
  개별 파일로 쓰고, 취합은 전체 job이 끝난 뒤 한 번에 수행.
- NFS 파일시스템에서 `rm: cannot remove ... .nfsXXXX: Device or
  resource busy` → `troubleshoot_nfs_busy.sh` 참고 (다른 세션이 해당
  디렉토리를 붙잡고 있을 수 있음).

### config.sh의 CONDA_SH 관련 (자주 겪는 실수)
- `setup.sh`가 `conda info --base`로 자동 감지해서 채워주지만, 이미
  `config.sh`가 존재하는 상태에서 `setup.sh`를 실행하면 건드리지
  않으므로 **직접 확인이 필요**합니다:
  ```bash
  grep "CONDA_SH=" scripts/config.sh
  ```
  값이 `/path/to/...`(placeholder) 그대로 남아있으면 모든 stage가
  `No such file or directory` 에러로 실패합니다.
- `sed`로 직접 고치다가 따옴표/이스케이프 처리가 꼬여서 엉뚱한 값(예:
  `config.sh` 자기 자신의 경로)이 들어가는 경우가 있었습니다. 수정 후
  반드시 아래로 재확인하세요:
  ```bash
  grep "CONDA_SH=" scripts/config.sh
  # CONDA_SH="/실제/conda/etc/profile.d/conda.sh" 형태인지 확인
  ls -la "$(conda info --base)/etc/profile.d/conda.sh"   # 파일 실존 확인
  ```
- `CONDA_SH` 값과 파일 경로가 맞는데도 `CondaError: Run 'conda init'
  before 'conda activate'`가 뜨면, 그 서버 셸에 `conda init`이 아직
  적용 안 된 것입니다:
  ```bash
  conda init bash
  source ~/.bashrc   # 또는 새 터미널 세션 열기
  ```

### samples.txt / raw fastq 파일명
- bcl2fastq 기본 명명 규칙(예: `01-1_S41_L002_R1_001.fastq.gz`)이
  실제 sample_id와 다르면 `sample_map.tsv`로 매핑해야 하며, 매핑 없이
  진행하면 raw fastq를 못 찾고 전부 SKIP됩니다.
- `samples.txt`가 개행 없이 끝나면 `while read` 루프가 마지막 줄을
  못 읽어 마지막 샘플이 계속 누락됩니다.

## 알려진 한계 / 주의 (본 대두 프로젝트 기준)

- 12품종 중 CMJ076 rep1: RNA degradation으로 라이브러리 복잡도가 낮아
  중복필터 후 unique read가 그룹 내 최저 → DE 검출력이 낮음. 추가
  시퀀싱으로도 개선이 제한적일 수 있음(QC 단계에서 제외 검토 대상).
- CMJ032 rep2: W582와 이상 상관을 보여 sample swap 의심 → 클러스터링
  QC(PtR 상관 히트맵)에서 재점검 필요.
- 같은 종 + 종자 조직(저장단백질 고발현) 특성상 PCA에서 품종 간 분리가
  흐릿하게 나오는 것은 자연스러운 현상. replicate 응집도와 outlier
  여부를 우선 확인할 것.

## 주요 산출물 (modules/ 최종)

| 파일 | 내용 |
|------|------|
| `genes.counts.matrix` | LOC 식별자 raw count (edgeR DE 입력) |
| `expression.TMM.matrix` | LOC 식별자 TMM 정규화 발현값 |
| `gene_id2loc.tsv` | gene_id → LOC 매핑표 |
| `loc2name.tsv` | LOC → gene name (해석용, 선택) |
| `samples.file` | 그룹(품종) ↔ 샘플 매핑 |
| `diffExpr.P*_C*.matrix.RData` | 클러스터링용 R 세션 (파일명은 자동 감지해서 사용) |
| `*.edgeR.DE_results` | pairwise DE 결과 |
| `diffExpr.*.genes_vs_samples_heatmap.pdf` | 유전자×샘플 heatmap |
| `clusters_fixed_P_30.*` | 유전자 클러스터 |

## 검증 체크리스트

`module_validate.sh`가 자동으로 확인하는 항목:
1. count matrix 식별자 중복 없음
2. 모든 라인의 컬럼 수 일치
3. **counts ⊆ TMM** (count matrix의 모든 유전자가 TMM matrix에도 존재) — 가장 중요
4. `samples.file`의 샘플명이 matrix 헤더 컬럼명과 정확히 일치

## 라이선스

별도 라이선스 파일이 없습니다. 조직/실험실 정책에 맞는 라이선스
(예: MIT, 또는 비공개 유지)를 저장소 루트의 `LICENSE` 파일로 추가하세요.

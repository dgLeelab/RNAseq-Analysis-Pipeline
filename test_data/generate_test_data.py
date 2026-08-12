#!/usr/bin/env python3
# =========================================================
# generate_test_data.py
#
# 외부 다운로드 없이, 파이프라인 전체를 빠르게(수 분 내) 검증할 수 있는
# 작은 합성(synthetic) 테스트 데이터셋을 생성합니다.
#
# 생성물:
#   ref/genome.fasta   - 2개의 nuclear 염색체(~3kb) + 1개의 organelle 유사
#                         scaffold(chrC, ~1kb) - module00의 organelle 제외
#                         로직을 테스트하기 위함
#   ref/genome.gtf      - 6개 유전자(각 2 exon), 표준 GTF 포맷
#                         (gene_id "GeneX"; 형태 - module07의 GTF 자동감지 테스트)
#   00_raw/{sample}_R1/R2.fastq.gz - 2그룹 x 2반복 = 4샘플, 샘플당 500 read pair
#                         (reverse-stranded 라이브러리로 시뮬레이션 ->
#                          Stage 5 strand 자동판별 테스트)
#   samples.txt
#
# 사용: python3 generate_test_data.py [출력 디렉토리, 기본: 현재 디렉토리]
# =========================================================

import os
import sys
import random
import gzip

random.seed(42)

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
REF_DIR = os.path.join(OUT_DIR, "ref")
RAW_DIR = os.path.join(OUT_DIR, "00_raw")
os.makedirs(REF_DIR, exist_ok=True)
os.makedirs(RAW_DIR, exist_ok=True)

BASES = "ACGT"


def random_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


def revcomp(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp[b] for b in reversed(seq))


def wrap_fasta(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


# ---------------------------------------------------------
# 1) genome fasta 생성 (nuclear 2개 + organelle 유사 1개)
# ---------------------------------------------------------
chr1_len = 3000
chr2_len = 3000
chrC_len = 1000  # organelle 유사 (module00 제외 테스트용)

chr1_seq = random_seq(chr1_len)
chr2_seq = random_seq(chr2_len)
chrC_seq = random_seq(chrC_len)

genome_fasta_path = os.path.join(REF_DIR, "genome.fasta")
with open(genome_fasta_path, "w") as f:
    f.write(">chr1 synthetic nuclear scaffold\n")
    f.write(wrap_fasta(chr1_seq) + "\n")
    f.write(">chr2 synthetic nuclear scaffold\n")
    f.write(wrap_fasta(chr2_seq) + "\n")
    f.write(">chrC synthetic organelle-like scaffold (for module00 exclusion test)\n")
    f.write(wrap_fasta(chrC_seq) + "\n")

# ---------------------------------------------------------
# 2) GTF 생성: chr1에 4개 유전자, chr2에 2개 유전자, chrC에 1개 유전자
#    각 유전자는 2 exon (100bp + 100bp, 사이 50bp intron)
# ---------------------------------------------------------
genes = []  # (gene_id, chrom, strand, exon1_start, exon1_end, exon2_start, exon2_end)

def place_gene(gene_id, chrom, chrom_len, strand, start):
    e1s = start
    e1e = e1s + 99
    e2s = e1e + 51
    e2e = e2s + 99
    assert e2e < chrom_len, f"gene {gene_id} 배치 실패 (chrom 길이 초과)"
    return (gene_id, chrom, strand, e1s, e1e, e2s, e2e)

genes.append(place_gene("GeneA1", "chr1", chr1_len, "+", 100))
genes.append(place_gene("GeneA2", "chr1", chr1_len, "-", 600))
genes.append(place_gene("GeneA3", "chr1", chr1_len, "+", 1200))
genes.append(place_gene("GeneA4", "chr1", chr1_len, "-", 1800))
genes.append(place_gene("GeneB1", "chr2", chr2_len, "+", 300))
genes.append(place_gene("GeneB2", "chr2", chr2_len, "-", 1000))
genes.append(place_gene("GeneC1", "chrC", chrC_len, "+", 100))  # organelle 유전자

gtf_path = os.path.join(REF_DIR, "genome.gtf")
with open(gtf_path, "w") as f:
    for gid, chrom, strand, e1s, e1e, e2s, e2e in genes:
        tid = f"{gid}.1"
        gene_attr = f'gene_id "{gid}"; transcript_id ""; gene_biotype "protein_coding";'
        tx_attr = f'gene_id "{gid}"; transcript_id "{tid}"; gene_biotype "protein_coding";'
        exon_attr_base = f'gene_id "{gid}"; transcript_id "{tid}"; gene_biotype "protein_coding";'
        # gene 라인 (transcript_id 빈 값 -> module00의 gene라인 제거 로직 테스트)
        f.write(f"{chrom}\ttest\tgene\t{e1s+1}\t{e2e+1}\t.\t{strand}\t.\t{gene_attr}\n")
        f.write(f"{chrom}\ttest\ttranscript\t{e1s+1}\t{e2e+1}\t.\t{strand}\t.\t{tx_attr}\n")
        f.write(f"{chrom}\ttest\texon\t{e1s+1}\t{e1e+1}\t.\t{strand}\t.\t{exon_attr_base} exon_number \"1\";\n")
        f.write(f"{chrom}\ttest\texon\t{e2s+1}\t{e2e+1}\t.\t{strand}\t.\t{exon_attr_base} exon_number \"2\";\n")

# ---------------------------------------------------------
# 3) LOC 매핑 테스트용: GeneID Dbxref 포함한 GFF3 버전도 별도 생성
#    (module07의 GTF/GFF3 자동감지 테스트 - 이 파일은 선택적으로 사용)
# ---------------------------------------------------------
gff3_path = os.path.join(REF_DIR, "genome.gff3")
with open(gff3_path, "w") as f:
    f.write("##gff-version 3\n")
    for i, (gid, chrom, strand, e1s, e1e, e2s, e2e) in enumerate(genes, start=1):
        fake_geneid = 900000 + i
        f.write(
            f"{chrom}\ttest\tgene\t{e1s+1}\t{e2e+1}\t.\t{strand}\t.\t"
            f"ID=gene-{gid};gene_id={gid};Dbxref=GeneID:{fake_geneid};gbkey=Gene;gene={gid}\n"
        )
        tid = f"{gid}.1"
        f.write(
            f"{chrom}\ttest\tmRNA\t{e1s+1}\t{e2e+1}\t.\t{strand}\t.\t"
            f"ID=rna-{tid};Parent=gene-{gid};gbkey=mRNA;gene={gid}\n"
        )
        f.write(
            f"{chrom}\ttest\texon\t{e1s+1}\t{e1e+1}\t.\t{strand}\t.\t"
            f"ID=exon-{tid}-1;Parent=rna-{tid}\n"
        )
        f.write(
            f"{chrom}\ttest\texon\t{e2s+1}\t{e2e+1}\t.\t{strand}\t.\t"
            f"ID=exon-{tid}-2;Parent=rna-{tid}\n"
        )

# ---------------------------------------------------------
# 4) paired-end fastq 시뮬레이션 (2그룹 x 2반복, reverse-stranded)
#    유전자 발현량에 차등을 둬서(GroupA는 GeneA*, GroupB는 GeneB* 고발현)
#    DE 분석 단계까지 의미있는 결과가 나오도록 함
# ---------------------------------------------------------
READ_LEN = 75
N_READS_PER_SAMPLE = 500

def get_transcript_seq(chrom_seq, e1s, e1e, e2s, e2e, strand):
    exon1 = chrom_seq[e1s:e1e + 1]
    exon2 = chrom_seq[e2s:e2e + 1]
    tx = exon1 + exon2
    if strand == "-":
        tx = revcomp(tx)
    return tx

chrom_seqs = {"chr1": chr1_seq, "chr2": chr2_seq, "chrC": chrC_seq}

transcripts = {}
for gid, chrom, strand, e1s, e1e, e2s, e2e in genes:
    transcripts[gid] = get_transcript_seq(chrom_seqs[chrom], e1s, e1e, e2s, e2e, strand)

# 그룹별 상대 발현량 가중치 (GroupA: GeneA* 고발현, GroupB: GeneB* 고발현, GeneC1은 낮게 유지)
def weights_for_group(group):
    w = {}
    for gid in transcripts:
        if group == "groupA":
            w[gid] = 10 if gid.startswith("GeneA") else 2
        else:
            w[gid] = 10 if gid.startswith("GeneB") else 2
        if gid.startswith("GeneC"):
            w[gid] = 1  # organelle 유전자는 항상 낮게
    return w


def add_seq_errors(seq, error_rate=0.01):
    seq = list(seq)
    for i in range(len(seq)):
        if random.random() < error_rate:
            seq[i] = random.choice(BASES)
    return "".join(seq)


def simulate_sample(sample_id, group, r1_path, r2_path):
    weights = weights_for_group(group)
    gene_ids = list(transcripts.keys())
    gene_weights = [weights[g] for g in gene_ids]

    with gzip.open(r1_path, "wt") as r1f, gzip.open(r2_path, "wt") as r2f:
        for i in range(N_READS_PER_SAMPLE):
            gid = random.choices(gene_ids, weights=gene_weights, k=1)[0]
            tx = transcripts[gid]
            if len(tx) < READ_LEN:
                continue
            # fragment 시작 위치 랜덤 (paired-end, insert size ~ READ_LEN*2 가정 -> 여기선 전체 tx 사용)
            max_start = max(0, len(tx) - READ_LEN)
            start = random.randint(0, max_start)
            r1_seq = tx[start:start + READ_LEN]
            # reverse-stranded 라이브러리 시뮬레이션: R2가 sense, R1이 antisense
            r2_seq = revcomp(r1_seq)

            r1_seq = add_seq_errors(r1_seq)
            r2_seq = add_seq_errors(r2_seq)

            qual = "I" * READ_LEN
            read_name = f"@{sample_id}:{i}:{gid}"
            r1f.write(f"{read_name}/1\n{r1_seq}\n+\n{qual}\n")
            r2f.write(f"{read_name}/2\n{r2_seq}\n+\n{qual}\n")


samples = [
    ("groupA_1", "groupA"),
    ("groupA_2", "groupA"),
    ("groupB_1", "groupB"),
    ("groupB_2", "groupB"),
]

samples_txt_path = os.path.join(OUT_DIR, "samples.txt")
with open(samples_txt_path, "w") as f:
    for sample_id, group in samples:
        r1_path = os.path.join(RAW_DIR, f"{sample_id}_R1.fastq.gz")
        r2_path = os.path.join(RAW_DIR, f"{sample_id}_R2.fastq.gz")
        simulate_sample(sample_id, group, r1_path, r2_path)
        f.write(sample_id + "\n")

print(f"[generate_test_data] 완료: {OUT_DIR}")
print(f"  - ref/genome.fasta  ({chr1_len + chr2_len + chrC_len} bp, chr1/chr2/chrC)")
print(f"  - ref/genome.gtf    (7개 유전자: GeneA1-4, GeneB1-2, GeneC1)")
print(f"  - ref/genome.gff3   (동일 유전자, GFF3 포맷 - module07 자동감지 테스트용)")
print(f"  - 00_raw/           (4개 샘플 x R1/R2, 샘플당 {N_READS_PER_SAMPLE} read pairs)")
print(f"  - samples.txt")
print(f"  organelle scaffold: chrC (module00 제외 테스트용)")

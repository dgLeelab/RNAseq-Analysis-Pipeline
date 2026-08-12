#!/usr/bin/env python3
# =========================================================
# module06_TPM_Extraction.py
#
# 각 샘플 폴더의 gene_abundances.tsv(StringTie -A 출력)에서 TPM 컬럼만
# 추출해서 전체 샘플을 합친 expression.TPM.matrix를 생성합니다.
#
# StringTie -A 출력 컬럼(탭 구분):
#   Gene ID, Gene Name, Reference, Strand, Start, End, Coverage, FPKM, TPM
#
# 실행 위치: module04_stringtie.sh 결과 상위 디렉토리
#   (각 {sample}/gene_abundances.tsv 가 있는 곳)
# 사용: python3 module06_TPM_Extraction.py [출력파일명(기본: expression.TPM.matrix)]
# =========================================================

import sys
import os
import glob
import csv
from collections import defaultdict

out_file = sys.argv[1] if len(sys.argv) > 1 else "expression.TPM.matrix"

sample_dirs = sorted(
    d for d in os.listdir(".")
    if os.path.isdir(d) and os.path.isfile(os.path.join(d, "gene_abundances.tsv"))
)

if not sample_dirs:
    sys.stderr.write("[module06] ERROR: gene_abundances.tsv 를 가진 샘플 폴더를 찾지 못했습니다.\n")
    sys.exit(1)

sys.stderr.write(f"[module06] {len(sample_dirs)}개 샘플 발견: {sample_dirs}\n")

# gene_id -> {sample: tpm}
tpm_matrix = defaultdict(dict)
gene_order = []
seen_genes = set()

for sample in sample_dirs:
    path = os.path.join(sample, "gene_abundances.tsv")
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        try:
            gene_idx = header.index("Gene ID")
            tpm_idx = header.index("TPM")
        except ValueError:
            sys.stderr.write(f"[module06] ERROR: {path} 헤더가 예상과 다릅니다: {header}\n")
            sys.exit(1)

        for row in reader:
            if not row:
                continue
            gene_id = row[gene_idx]
            tpm = row[tpm_idx]
            tpm_matrix[gene_id][sample] = tpm
            if gene_id not in seen_genes:
                seen_genes.add(gene_id)
                gene_order.append(gene_id)

with open(out_file, "w", newline="") as fout:
    writer = csv.writer(fout, delimiter="\t")
    writer.writerow(["gene_id"] + sample_dirs)
    for gene_id in gene_order:
        row = [gene_id] + [tpm_matrix[gene_id].get(s, "0") for s in sample_dirs]
        writer.writerow(row)

sys.stderr.write(f"[module06] 완료: {out_file} ({len(gene_order)}개 유전자 x {len(sample_dirs)}개 샘플)\n")

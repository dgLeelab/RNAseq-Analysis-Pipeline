#!/usr/bin/env bash
# =========================================================
# module00_prep_annotation.sh
# StringTie(-G)용 annotation 전처리 (Trans-splicing organelle 유전자,
# 빈 transcript_id, gene_id 내 세미콜론 문제를 모두 해결)
#
# 배경 (실제 겪은 에러들):
#   1) organelle(chloroplast/mitochondria) trans-splicing 유전자
#      (예: mitochondria nad1)는 exon이 genome 상 역순/여러 조각으로
#      흩어져 있어 strand가 "?"로 표기되거나 end<start인 경우가 있어
#      StringTie가 파싱 중 죽음 ("invalid strand for mRNA" 등)
#   2) GTF의 gene feature 라인은 transcript_id ""(빈 값)으로 되어 있어
#      StringTie가 "no valid ID found for GFF record" 에러를 냄
#   3) 일부 gene_id/gene 값 자체에 세미콜론이 포함된 경우
#      (예: CYCB1;1) GTF 속성 구분자(;)와 충돌하여 파싱 에러 유발
#
# 사용: ./module00_prep_annotation.sh <원본 GTF 또는 GFF3> <출력 파일>
#       [organelle_scaffold_id ...]
# 예:   ./module00_prep_annotation.sh \
#         GCF_000004515.6_Glycine_max_v4.0_genomic.gff \
#         Gmax_v4.0_nuclear.gff \
#         NC_007942.1 NC_020455.1
# =========================================================
set -euo pipefail

INFILE="${1:?사용법: $0 <input.gtf|gff> <output.gff> [organelle_scaffold_id ...]}"
OUTFILE="${2:?출력 파일명을 지정하세요}"
shift 2 || true
ORGANELLE_IDS=("$@")

if [[ ${#ORGANELLE_IDS[@]} -eq 0 ]]; then
    echo "[module00] organelle scaffold ID가 지정되지 않았습니다."
    echo "           (chloroplast/mitochondria 등 trans-splicing 유전자가 있는 organelle genome)"
    echo "           대두 예시: NC_007942.1(chloroplast) NC_020455.1(mitochondria)"
fi

echo "[module00] 1) organelle scaffold 제거"
if [[ ${#ORGANELLE_IDS[@]} -gt 0 ]]; then
    python3 - "${INFILE}" /tmp/module00_step1.gff "${ORGANELLE_IDS[@]}" << 'PYEOF'
import sys
infile, outfile = sys.argv[1], sys.argv[2]
organelle_ids = set(sys.argv[3:])
with open(infile) as fin, open(outfile, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            skip = any(f" {oid}" in line or line.startswith(f"##sequence-region {oid}") for oid in organelle_ids)
            if not skip:
                fout.write(line)
            continue
        fields = line.split("\t")
        if len(fields) < 1 or fields[0] not in organelle_ids:
            fout.write(line)
PYEOF
else
    cp "${INFILE}" /tmp/module00_step1.gff
fi

echo "[module00] 2) gene feature 라인 제거 + gene_id/gene 값 내 세미콜론 치환"
python3 - /tmp/module00_step1.gff "${OUTFILE}" << 'PYEOF'
import re, sys
infile, outfile = sys.argv[1], sys.argv[2]

def fix_line(line):
    if line.startswith("#"):
        return line
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 9:
        return line
    if fields[2] == "gene":
        return None  # gene feature 라인 제거 (transcript_id "" 문제 회피)
    attrs = fields[8]
    # 큰따옴표로 감싸인 속성 값(GTF) 또는 속성 값 내부의 세미콜론을 언더스코어로 치환
    fixed_attrs = re.sub(r'"[^"]*"', lambda m: m.group(0).replace(";", "_"), attrs)
    fields[8] = fixed_attrs
    return "\t".join(fields) + "\n"

with open(infile) as fin, open(outfile, "w") as fout:
    for line in fin:
        fixed = fix_line(line)
        if fixed is not None:
            fout.write(fixed)
PYEOF

rm -f /tmp/module00_step1.gff

echo "[module00] 완료: ${OUTFILE}"
echo "[module00] 검증:"
for id in "${ORGANELLE_IDS[@]}"; do
    cnt=$(grep -c "${id}" "${OUTFILE}" || true)
    echo "  ${id} 잔존 라인: ${cnt} (0 이어야 정상)"
done
cnt_empty_tid=$(grep -c 'transcript_id ""' "${OUTFILE}" || true)
echo "  transcript_id \"\" (빈값) 잔존: ${cnt_empty_tid} (0 이어야 정상)"

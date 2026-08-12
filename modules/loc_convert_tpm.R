#!/usr/bin/env Rscript
# =========================================================
# loc_convert_tpm.R
# TPM matrix의 gene_id를 LOC 식별자로 치환
#
# 실행: R_ENVIRON_USER=/dev/null Rscript loc_convert_tpm.R \
#         [gene_id2loc.tsv] [expression.TPM.matrix] [expression.TPM.matrix.loc]
# =========================================================
args <- commandArgs(trailingOnly = TRUE)
map_file  <- if (length(args) >= 1) args[1] else "gene_id2loc.tsv"
tpm_file  <- if (length(args) >= 2) args[2] else "expression.TPM.matrix"
out_file  <- if (length(args) >= 3) args[3] else "expression.TPM.matrix.loc"

map <- read.table(map_file, sep="\t", header=F, quote="", comment.char="", stringsAsFactors=F)
m <- setNames(map[,2], map[,1])
d <- read.table(tpm_file, header=T, sep="\t", quote="", comment.char="", check.names=F, stringsAsFactors=F)
d[,1] <- ifelse(d[,1] %in% names(m), m[d[,1]], d[,1])
cat("duplicate count:", sum(duplicated(d[,1])), "\n")
write.table(d, out_file, sep="\t", quote=F, row.names=F)

#!/usr/bin/env Rscript
# =========================================================
# loc_convert.R
# gene_count_matrix.csv 의 gene_id 를 LOC 식별자로 치환
#
# 실행: R_ENVIRON_USER=/dev/null Rscript loc_convert.R \
#         [gene_id2loc.tsv] [gene_count_matrix.csv] [gene_count_matrix.loc.csv]
# R_ENVIRON_USER=/dev/null 필수 (.Renviron의 R_LIBS가 S7 패키지와 충돌)
# =========================================================
args <- commandArgs(trailingOnly = TRUE)
map_file    <- if (length(args) >= 1) args[1] else "gene_id2loc.tsv"
count_file  <- if (length(args) >= 2) args[2] else "gene_count_matrix.csv"
out_file    <- if (length(args) >= 3) args[3] else "gene_count_matrix.loc.csv"

map <- read.table(map_file, sep="\t", header=F, quote="", comment.char="", stringsAsFactors=F)
m <- setNames(map[,2], map[,1])
d <- read.csv(count_file, header=T, check.names=F, stringsAsFactors=F)
gid <- sub("\\|.*", "", d[,1])
d[,1] <- ifelse(gid %in% names(m), m[gid], gid)
cat("duplicate count:", sum(duplicated(d[,1])), "\n")
write.csv(d, out_file, row.names=F, quote=F)

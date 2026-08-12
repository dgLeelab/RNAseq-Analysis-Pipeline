#!/usr/bin/env Rscript
# =========================================================
# run_dupradar.R
# 사용: Rscript run_dupradar.R <markdup.bam> <annotation.gtf> <out_prefix>
# 실행 env: dupradar_env
# =========================================================

suppressMessages(library(dupRadar))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: run_dupradar.R <markdup.bam> <annotation.gtf> <out_prefix> <stranded:0|1|2>")
}

bam_file    <- args[1]
gtf_file    <- args[2]
out_prefix  <- args[3]
stranded    <- as.integer(args[4])   # 0=unstranded, 1=fwd, 2=rev (Stage 5 RSeQC 자동판별 결과)
paired_end  <- TRUE
threads    <- 8

# ---- dupRadar 실행 --------------------------------------------
dm <- analyzeDuprates(
    bam        = bam_file,
    gtf        = gtf_file,
    stranded   = stranded,
    paired     = paired_end,
    threads    = threads
)

# ---- 결과 저장 --------------------------------------------------
write.table(
    dm,
    file      = paste0(out_prefix, "_dupmatrix.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
)

# duplication rate vs expression 시각화 (2D density plot)
pdf(paste0(out_prefix, "_dupradar_plot.pdf"), width = 6, height = 6)
duprateExpDensPlot(DupMat = dm)
title(basename(out_prefix))
dev.off()

cat("Done:", out_prefix, "\n")

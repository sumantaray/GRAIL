# Rscript scripts/simulate.R <output folder> [simulation seed]
# Reconstructed Table 1 designs. Keep the measured cell expression and branch labels.
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
root <- args[[1]]
seed <- if (length(args) > 1) as.integer(args[[2]]) else 42L
if (any(file.exists(file.path(root, paste0("P", 1:7), "matrix.mtx"))))
    stop("Simulation files already exist; choose a fresh output directory.")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
session <- capture.output(sessionInfo())
session <- sub("BLAS/LAPACK: .*/", "BLAS/LAPACK: ", session)
writeLines(trimws(session, which = "right"), file.path(root, "R-session.txt"))
for (i in seq_len(7)) {
    run_seed <- seed + i - 1L
    target <- switch(as.character(i), "2" = c(100, 800), "6" = c(300, 350),
                     "7" = c(500, 300), c(300, 400))
    opt <- list(nGenes = 3000L, batchCells = 1400L, seed = run_seed,
                group.prob = c(0.3, 0.3, 0.4), path.from = c(0L, 1L, 2L),
                de.prob = 0.4, dropout.type = "experiment", dropout.mid = 0, verbose = FALSE)
    if (i == 2) opt$group.prob <- c(0.1, 0.1, 0.8)
    if (i == 3) { opt$de.prob <- 0.1; opt$de.facLoc <- 0.2; opt$de.facScale <- 0.4 }
    if (i == 4) { opt$de.prob <- 0.7; opt$de.facLoc <- 1.0; opt$de.facScale <- 1.2 }
    if (i == 5) opt$dropout.mid <- 2
    if (i == 6) { opt$group.prob <- c(0.3, 0.35, 0.35); opt$path.from <- c(0L, 1L, 1L) }
    if (i == 7) { opt$batchCells <- c(500L, 200L, 300L); opt$batch.facLoc <- 1.5 }
    sce <- do.call(splatSimulatePaths, opt)
    meta <- as.data.frame(colData(sce))
    path <- as.integer(gsub("[^0-9]", "", as.character(meta$Group)))
    batch <- as.integer(gsub("[^0-9]", "", as.character(meta$Batch)))
    first <- if (i == 7) which(batch == 1) else which(path == 1)
    last <- if (i == 7) which(batch == 3) else if (i == 6) which(path %in% c(2, 3)) else which(path == 3)
    if (length(first) < target[1] || length(last) < target[2]) stop("Insufficient simulated population.")
    set.seed(run_seed + 10000L)
    ia <- sort(sample(first, target[1]))
    ib <- sort(sample(last, target[2]))
    idx <- c(ia, ib)
    info <- data.frame(cell = colnames(sce)[idx], t = rep(c(0,1), target),
                       simulator_path = path[idx], batch = batch[idx])
    if (i == 6) info$branch <- c(rep(NA, length(ia)), paste0("branch_", path[ib] - 1L))
    out <- file.path(root, paste0("P", i))
    dir.create(out)
    writeMM(as(counts(sce)[,idx], "CsparseMatrix"), file.path(out, "matrix.mtx"))
    write.table(info, file.path(out, "cells.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    writeLines(rownames(sce), file.path(out, "genes.txt"))
    writeLines(capture.output(dput(opt)), file.path(out, "parameters.R"))
}

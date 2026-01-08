install.packages("lattice", repos = "https://cloud.r-project.org") # might be necessary
install.packages("https://cloud.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-3.tar.gz")#, repos = "https://cloud.r-project.org")
install.packages("https://cloud.r-project.org/src/contrib/Archive/MASS/MASS_7.3-59.tar.gz")#, repos = "https://cloud.r-project.org")
#install.packages("MASS", repos = "https://cloud.r-project.org")
#install.packages("Matrix", repos = "https://cloud.r-project.org")
install.packages("mgcv", repos = "https://cloud.r-project.org")
install.packages("dplyr", repos = "https://cloud.r-project.org")
library(dplyr)
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("bsseq")
BiocManager::install("DMRcate")
BiocManager::install("BiocParallel")
library(bsseq)
library(DMRcate)
library(BiocParallel)

wgbs_counts_t <- readRDS("/researchers/caitlin.page/wgbs_counts_t.rds")

wgbs_counts_t <- wgbs_counts_t %>% filter(!seqnames %in% c("chrX", "chrY", "chrM"))
m_df <- wgbs_counts_t[,c(5,7,9,14,16,18)]
colnames(m_df) <- colnames(m_df) %>% gsub("_M", "", .)

coverage_df <- wgbs_counts_t[,c(6,8,10,15,17,19)]
colnames(coverage_df) <- colnames(coverage_df) %>% gsub("_C", "", .)

sample_stuff <- cbind(Type = c("cd4", "cd4", "cd4", "cd8", "cd8", "cd8")) %>% data.frame()
rownames(sample_stuff) <- colnames(m_df)

## make BSseq object
bseq_cd4_vs_cd8_all <- BSseq(M = as.matrix(m_df),
                             Cov = as.matrix(coverage_df), pData = sample_stuff,
                             chr = wgbs_counts_t$seqnames, pos = as.numeric(wgbs_counts_t$start),
                             sampleNames = colnames(m_df))

## BSmooth by sample - so it doesn't crash
bsmooth1 <- BSmooth(BSseq = bseq_cd4_vs_cd8_all[,1], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)

bsmooth2 <- BSmooth(BSseq = bseq_cd4_vs_cd8_all[,2], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)
bsmooth3 <- BSmooth(BSseq = bseq_cd4_vs_cd8_all[,3], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)

bsmooth4 <- BSmooth(BSseq = bseq_cd4_vs_cd8_all[,4], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)
bsmooth5 <- BSmooth(BSseq = bseq_cd4_vs_cd8_all[,5], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)
bsmooth6 <- BSmooth(BSseq = bseq_cd4_vs_cd8_all[,6], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)

bsmooth_obj <- combine(bsmooth1, bsmooth2, bsmooth3, bsmooth4, bsmooth5, bsmooth6)
## filter
bseq.cov <- getCoverage(bsmooth_obj)
keep <- which(rowSums(bseq.cov[, bsmooth_obj$Type == "cd4"] >= 2) >= 2 &
                rowSums(bseq.cov[, bsmooth_obj$Type == "cd8"] >= 2) >= 2)

bsmooth_obj <- bsmooth_obj[keep,]
saveRDS(bsmooth_obj, "../researchers/caitlin.page/bsmooth_obj_t.rds")

tissue <- factor(c("cd4", "cd4", "cd4", "cd8", "cd8", "cd8"), levels = c("cd4", "cd8"))
tissue
#Regular matrix design
design <- model.matrix(~tissue)
design

methdesign <- modelMatrixMeth(design)
methdesign

seq_annot <- sequencing.annotate(bsmooth_obj, methdesign, all.cov = TRUE,
                                 coef = 8, fdr=0.05)

seq_dmrcate.res <- dmrcate(seq_annot, C=2, min.cpgs = 5)

dmrcate_seq_dmr <- extractRanges(seq_dmrcate.res, genome = "hg19")
dmrcate_seq_dmr <- dmrcate_seq_dmr %>% data.frame()
dmrcate_seq_anno <- seq_annot@ranges %>% data.frame()
saveRDS(seq_annot, "/researchers/caitlin.page/dmrcate_seq_anno_t.rds")
saveRDS(dmrcate_seq_dmr, "/researchers/caitlin.page/dmrcate_seq_dmr_t.rds")

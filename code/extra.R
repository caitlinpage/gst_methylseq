
library(DMRcate)
library(dplyr)

library(bsseq)
library(BiocParallel)
library(limma)
library(edgeR)

wgbs_counts <- readRDS("/researchers/caitlin.page/phd_gst/output/wgbs_counts.rds")

beta_all <- wgbs_counts
beta_all <- beta_all %>% filter(!seqnames %in% c("chrX", "chrY", "chrM"))
m_df <- beta_all[,c(5,7,9,14,16,18)]
colnames(m_df) <- colnames(m_df) %>% gsub("_M", "", .)

coverage_df <- beta_all[,c(6,8,10,15,17,19)]
colnames(coverage_df) <- colnames(coverage_df) %>% gsub("_C", "", .)

sample_stuff <- cbind(Type = c("nk", "nk", "nk", "b", "b", "b")) %>% data.frame()
rownames(sample_stuff) <- colnames(m_df)

## make BSseq object
bseq_nk_vs_b_all <- BSseq(M = as.matrix(m_df),
                          Cov = as.matrix(coverage_df), pData = sample_stuff,
                          chr = beta_all$seqnames, pos = as.numeric(beta_all$start),
                          sampleNames = colnames(m_df))

bsmooth_obj_1 <- BSmooth(BSseq = bseq_nk_vs_b_all[1,], BPPARAM = MulticoreParam(workers = 12), verbose = TRUE)
bsmooth_obj_2 <- BSmooth(BSseq = bseq_nk_vs_b_all[2,], BPPARAM = MulticoreParam(workers = 12), verbose = TRUE)
bsmooth_obj_3 <- BSmooth(BSseq = bseq_nk_vs_b_all[3,], BPPARAM = MulticoreParam(workers = 12), verbose = TRUE)
bsmooth_obj_4 <- BSmooth(BSseq = bseq_nk_vs_b_all[4,], BPPARAM = MulticoreParam(workers = 12), verbose = TRUE)
bsmooth_obj_5 <- BSmooth(BSseq = bseq_nk_vs_b_all[5,], BPPARAM = MulticoreParam(workers = 12), verbose = TRUE)
bsmooth_obj_6 <- BSmooth(BSseq = bseq_nk_vs_b_all[6,], BPPARAM = MulticoreParam(workers = 12), verbose = TRUE)

bsmooth_obj <- bsseq::combine(bsmooth_obj_1, bsmooth_obj_2, bsmooth_obj_3, bsmooth_obj_4, bsmooth_obj_5, bsmooth_obj_6)

bseq.cov <- getCoverage(bsmooth_obj)
keep <- which(rowSums(bseq.cov[, bsmooth_obj$Type == "nk"] >= 2) >= 2 &
                rowSums(bseq.cov[, bsmooth_obj$Type == "b"] >= 2) >= 2)
length(keep)

bsmooth_obj <- bsmooth_obj[keep,]

tissue <- factor(pData(bsmooth_obj)$Type)

#Regular matrix design
design <- model.matrix(~tissue)

methdesign <- edgeR::modelMatrixMeth(design)

seq_annot <- sequencing.annotate(bsmooth_obj, methdesign, all.cov = TRUE,
                                 coef = 8, fdr=0.05)

dmrcate.res <- dmrcate(seq_annot, C=2, min.cpgs = 5)

dmrcate_seq_dmr <-  extractRanges(dmrcate.res, genome="hg19")

saveRDS(seq_annot, "/researchers/caitlin.page/phd_gst/output/dmrcate_seq.rds")

saveRDS(dmrcate_seq_dmr, "/researchers/caitlin.page/phd_gst/output/dmrcate_dmrs.rds")

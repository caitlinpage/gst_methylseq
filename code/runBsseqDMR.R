library(rtracklayer)
library(plyranges)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(bsseq)
library(BiocParallel)

# bsseq vignette
# https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq_analysis.html

# get the positions of cg sites
seq_names <- seqnames(BSgenome.Hsapiens.UCSC.hg19)[1:25]
anno_seq <- lapply(seq_names, function(x) {
  cbind(data.frame(matchPattern("CG", BSgenome.Hsapiens.UCSC.hg19[[x]])),
        seqnames = x
  )
}) %>%
  dplyr::bind_rows()
anno_seq <- anno_seq %>% mutate(pos = paste0(seqnames, "-", start)) %>%
  relocate(pos, seqnames)


# NK counts
fname <- "data/wgbs/GSM5652299_Blood-NK-Z000000TM.beta"
N <- file.info(fname)$size
beta_nk <- data.frame(matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE))
fname <- "data/wgbs/GSM5652300_Blood-NK-Z000000U1.beta"
N <- file.info(fname)$size
beta_nk <- cbind(beta_nk, data.frame(matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)))
fname <- "data/wgbs/GSM5652301_Blood-NK-Z000000UF.beta"
N <- file.info(fname)$size
beta_nk <- cbind(beta_nk, data.frame(matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)))
colnames(beta_nk) <- c("M_TM", "C_TM", "M_U1", "C_U1", "M_UF", "C_UF")
beta_nk <- cbind(anno_seq[,1:4], beta_nk)

# B counts
fname <- "data/wgbs/GSM5652316_Blood-B-Z000000TX.beta"
N <- file.info(fname)$size
beta_b <- data.frame(matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE))
fname <- "data/wgbs/GSM5652317_Blood-B-Z000000UB.beta"
N <- file.info(fname)$size
beta_b <- cbind(beta_b, data.frame(matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)))
fname <- "data/wgbs/GSM5652318_Blood-B-Z000000UR.beta"
N <- file.info(fname)$size
beta_b <- cbind(beta_b, data.frame(matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)))
colnames(beta_b) <- c("M_TX", "C_TX", "M_UB", "C_UB", "M_UR", "C_UR")
beta_b <- cbind(anno_seq[,1:4], beta_b)

# add unmeth counts
beta_nk <- beta_nk %>% mutate(Un_TM = C_TM - M_TM, Un_U1 = C_U1 - M_U1, Un_UF = C_UF - M_UF)
beta_b <- beta_b %>% mutate(Un_TX = C_TX - M_TX, Un_UB = C_UB - M_UB, Un_UR = C_UR - M_UR)

# combine
colnames(beta_b)[5:13] <- paste0("b_", colnames(beta_b)[5:13])
colnames(beta_nk)[5:13] <- paste0("nk_", colnames(beta_nk)[5:13])

beta_all <- cbind(beta_nk, beta_b[,5:13])

saveRDS(beta_all, "output/wgbs_counts.rds")

#####

# bsseq dmr analysis
## set up
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

## BSmooth by sample - so it doesn't crash
bsmooth1 <- BSmooth(BSseq = bseq_nk_vs_b_all[,1], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)

bsmooth2 <- BSmooth(BSseq = bseq_nk_vs_b_all[,2], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)
bsmooth3 <- BSmooth(BSseq = bseq_nk_vs_b_all[,3], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)

bsmooth4 <- BSmooth(BSseq = bseq_nk_vs_b_all[,4], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)
bsmooth5 <- BSmooth(BSseq = bseq_nk_vs_b_all[,5], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)
bsmooth6 <- BSmooth(BSseq = bseq_nk_vs_b_all[,6], BPPARAM = MulticoreParam(workers = 6), verbose = TRUE)

bseq_nk_vs_b_all_fit <- bsseq::combine(bsmooth1, bsmooth2, bsmooth3, bsmooth4, bsmooth5, bsmooth6)

## filter
bseq.cov <- getCoverage(bseq_nk_vs_b_all_fit)
keep <- which(rowSums(bseq.cov[, bseq_nk_vs_b_all_fit$Type == "nk"] >= 2) >= 2 &
                rowSums(bseq.cov[, bseq_nk_vs_b_all_fit$Type == "b"] >= 2) >= 2)
length(keep)

bseq_nk_vs_b_all_fit <- bseq_nk_vs_b_all_fit[keep,]

## testing
bseq.tstat <- BSmooth.tstat(bseq_nk_vs_b_all_fit, group1 = c("nk_TM", "nk_U1", "nk_UF"),
                            group2 = c("b_TX", "b_UB", "b_UR"),
                            estimate.var = "group2",
                            local.correct = TRUE,
                            verbose = TRUE, mc.cores = 6)

bseq_all_res <- bseq.tstat@stats %>% data.frame() %>%
  mutate(seqnames = data.frame(bseq.tstat@gr)$seqnames,
         start = data.frame(bseq.tstat@gr)$start,
         end = start, position = paste0(seqnames, "-", start)) %>%
  relocate(position)
saveRDS(bseq_all_res, "output/bsseq_res.rds")

## dmrs
bsseq_dmrs <- dmrFinder(bseq.tstat)
# filter the dmrs
bsseq_dmrs <- subset(bsseq_dmrs, n >= 3 & abs(meanDiff) >= 0.1)
colnames(bsseq_dmrs)[1] <- "seqnames"
bsseq_dmrs <- bsseq_dmrs %>%
  mutate(position = paste0(seqnames, "-", start)) %>% relocate(position)
bsseq_dmrs$rank <- 1:nrow(bsseq_dmrs)
saveRDS(bsseq_dmrs, "output/bsseq_dmrs.rds")

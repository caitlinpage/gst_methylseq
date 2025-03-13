
library(dplyr)


## Run edgeR method
library(edgeR)
#https://f1000research.com/articles/6-2055

#anno_overlap



#beta_all <- wgbs_counts[,c(1:4)]

edge_matrix <- wgbs_counts[,c(5,11,7,12,9,13,14,20,16,21,18,22)]

colnames(edge_matrix) <- c(
  "NK1.Me", "NK1.Un", "NK2.Me", "NK2.Un", "NK3.Me", "NK3.Un",
  "B1.Me", "B1.Un", "B2.Me", "B2.Un", "B3.Me", "B3.Un")
rownames(edge_matrix) <- wgbs_counts$pos
#edge_matrix <- edge_matrix[anno_overlap$position,]

edge_matrix <- edge_matrix[!grepl("chrX", rownames(edge_matrix)),]
edge_matrix <- edge_matrix[!grepl("chrY", rownames(edge_matrix)),]

sample <- gl(6,2,12)
methylation <- gl(2,1,12, labels=c("Me","Un"))
condition <- gl(2,6,12, labels=c("NK", "B"))
edge_sample_info <- cbind(sample, methylation, condition)
rownames(edge_sample_info) <- colnames(edge_matrix)
data.frame(edge_sample_info)

#anno_overlap

edge_dge <- DGEList(as.matrix(edge_matrix), samples = data.frame(edge_sample_info),
                    group = c(
                      "NK", "NK", "NK", "NK", "NK", "NK",
                      "B", "B", "B", "B", "B", "B"))#, genes = filter(anno_overlap, !seqnames %in% c("chrX", "chrY"))[,c("position", "seqnames", "start", "Name")])

edge_coverage <- edge_dge$counts[, c(1,3,5,7,9,11)] + edge_dge$counts[, c(2,4,6,8,10,12)]
head(edge_coverage)
keep <- rowSums(edge_coverage >= 8) == 6

edge_dge <- edge_dge[keep,,keep.lib.sizes=FALSE]

edge_dge$samples

edge_total_lib_size <- filter(edge_dge$samples, methylation==1)$lib.size +
  filter(edge_dge$samples, methylation==2)$lib.size

edge_dge$samples$lib.size <- rep(edge_total_lib_size, each=2)
edge_dge$samples

#edge_Me <- edge_dge$counts[, c(1,3,5,7,9,11)]
#edge_Un <- edge_dge$counts[, c(2,4,6,8,10,12)]
#edge_M <- log2(edge_Me + 2) - log2(edge_Un + 2)

#plotMDS(edge_M)

edge_samples <- cbind(Group = c("NK", "NK", "NK", "B", "B", "B")) %>% data.frame()
rownames(edge_samples) <- c("NK1", "NK2", "NK3", "B1", "B2", "B3")
edge_samples <- edge_samples %>% mutate(Group = as.factor(Group))
#edge_samples

edge_design <- model.matrix(~0 + Group, data = edge_samples)
colnames(edge_design) <- c("B", "NK")

edge_design <- modelMatrixMeth(edge_design)

edge_dge <- estimateDisp(edge_dge, edge_design, trend="none")
#edge_dge$common.dispersion
#edge_dge$prior.df

edge_fit <- glmFit(edge_dge, edge_design)
edge_contrasts <- makeContrasts(NK - B,
                                levels=edge_design)
edge_lrt <- glmLRT(edge_fit, contrast=edge_contrasts)

#plotMD(edge_lrt)

# by chrom analysis

#unfactor(unique(edge_dge$genes$seqnames))

#ChrIndices <- list()
#for (a in unfactor(unique(edge_dge$genes$seqnames))) ChrIndices[[a]] <- which(edge_dge$genes$seqnames==a)

#fry(edge_dge, index=ChrIndices, design=edge_design, contrast=edge_contrasts)


#the mixed vals are signif but not the indiv ones so i think that means there's nothign really going on

############
# the results


edge_res <- topTags(edge_lrt, n = Inf)

edge_res <- edge_res$table %>% mutate(signif.05 = ifelse(FDR <= 0.05, TRUE, FALSE))



###
# save file

write.table(edge_res, "~/Desktop/r/phd/2_ontology_bias/output/seq_edger_res.txt")


#####################################################

## run DMRcate
# start with bsseq smoothed object
library(DMRcate)
library(bsseq)
bsmooth_obj <- readRDS("/researchers/caitlin.page/phd_gst/output/bsmooth_obj.rds")
bseq.cov <- getCoverage(bsmooth_obj)
keep <- which(rowSums(bseq.cov[, bsmooth_obj$Type == "nk"] >= 2) >= 2 &
                rowSums(bseq.cov[, bsmooth_obj$Type == "b"] >= 2) >= 2)
length(keep)

bsmooth_obj <- bsmooth_obj[keep,]

tissue <- factor(pData(bsmooth_obj)$Type)
tissue
#Regular matrix design
design <- model.matrix(~tissue)

methdesign <- edgeR::modelMatrixMeth(design)

seq_annot <- sequencing.annotate(bsmooth_obj, methdesign, all.cov = TRUE,
                                 coef = 8, fdr=0.05)

dmrcate.res <- dmrcate(seq_annot, C=2, min.cpgs = 5)

dmrcate_seq_dmr <-  extractRanges(dmrcate.res, genome="hg19")

saveRDS(seq_annot, "/researchers/caitlin.page/phd_gst/output/dmrcate_seq_anno.rds")
saveRDS(dmrcate_seq_dmr, "/researchers/caitlin.page/phd_gst/output/dmrcate_seq_dmr.rds")


######################################################

## run DSS

library(DSS)
library(bsseq)
library(BiocParallel)
# start with bsseq object (unsure if smoothed)
bseq_nk_vs_b_all <- readRDS("~/Desktop/r/phd/2_ontology_bias/output/bseq_obj_all.rds")

dss_dmlTest_all <- DMLtest(bseq_nk_vs_b_all,
                                  group1=c("nk_TM", "nk_U1", "nk_UF"),
                                  group2 = c("b_TX", "b_UB", "b_UR"),
                               smoothing = TRUE,
                               ncores = 8)
dss_dmrs_pos_filt <- callDMR(dss_dmlTest_nk_vs_b_pos_filt, p.threshold=0.01)
dss_dmrs <- callDMR(dss_dmlTest_nk_vs_b, p.threshold=0.01)

#####################################################

## run MethylKit
# abandoned - was too annoying
# also tried in [1_methods/analysis/2_trumethbs.Rmd]
# kept having timeout issues - would require cluster (which I didn't use)


## run RnBeads

## run betaHMM

## run bumphunter

## run other things (and fast)

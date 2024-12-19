library(plyranges)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ChIPseqSpikeInFree)
library(biomaRt)


plotBiasGenomeBins <- function(counts, dmr, bin_size=2000) {
  genome_bins <- GenerateBins("hg19", binSize = bin_size)
  colnames(genome_bins) <- c("seqnames", "start", "end")
  genome_bins <- genome_bins %>% mutate(bin_pos = paste0(seqnames, "-", start), num = 1:n())

  counts_per_bin <- find_overlaps(as_granges(genome_bins), as_granges(counts)) %>% data.frame()

  cg_per_bin <- counts_per_bin[,6:7] %>% group_by(num, bin_pos) %>% summarise(n_region = n())
  cg_per_bin <- cg_per_bin %>% ungroup() %>% data.frame()

  genome_bins <- genome_bins %>% mutate(num_cg = cg_per_bin[match(.$num, cg_per_bin$num), "n_region"])
  genome_bins[is.na(genome_bins)] <- 0
  genome_bins <- genome_bins %>% filter(!seqnames %in% c("chrX", "chrY", "chrM"))


  bin_dmr_overlap <- find_overlaps(as_granges(genome_bins), as_granges(dmr)) %>% data.frame()
  bin_dmr_overlap <- bin_dmr_overlap %>% group_by(bin_pos) %>%
    summarise(n_dmr = n()) %>% ungroup() %>% data.frame()
  genome_bins <- genome_bins %>% mutate(n_dmr = bin_dmr_overlap[match(.$bin_pos, bin_dmr_overlap$bin_pos), "n_dmr"])

  genome_bins[is.na(genome_bins)] <- 0

  genome_bins <- genome_bins %>% mutate(has_dmr = ifelse(n_dmr == 0, FALSE, TRUE))

  plot <- genome_bins %>%
    group_by(num_cg) %>%
    mutate(num_bins_same_cg = n()) %>%
    group_by(num_cg, has_dmr) %>%
    mutate(num_per_dmr = n(), prop = n()/num_bins_same_cg) %>%
    distinct(num_cg, num_bins_same_cg, num_per_dmr, prop, has_dmr) %>%
    group_by(num_cg) %>%
    .[order(.$num_cg),] %>%
    mutate(prop = ifelse(has_dmr == FALSE, 0, prop)) %>%
    mutate(num = n()) %>%
    mutate(has_dmr = ifelse(num == 1, TRUE, has_dmr)) %>% filter(has_dmr == TRUE) %>%
    ggplot(aes(x = num_cg, y = prop)) +
      geom_point(alpha = 0.4) +
      geom_smooth()

  list(genome_bins, plot)
}


sourceGenes <- function(counts, gene_source = c("biomaRt", "ExperimentHub"), gene_feature=c("gene", "promoter")) {
  gene_source <- match.arg(tolower(gene_source), c("biomart", "experimenthub"))
  gene_feature <- match.arg(tolower(gene_feature), c("gene", "promoter"))

  if(gene_source == "biomart") {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "GRCh37")
    genes <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id",
                     "chromosome_name", "start_position", "end_position", "strand",
                     "transcription_start_site" , "transcript_start", "transcript_end", "transcript_length"),
      filters = "chromosome_name",
      values = gsub("chr", "", unique(counts$seqnames)[1:22]),
      mart = ensembl
    )
    genes <- genes %>% .[order(.$chromosome_name, .$start_position), ]
    colnames(genes) <- c("ensembl_gene_id", "gene_name", "ensembl_transcript_id", "seqnames", "start", "end", "strand", "TSS", "t_start", "t_end", "t_length")
    genes <- genes %>% group_by(ensembl_gene_id) %>% filter(t_length == max(t_length)) %>% ungroup()
    genes$seqnames <- paste0("chr", as.character(genes$seqnames))
    genes <- genes %>% group_by(ensembl_gene_id) %>% mutate(num_t = 1:dplyr::n()) %>% filter(num_t == 1) %>% ungroup()
  } else if (gene_source == "experimenthub") {
    eh <- ExperimentHub()
    genes <- eh[["EH3132"]]
  }
  genes <- genes %>% data.frame() %>% ungroup()
  # ensure ensembl_id name consistency
  colnames(genes)[grepl("ENSG", genes[1,])] <- "ensembl_gene_id"
  # double check only 1 ensembl id per gene entry
  genes$width <- genes$end - genes$start + 1
  genes <- genes %>% group_by(ensembl_gene_id) %>% filter(width == max(width)) %>% ungroup()
  # add random identifier
  genes$gene_num <- 1:nrow(genes)
  #is it promoter???
  if(length(gene_feature) == 1 & gene_feature == "promoter") {
    genes$start <- ifelse(genes$strand == 1, genes$start - 2000, genes$end - 500)
    genes$end <- ifelse(genes$strand == 1, genes$start + 2500, genes$end + 2000)
    genes <- genes[,!colnames(genes) == "width"]
  }
  # add number of cpgs overlapping gene
  gene_cpg <- find_overlaps(as_granges(genes), as_granges(counts)) %>% data.frame()
  gene_cpg <- gene_cpg %>% group_by(gene_num) %>% summarise(n_cpg = dplyr::n()) %>% ungroup() %>% data.frame()

  genes$num_cg <- gene_cpg[match(genes$gene_num, gene_cpg$gene_num), "n_cpg"]
  genes[is.na(genes)] <- 0
  genes
}

annoGeneDmr <- function(counts, dmrs, gene_source = c("biomart", "experimenthub", "self"), genes) { #why do I have the gene source here?? why???
  gene_source <- match.arg(tolower(gene_source), c("biomart", "experimenthub", "self"))
  genes <- data.frame(genes)
  # ensure ensembl_id name consistency
  colnames(genes)[grepl("ENSG", genes[1,])] <- "ensembl_gene_id"
  # double check only 1 ensembl id per gene entry
  genes$width <- genes$end - genes$start + 1
  genes <- genes %>% group_by(ensembl_gene_id) %>% filter(width == max(width)) %>% ungroup()
  # add random identifier
  genes$gene_num <- 1:nrow(genes)
  # add number of cpgs overlapping gene
  gene_cpg <- find_overlaps(as_granges(genes), as_granges(counts)) %>% data.frame()
  gene_cpg <- gene_cpg %>% group_by(gene_num) %>% summarise(n_cpg = n()) %>% ungroup() %>% data.frame()

  genes$num_cg <- gene_cpg[match(genes$gene_num, gene_cpg$gene_num), "n_cpg"]
  genes[is.na(genes)] <- 0
  if(!"TSS" %in% colnames(genes)) {
    genes$TSS <- ifelse(genes$strand == "+", genes$start, genes$end) # above I have strand for 1/-1
  }
  # assume dmrs are already in order of most to least signif
  dmrs$rank <- 1:nrow(dmrs)
  # add if gene overlaps dmr
  overlap_dmr <- find_overlaps(as_granges(dmrs), as_granges(genes)) %>% data.frame() %>%
    mutate(combo = paste0(ensembl_gene_id, "-DMR", rank))
  anno_dmr <- pair_nearest(as_granges(dmrs), as_granges(genes)) %>% data.frame() %>%
    mutate(combo = paste0(ensembl_gene_id, "-DMR", rank))
  #dist gene to dmr
  anno_dmr <- anno_dmr %>%
    mutate(tss_start = TSS - granges.x.start, tss_end = TSS - granges.x.end) %>%
    mutate(tss_rel_peak = ifelse(combo %in% overlap_dmr$combo, "Overlap", NA),
            tss_rel_peak = case_when(is.na(tss_rel_peak) & tss_start > 0 & tss_end > 0 & tss_end < tss_start ~ "Peak_upstream",
                                     is.na(tss_rel_peak) & tss_start < 0 & tss_end < 0 &  tss_end < tss_start ~ "Peak_downstream", TRUE ~ tss_rel_peak),
           tss_rel_peak = case_when(tss_rel_peak == "Peak_upstream" & granges.y.strand == "-" ~ "Peak_downstream",
                                   tss_rel_peak == "Peak_downstream" & granges.y.strand == "-" ~ "Peak_upstream", TRUE ~ tss_rel_peak),
          min_dist = case_when(tss_rel_peak == "Overlap" ~ 0, tss_rel_peak == "Peak_downstream" & granges.y.strand == "+" ~ tss_start,
                                tss_rel_peak == "Peak_downstream" & granges.y.strand == "-" ~ tss_end, tss_rel_peak == "Peak_upstream" & granges.y.strand == "+" ~ tss_end,
                                tss_rel_peak == "Peak_upstream" & granges.y.strand == "-" ~ tss_start))
  anno_dmr

}

plotBiasGrouped <- function(genes, anno_genes, regression_line=FALSE, log2_scale=FALSE) {
  genes$width <- genes$end - genes$start + 1
  genes$has_dmr <- ifelse(genes$ensembl_gene_id %in% filter(anno_genes, abs(min_dist) == 0)$ensembl_gene_id, TRUE, FALSE)
  genes <- genes[order(.$num_cg),]
  genes$bin_group <- rep(1:ceiling(nrow(genes)/100), each = 100)[1:nrow(genes)]
  genes <- genes %>%
    group_by(bin_group) %>%
    mutate(num_cg_bin_group = mean(num_cg)) %>%
    group_by(num_cg_bin_group) %>%
    mutate(num_bins_group_same_cg = dplyr::n()) %>%
    group_by(num_cg_bin_group, has_dmr) %>%
    mutate(num_per_dmr = n(), prop = dplyr::n()/num_bins_group_same_cg) %>%
    distinct(num_cg_bin_group, num_bins_group_same_cg, num_per_dmr, prop, has_dmr) %>%
    group_by(num_cg_bin_group) %>%
    .[order(.$num_cg_bin_group),] %>%
    mutate(prop = ifelse(has_dmr == FALSE, 0, prop)) %>%
    mutate(num = dplyr::n()) %>%
    mutate(has_dmr = ifelse(num == 1, TRUE, has_dmr)) %>% filter(has_dmr == TRUE)
  plot <- genes %>%
    ggplot(aes(x = num_cg_bin_group, y = prop)) +
    geom_point(alpha = 0.4)
  if(log2_scale == TRUE) {
    plot <- genes %>%
      ggplot(aes(x = log2(num_cg_bin_group), y = prop)) +
      geom_point(alpha = 0.4)
  }
  if(regression_line == TRUE) {
    plot <- plot + geom_smooth()
  }
  plot
}

plotBiasGroupedMedian <- function(genes, anno_genes, regression_line=FALSE, log2_scale=FALSE) {
  genes <- genes %>% mutate(width = end - start + 1, has_dmr = ifelse(ensembl_gene_id %in% filter(anno_genes, abs(min_dist) == 0)$ensembl_gene_id, TRUE, FALSE)) %>%
    group_by(num_cg) %>%
    mutate(num_bins_same_cg = dplyr::n()) %>%
    .[order(.$num_cg),] %>%
    ungroup() %>%
    mutate(bin_group = rep(1:ceiling(nrow(genes)/100), each = 100)[1:nrow(genes)]) %>%
    group_by(bin_group) %>%
    mutate(num_cg_bin_group = median(num_cg)) %>%
    group_by(num_cg_bin_group) %>%
    mutate(num_bins_group_same_cg = dplyr::n()) %>%
    group_by(num_cg_bin_group, has_dmr) %>%
    mutate(num_per_dmr = dplyr::n(), prop = dplyr::n()/num_bins_group_same_cg) %>%
    distinct(num_cg_bin_group, num_bins_group_same_cg, num_per_dmr, prop, has_dmr) %>%
    group_by(num_cg_bin_group) %>%
    .[order(.$num_cg_bin_group),] %>%
    mutate(prop = ifelse(has_dmr == FALSE, 0, prop)) %>%
    mutate(num = dplyr::n()) %>%
    mutate(has_dmr = ifelse(num == 1, TRUE, has_dmr)) %>% filter(has_dmr == TRUE)
  plot <- genes %>%
    ggplot(aes(x = num_cg_bin_group, y = prop)) +
    geom_point(alpha = 0.4)
  if(log2_scale == TRUE) {
    plot <- genes %>%
      ggplot(aes(x = log2(num_cg_bin_group), y = prop)) +
      geom_point(alpha = 0.4)
  }
  if(regression_line == TRUE) {
    plot <- plot + geom_smooth()
  }
  plot
}

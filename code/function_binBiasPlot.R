library(plyranges)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ChIPseqSpikeInFree)



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

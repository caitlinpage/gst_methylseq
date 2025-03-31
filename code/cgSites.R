
#' Obtain CG genomic positions for given genome
#'
#' @param genome BSgenome object
#' @param x First seqnames index to return
#' @param y Last seqnames index to return. Is designed to prevent the inclusion of scaffolds, which are not part of the primary mapping.
#'
#' @returns data.frame object with the CG positions of the genome. Is compatible with GRanges as chromosome number is labelled as seqnames.
#' @export
#'
#' @examples
get_cg_sites <- function(genome, x, y) {
  seq_names <- GenomeInfoDb::seqnames(genome)[x:y]
  anno_seq <- lapply(seq_names, function(x) {
    cbind(data.frame(Biostrings::matchPattern("CG", genome[[x]])),
          seqnames = x
    )
  }) %>%
    dplyr::bind_rows()
  anno_seq <- anno_seq %>% dplyr::mutate(pos = paste0(seqnames, "-", start)) %>%
    dplyr::relocate(pos, seqnames)
  anno_seq
}

library(BSgenome.Hsapiens.UCSC.hg19)


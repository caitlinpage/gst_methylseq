# Analysis: 04_compareGeneOntology.Rmd
library(org.Hs.eg.db)
library(dplyr)
x <- org.Hs.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
entrez_map <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[entrez_map])

xx <- as.list(org.Hs.egGO2ALLEGS)
if(length(xx) > 0){
  # Gets the Entrez Gene identifiers for the top 2nd and 3nd GO identifiers
  goids <- xx[2:3]
  # Gets all the Entrez Gene identifiers for the first element of goids
  goids[[1]]
  # Evidence code for the mappings
  names(goids[[1]])
}

go_entrez <- data.frame(cbind(go = gsub("\\..*", "", names(unlist(xx))), entrez = unlist(xx)))
go_entrez

go_entrez <- go_entrez %>% mutate(num_cg = data.frame(test_out_num$freq)[match(.$entrez, data.frame(test_out_num$freq)$freq_make_table), "Freq"],
                                  num_cg_bin = data.frame(test_out$freq)[match(.$entrez, data.frame(test_out$freq)$freq_make_table), "Freq"])

go_entrez <- go_entrez %>% mutate(array_cat_signif = ifelse(go %in% rownames(sig_terms_array), TRUE, FALSE),
                                  array_cat_signif_bias = ifelse(go %in% rownames(sig_terms_array_biased), TRUE, FALSE),
                                  seq_cat_signif = ifelse(go %in% GO.wall$category, TRUE, FALSE),
                                  seq_cat_signif_bias = ifelse(go %in% GO.wall.null$category, TRUE, FALSE),
                                  mm_seq = ifelse(go %in% mm_hack$GOID, TRUE, FALSE),
                                  mm_seq_null = ifelse(go %in% mm_hack_null$GOID, TRUE, FALSE),
                                  mm_fn_bin = ifelse(go %in% rownames(test_b), TRUE, FALSE),
                                  mm_fn_bin_null = ifelse(go %in% test_b_null$GOID, TRUE, FALSE),
                                  mm_fn_num = ifelse(go %in% test_b_num$GOID, TRUE, FALSE),
                                  mm_fn_num_null = ifelse(go %in% test_b_num_null$GOID, TRUE, FALSE))
go_entrez <- go_entrez %>% mutate(mm_fn_num_filt_dmr = ifelse(go %in% filter(test_b_num_dmr, FDR < 0.01)$GOID, TRUE, FALSE))

go_entrez <- go_entrez %>% mutate(mm_fn_num_filt_dmr_null = ifelse(go %in% filter(test_b_num_dmr_null, FDR < 0.01)$GOID, TRUE, FALSE))

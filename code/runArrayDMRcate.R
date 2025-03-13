
## dmrcate array
library(DMRcate)
library(dplyr)

file_source <- c("natural_killer", "b_cell")
files_new <- character()
for(i in 1:2) {
  files <- list.files(paste0("../data/blood/blood_microarray/", file_source[i])) %>% gsub(".idat.gz", "", .)
  files <- files[!grepl("Red", files)]
  files <- files %>% gsub("_Grn", "", .)
  files <- paste0("../data/blood/blood_microarray/", file_source[i], "/", files)
  files_new <- c(files_new, files)
}

rgSet <- read.metharray(files_new)
rgSet

library(ExperimentHub)

eh <- ExperimentHub()
FlowSorted.Blood.EPIC <- eh[["EH1136"]]
tcell <- FlowSorted.Blood.EPIC[,colData(FlowSorted.Blood.EPIC)$CD4T==100 |
                                 colData(FlowSorted.Blood.EPIC)$CD8T==100]


blood_samples <- cbind(sample_num = 1:12) %>% data.frame()
blood_samples <- blood_samples %>% mutate(sample_group = rep(file_source, times = 6) %>% .[order(., decreasing = TRUE)])
blood_samples <- blood_samples %>% mutate(sample_id = substring(files_new, first = stringr::str_locate(files_new, "G")[,1])
)
blood_samples <- blood_samples %>% mutate(short_group = ifelse(sample_group == "natural_killer", "NK", "B"))
blood_samples <- blood_samples %>% mutate(id_group = paste0(short_group, "_", substring(sample_id, first = 1, last = stringr::str_locate(sample_id, "_")[,1] - 1)))

sampleNames(rgSet) <- blood_samples$id_group

detP <- detectionP(rgSet)

remove <- apply(detP, 1, function (x) any(x > 0.01))
rgSet <- preprocessFunnorm(rgSet)

rgSet <- rgSet[!seqnames(rgSet) %in% c("chrX", "chrY"),]
rgSet <- rgSet[!rownames(rgSet) %in% names(which(remove)),]

array_m <- getM(rgSet)

noSNPs <- rmSNPandCH(array_m, dist=2, mafcut=0.05)

rgSet <- rgSet[rownames(noSNPs),]
colnames(noSNPs) <- colnames(rgSet)
assays(rgSet)[["M"]] <- noSNPs
assays(rgSet)[["Beta"]] <- ilogit2(noSNPs)

type <- factor(blood_samples$short_group)
design <- model.matrix(~type)
myannotation <- cpg.annotate("array", rgSet, arraytype = "EPICv1",
                             analysis.type="differential", design=design, coef=2)
myannotation

dmrcate_array_res <- dmrcate(myannotation, lambda=1000, C=2)

dmrcate_array_res <- extractRanges(dmrcate_array_res, genome = "hg19")

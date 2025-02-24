# copied from MissMethyl

.getGO <- function(){
  if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("org.Hs.eg.db package required but not installed.")
  egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
  GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
  d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
  GeneID.PathID <- GeneID.PathID[d, ]
  GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,
                                                      keys=unique(GeneID.PathID$go_id),
                                                      columns=c("GOID","ONTOLOGY","TERM"),
                                                      keytype="GOID"))
  go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)

  list(idList=go, idTable=GOID.TERM)
}


.estimatePWF <- function(D,bias)
  # An alternative to goseq function nullp, which is transformation invariant
  # Belinda Phipson and Gordon Smyth
  # 6 March 2015
{
  prior.prob <- bias
  o <- order(bias)
  prior.prob[o] <- limma::tricubeMovingAverage(D[o],span=0.5)
  prior.prob
}
test_pwf <- .estimatePWF(D=test.de,bias=as.vector(freq_genes))

.plotBias <- function(D,bias)
  # Plotting function to show gene level CpG density bias
  # Belinda Phipson
  # 5 March 2015
{
  o <- order(bias)
  splitf <- rep(1:100,each=200)[1:length(bias)]
  avgbias <- tapply(bias[o],factor(splitf),mean)
  sumDM <- tapply(D[o],factor(splitf),sum)
  propDM <- sumDM/table(splitf)
  graphics::par(mar=c(5,5,2,2))
  graphics::plot(avgbias,as.vector(propDM),
                 xlab="Number of CpGs per gene (binned)",
                 ylab="Proportion Differential Methylation",cex.lab=1.5,
                 cex.axis=1.2)
  graphics::lines(stats::lowess(avgbias,propDM),col=4,lwd=2)
}

# get the genes and dmrs organised
run_miss_methyl_a <- function(biomart_genes, test_anno) { # add something about extending the bias to include upstream of gene maybe??
  #get entrez gene ids and gene symbols
  entrez_genes <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys=AnnotationDbi::keys(org.Hs.eg.db),
                                        columns=c("ENTREZID","SYMBOL"),
                                        keytype="ENTREZID")
  biomart_genes$entrez_id <- entrez_genes[match(biomart_genes$gene_name, entrez_genes$SYMBOL), "ENTREZID"]
  biomart_genes_filt <- biomart_genes[!biomart_genes$entrez_id %in% NA,]
  biomart_genes_filt <- biomart_genes_filt %>% group_by(gene_name) %>% filter(num_cg == max(num_cg)) %>% ungroup() #this line will be slow

  test_anno$entrez_id <- entrez_genes[match(test_anno$gene_name, entrez_genes$SYMBOL), "ENTREZID"]
  test_anno_filt <- test_anno[!test_anno$entrez_id %in% NA,]
  test_anno_filt <- test_anno_filt[test_anno_filt$tss_rel_peak %in% "Overlap",]
  biomart_genes_filt$bias <- biomart_genes_filt$num_cg


  # just entrez genes and number of overlapping CpG
  test_uni_freq <- biomart_genes_filt %>% .[order(.$entrez_id),] %>% group_by(entrez_id) %>% distinct(entrez_id, bias) %>% ungroup() %>% data.frame() # probs slow
  test_uni_freq <- test_uni_freq[test_uni_freq$bias > 0,]
  test_uni_freq

  freq_make_table <- rep(test_uni_freq$entrez_id, test_uni_freq$bias) #this has replaced the loop and sped it up

  freq_make_table <- table(freq_make_table)
  biomart_genes_filt <- biomart_genes_filt[biomart_genes_filt$bias > 0,]
  test_anno_filt <- test_anno_filt[test_anno_filt$ensembl_gene_id %in% biomart_genes_filt$ensembl_gene_id,]
  biomart_genes_filt$is_de <- ifelse(biomart_genes_filt$entrez_id %in% test_anno_filt$entrez_id, 1, 0)
  biomart_genes_filt$is_de <- as.integer(biomart_genes_filt$is_de)
  out <- list(
    #fract.counts = biomart_genes_filt %>%
    # mutate(sigid = entrez_id) %>% distinct(sigid) %>% .[order(.$sigid),] %>% data.frame() %>%
    # mutate(frac = as.double(1)),
    sig.eg = test_anno_filt %>% distinct(entrez_id) %>% .[order(.$entrez_id),],
    universe = biomart_genes_filt %>% distinct(entrez_id) %>%
      .[order(.$entrez_id),] %>% .$entrez_id,
    freq = freq_make_table,
    equiv = as.array(freq_make_table), # column names wrong but maybe not too bad?
    de = biomart_genes_filt[order(biomart_genes_filt$entrez_id),] %>% .$is_de
  )
  out
}

# does all the testing
run_miss_methyl_b <- function(out, method, plot_bias=FALSE, ...) {
  sorted.eg.sig <- out$sig.eg #entrez ids for sig cpgs - genes that have sig cpgs
  eg.universe <- out$universe # entrez ids all genes
  freq_genes <- out$freq # num cpgs per gene
  test.de <- out$de #gene de or not
  #frac <- out$fract.counts
  equiv <- out$equiv
  ## setting up stuff

  # various go filters
  # Check collection is a list with character vectors
  go <- .getGO()
  collection <- go$idList
  #only if it's not a list
  #collection <- list(collection=go$idList)
  collection <- lapply(collection, as.character)
  # Make sure gene set collections don't have any NAs
  collection <- lapply(collection, function(x) x[!is.na(x)])
  # Remove genes that are NOT in the universe from collections
  collection <- lapply(collection, function(x) x[x %in% eg.universe])
  # Remove collections with no genes left after universe filter
  inUniv <- sapply(collection, function(x) length(x) > 0)
  collection <- collection[inUniv]

  # set up testing
  results <- matrix(NA,ncol=4,nrow=length(collection))
  colnames(results) <- c("N", "DE", "P.DE","FDR")
  rownames(results) <- names(collection)
  results[,"N"] <- unlist(lapply(collection,length))
  SigGenesInSet <- rep(NA,length(collection)) # this one is optional

  results[,"DE"] <- unlist(lapply(collection, function(x)
    sum((sorted.eg.sig %in% x))))
  Nuniverse <- length(eg.universe)
  m <- length(sorted.eg.sig)

  # speed up testing by removing categories with 0 de genes (won't be signif)
  # so don't test, and then add them back at the end
  results <- data.frame(results)
  results_null <- results[results$DE == 0,]
  results_null$P.DE <- 1
  results_null$FDR <- 1

  results <- results[results$DE != 0,]
  collection <- collection[rownames(results)]
  results <- as.matrix(results)

  # testing

  if(method != "None") {
    test_pwf <- .estimatePWF(D=test.de,bias=as.vector(freq_genes))
    plot <- .plotBias(D = test.de, bias = as.vector(freq_genes))
  }
  if(method == "Wallenius") {
    results[, c("P.DE")] <-
      t(sapply(1:length(collection), function(i) {
        InSet <- eg.universe %in% collection[[i]] # genes in category
        pw.red <- sum(test_pwf[InSet], na.rm = TRUE) / results[i, "N"] # so i don't know why I get na's - or do i?
        pw.white <- sum(test_pwf[!InSet], na.rm = TRUE) / (Nuniverse - results[i, "N"]) #anyway the na removal is not in original code
        odds <- pw.red / pw.white
        p_val <- BiasedUrn::pWNCHypergeo(results[i,"DE"],
                                         results[i,"N"],
                                         Nuniverse-results[i,"N"],
                                         m,odds,lower.tail=FALSE) +
          BiasedUrn::dWNCHypergeo(results[i,"DE"],
                                  results[i,"N"],
                                  Nuniverse-results[i,"N"],
                                  m,odds)
        #the optional thing
        # Get gene symbols of significant genes
        SigGenesEntrezID <- sorted.eg.sig[sorted.eg.sig %in% collection[[i]]]
        SigGenesSymbol <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
                                                                 keys = SigGenesEntrezID,
                                                                 columns = "SYMBOL"))
        SigGenesInSet[i] <- paste(SigGenesSymbol$SYMBOL,collapse=",")
        if(p_val==0){
          message("Pvalue of exactly zero detected. Performing hypergeometric
                test for gene set ", rownames(results)[i])
          p_val <- stats::phyper(q=results[i,"DE"]-0.5,m=m,
                                 n=Nuniverse-m,k=results[i,"N"],
                                 lower.tail=FALSE)
        }
        #  }
        c(p_val)  # return the values as a vector
      })) #chatgpt translated loop to apply
  } else if(method == "None") {
    results[, c("P.DE")] <-
      t(sapply(1:length(collection), function(i) {
        p_val <- stats::phyper(q=results[i,"DE"]-0.5,m=m,
                               n=Nuniverse-m,k=results[i,"N"],
                               lower.tail=FALSE)
        #the optional thing
        # Get gene symbols of significant genes
        SigGenesEntrezID <- sorted.eg.sig[sorted.eg.sig %in% collection[[i]]]
        SigGenesSymbol <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
                                                                 keys = SigGenesEntrezID,
                                                                 columns = "SYMBOL"))
        SigGenesInSet[i] <- paste(SigGenesSymbol$SYMBOL,collapse=",")
        c(p_val)
      }))
  } else if(method == "Fishers") {
    results[, c("P.DE")] <-
      t(sapply(1:length(collection), function(i) {
        InSet <- eg.universe %in% collection[[i]]
        pw.red <- sum(test_pwf[InSet], na.rm = TRUE) / results[i, "N"] # so i don't know why I get na's - or do i?
        pw.white <- sum(test_pwf[!InSet], na.rm = TRUE) / (Nuniverse - results[i, "N"]) #anyway the na removal is not in original code
        odds <- pw.red / pw.white
        p_val <- BiasedUrn::pFNCHypergeo(results[i,"DE"],
                                         results[i,"N"],
                                         Nuniverse-results[i,"N"],
                                         m,odds,lower.tail=FALSE) +
          BiasedUrn::dFNCHypergeo(results[i,"DE"],
                                  results[i,"N"],
                                  Nuniverse-results[i,"N"],
                                  m,odds)
        #the optional thing
        # Get gene symbols of significant genes
        SigGenesEntrezID <- sorted.eg.sig[sorted.eg.sig %in% collection[[i]]]
        SigGenesSymbol <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
                                                                 keys = SigGenesEntrezID,
                                                                 columns = "SYMBOL"))
        SigGenesInSet[i] <- paste(SigGenesSymbol$SYMBOL,collapse=",")
        if(p_val==0){
          message("Pvalue of exactly zero detected. Performing hypergeometric
                test for gene set ", rownames(results)[i])
          p_val <- stats::phyper(q=results[i,"DE"]-0.5,m=m,
                                 n=Nuniverse-m,k=results[i,"N"],
                                 lower.tail=FALSE)
        }
        c(p_val)  # return the values as a vector
      }))
  }


  results <- results %>% data.frame()
  results[,"FDR"] <- stats::p.adjust(results[,"P.DE"],method="BH")
  results[,"DE"] <- floor(results[,"DE"])
  results <- rbind(results, results_null)
  results <- data.frame(results, SigGenesInSet)

  results <- merge(go$idTable,results,by.x="GOID",by.y="row.names")
  rownames(results) <- results$GOID
  results <- results %>% .[order(.$FDR),]
  list(results, plot)

}

# main wrapper function
run_gst_seq <- function(dmrs, tested_cpgs,
                        gene_source = c("biomaRt", "ExperimentHub"), gene_feature=c("gene", "promoter"),
                        method = c("Wallenius", "Fishers", "None"), plot_bias=TRUE, ...) {
  gene_source <- match.arg(tolower(gene_source), c("biomart", "experimenthub"))
  gene_feature <- match.arg(tolower(gene_feature), c("gene", "promoter"))
  genes <- sourceGenes(counts = tested_cpgs, gene_source = gene_source, gene_feature = gene_feature)
  if(nrow(genes) == 0) {
    stop("sourceGenes failed")
  }
  if(nrow(dmrs) > 30000) {
    stop("Please appropriately filter DMRs to most significant") #idk
  }
  anno <- annoGeneDmr(counts=tested_cpgs, dmrs = dmrs, gene_source = gene_source, genes = genes)
  if(nrow(anno) == 0) {
    stop("anno failed")
  }
  out <- run_miss_methyl_a(genes, anno)
  results <- run_miss_methyl_b(out, correct_bias=correct_bias, plot_bias=plot_bias, ...)

  list(res = results, genes = genes, anno = anno, out = out)
}

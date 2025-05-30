reversemapping=function(map){
  tmp=unlist(map,use.names=FALSE)
  names(tmp)=rep(names(map),times=as.numeric(summary(map)[,1]))
  return(split(names(tmp),as.vector(tmp)))
}

goseq_mod <- function(pwf,genome,id,gene2cat=NULL,test.cats=c("GO:CC","GO:BP","GO:MF"),method="Wallenius",repcnt=2000,use_genes_without_cat=FALSE){
  ################# Input pre-processing and validation ###################
  #Do some validation of input variables
  if(any(!test.cats%in%c("GO:CC","GO:BP","GO:MF","KEGG"))){
    stop("Invalid category specified.  Valid categories are GO:CC, GO:BP, GO:MF or KEGG")
  }
  if((missing(genome) | missing(id))){
    if(is.null(gene2cat)){
      stop("You must specify the genome and gene ID format when automatically fetching gene to GO category mappings.")
    }
    #If we're using user specified mappings, this obviously isn't a problem
    genome='dummy'
    id='dummy'
  }
  if(!any(method%in%c("Wallenius","Sampling","Hypergeometric", "Fishers"))){
    stop("Invalid calculation method selected.  Valid options are Wallenius, Sampling & Hypergeometric.")
  }
  if(!is.null(gene2cat) && (!is.data.frame(gene2cat) & !is.list(gene2cat))){
    stop("Was expecting a dataframe or a list mapping categories to genes.  Check gene2cat input and try again.")
  }

  #Factors are evil
  pwf=unfactor(pwf)
  gene2cat=unfactor(gene2cat)

  ###################### Data fetching and processing ########################
  if(is.null(gene2cat)){
    #When we fetch the data using getgo it will be in the list format
    message("Fetching GO annotations...")
    gene2cat=getgo(rownames(pwf),genome,id,fetch.cats=test.cats)
    names(gene2cat)=rownames(pwf)
    #Do the two rebuilds to remove any nulls
    cat2gene=reversemapping(gene2cat)
    gene2cat=reversemapping(cat2gene)
  }else{
    #The gene2cat input accepts a number of formats, we need to check each of them in term
    message("Using manually entered categories.")
    #The options are a flat mapping (that is a data frame or matrix) or a list, where the list can be either gene->categories or category->genes
    if(class(gene2cat)!="list"){
      #it's not a list so it must be a data.frame, work out which column contains the genes
      genecol_sum=as.numeric(apply(gene2cat,2,function(u){sum(u%in%rownames(pwf))}))
      genecol=which(genecol_sum!=0)
      if(length(genecol)>1){
        genecol=genecol[order(-genecol_sum)[1]]
        warning(paste("More than one possible gene column found in gene2cat, using the one headed",colnames(gene2cat)[genecol]))
      }
      if(length(genecol)==0){
        genecol=1
        warning(paste("Gene column could not be identified in gene2cat conclusively, using the one headed",colnames(gene2cat)[genecol]))
      }
      othercol=1
      if(genecol==1){othercol=2}
      #Now put it into our delicious listy format
      gene2cat=split(gene2cat[,othercol],gene2cat[,genecol])
      #Do the appropriate builds
      cat2gene=reversemapping(gene2cat)
      gene2cat=reversemapping(cat2gene)
    }
    #!!!!
    #The following conditional has been flagged as a potential issue when using certain
    #types of input where the category names are the same as gene names (which seems like
    #something you should avoid anyway...).  Leave it for now
    #!!!!
    #We're now garunteed to have a list (unless the user screwed up the input) but it could
    #be category->genes rather than the gene->categories that we want.
    if(sum(unique(unlist(gene2cat,use.names=FALSE))%in%rownames(pwf))>sum(unique(names(gene2cat))%in%rownames(pwf))){
      gene2cat=reversemapping(gene2cat)
    }
    #Alright, we're garunteed a list going in the direction we want now.  Throw out genes which we will not use
    gene2cat=gene2cat[names(gene2cat)%in%rownames(pwf)]

    #Rebuild because it's a fun thing to do
    cat2gene=reversemapping(gene2cat)
    gene2cat=reversemapping(cat2gene)

    ## make sure we remove duplicate entries .. e.g. see
    ## http://permalink.gmane.org/gmane.science.biology.informatics.conductor/46876
    cat2gene=lapply(cat2gene,function(x){unique(x)})
    gene2cat=lapply(gene2cat,function(x){unique(x)})
  }

  nafrac=(sum(is.na(pwf$pwf))/nrow(pwf))*100
  if(nafrac>50){
    warning(paste("Missing length data for ",round(nafrac),"% of genes.  Accuarcy of GO test will be reduced.",sep=''))
  }
  #Give the genes with unknown length the weight used by the median gene (not the median weighting!)
  pwf$pwf[is.na(pwf$pwf)]=pwf$pwf[match(sort(pwf$bias.data[!is.na(pwf$bias.data)])[ceiling(sum(!is.na(pwf$bias.data))/2)],pwf$bias.data)]

  ###################### Calculating the p-values ########################
  # Remove all the genes with unknown GOterms
  unknown_go_terms=nrow(pwf)-length(gene2cat)
  if((!use_genes_without_cat) && unknown_go_terms>0 ){
    message(paste("For",unknown_go_terms,"genes, we could not find any categories. These genes will be excluded."))
    message("To force their use, please run with use_genes_without_cat=TRUE (see documentation).")
    message("This was the default behavior for version 1.15.1 and earlier.")
    pwf=pwf[rownames(pwf) %in% names(gene2cat),]
  }
  #A few variables are always useful so calculate them
  cats=names(cat2gene)
  DE=rownames(pwf)[pwf$DEgenes==1]
  num_de=length(DE)
  num_genes=nrow(pwf)
  pvals=data.frame(category=cats,over_represented_pvalue=NA,under_represented_pvalue=NA,stringsAsFactors=FALSE,numDEInCat=NA,numInCat=NA)
  if(method=="Sampling"){
    #We need to know the number of DE genes in each category, make this as a mask that we can use later...
    num_DE_mask=rep(0,length(cats))
    a=table(unlist(gene2cat[DE],FALSE,FALSE))
    num_DE_mask[match(names(a),cats)]=as.numeric(a)
    num_DE_mask=as.integer(num_DE_mask)
    #We have to ensure that genes not associated with a category are included in the simulation, to do this they need an empty entry in the gene2cat list
    gene2cat=gene2cat[rownames(pwf)]
    names(gene2cat)=rownames(pwf)
    message("Running the simulation...")
    #Now do the actual simulating
    lookup=matrix(0,nrow=repcnt,ncol=length(cats))
    for(i in 1:repcnt){
      #A more efficient way of doing weighted random sampling without replacment than the built in function
      #The order(runif...)[1:n] bit picks n genes at random, weighting them by the PWF
      #The table(as.character(unlist(...))) bit then counts the number of times this random set occured in each category
      a=table(as.character(unlist(gene2cat[order(runif(num_genes)^(1/pwf$pwf),decreasing=TRUE)[1:num_de]],FALSE,FALSE)))
      lookup[i,match(names(a),cats)]=a
      pp(repcnt)
    }
    message("Calculating the p-values...")
    #The only advantage of the loop is it uses less memory...
    #for(i in 1:length(cats)){
    #	pvals[i,2:3]=c((sum(lookup[,i]>=num_DE_mask[i])+1)/(repcnt+1),(sum(lookup[,i]<=num_DE_mask[i])+1)/(repcnt+1))
    #	pp(length(cats))
    #}
    pvals[,2]=(colSums(lookup>=outer(rep(1,repcnt),num_DE_mask))+1)/(repcnt+1)
    pvals[,3]=(colSums(lookup<=outer(rep(1,repcnt),num_DE_mask))+1)/(repcnt+1)
  }
  if(method=="Wallenius"){
    message("Calculating the p-values...")
    #All these things are just to make stuff run faster, mostly because comparison of integers is faster than string comparison
    degenesnum=which(pwf$DEgenes==1)
    #Turn all genes into a reference to the pwf object
    cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
    #This value is used in every calculation, by storing it we need only calculate it once
    alpha=sum(pwf$pwf)
    #Each category will have a different weighting so needs its own test
    pvals[,2:3]=t(sapply(cat2genenum,function(u){
      #The number of DE genes in this category
      num_de_incat=sum(degenesnum%in%u)
      #The total number of genes in this category
      num_incat=length(u)
      #This is just a quick way of calculating weight=avg(PWF within category)/avg(PWF outside of category)
      avg_weight=mean(pwf$pwf[u])
      weight=(avg_weight*(num_genes-num_incat))/(alpha-num_incat*avg_weight)
      if(num_incat==num_genes){ weight=1 } #case for the root GO terms
      #Now calculate the sum of the tails of the Wallenius distribution (the p-values)
      c(dWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight)
        +pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight,lower.tail=FALSE),
        pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight))
    }))
  }
  if(method=="Fishers"){
    message("Calculating the p-values...")
    #All these things are just to make stuff run faster, mostly because comparison of integers is faster than string comparison
    degenesnum=which(pwf$DEgenes==1)
    #Turn all genes into a reference to the pwf object
    cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
    #This value is used in every calculation, by storing it we need only calculate it once
    alpha=sum(pwf$pwf)
    #Each category will have a different weighting so needs its own test
    pvals[,2:3]=t(sapply(cat2genenum,function(u){
      #The number of DE genes in this category
      num_de_incat=sum(degenesnum%in%u)
      #The total number of genes in this category
      num_incat=length(u)
      #This is just a quick way of calculating weight=avg(PWF within category)/avg(PWF outside of category)
      avg_weight=mean(pwf$pwf[u])
      weight=(avg_weight*(num_genes-num_incat))/(alpha-num_incat*avg_weight)
      if(num_incat==num_genes){ weight=1 } #case for the root GO terms
      #Now calculate the sum of the tails of the Wallenius distribution (the p-values)
      c(dFNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight)
        +pFNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight,lower.tail=FALSE),
        pFNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight))
    }))
  }
  if(method=="Hypergeometric"){
    message("Calculating the p-values...")
    #All these things are just to make stuff run faster, mostly because comparison of integers is faster than string comparison
    degenesnum=which(pwf$DEgenes==1)
    #Turn all genes into a reference to the pwf object
    cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
    #Simple hypergeometric test, one category at a time
    pvals[,2:3]=t(sapply(cat2genenum,function(u){
      #The number of DE genes in this category
      num_de_incat=sum(degenesnum%in%u)
      #The total number of genes in this category
      num_incat=length(u)
      #Calculate the sum of the tails of the hypergeometric distribution (the p-values)
      c(dhyper(num_de_incat,num_incat,num_genes-num_incat,num_de)+phyper(num_de_incat,num_incat,num_genes-num_incat,num_de,lower.tail=FALSE),phyper(num_de_incat,num_incat,num_genes-num_incat,num_de))
    }))
  }
  #Populate the count columns...
  degenesnum=which(pwf$DEgenes==1)
  cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
  pvals[,4:5]=t(sapply(cat2genenum,function(u){
    c(sum(degenesnum%in%u),length(u))
  }))

  #Finally, sort by p-value
  pvals=pvals[order(pvals$over_represented_pvalue),]

  # Supplement the table with the GO term name and ontology group
  # but only if the enrichment categories are actually GO terms
  if(any(grep("^GO:",pvals$category))){
    GOnames=select(GO.db,keys=pvals$category,columns=c("TERM","ONTOLOGY"))[,2:3]
    colnames(GOnames)<-tolower(colnames(GOnames))
    pvals=cbind(pvals,GOnames)
  }

  # And return
  return(pvals)
}

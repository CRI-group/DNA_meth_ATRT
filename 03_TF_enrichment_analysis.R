
#helper function to retrieve given data from GRanges based on a variable								
subsetter = function(gr, cname) {
  return(mcols(gr)[[cname]])
}


#keep only somatic chromosomes and sort the Granges
remove_unwanted_chromosomes_inGranges = function(granges){
  
  cat("\n")
  cat("Preprocessing GRanges. Original levels: \n")
  print(seqnames(granges))
  print(seqlevels(granges))
  
  filtered.granges = keepSeqlevels(granges, paste0("chr", c(1:22)), pruning.mode="coarse")
  
  cat("Preprocessing GRanges. Remaining levels: \n")
  print(seqnames(filtered.granges))
  print(seqlevels(filtered.granges))
  
  cat("\n")
  
  #sorting the ranges
  filtered.granges <- sortSeqlevels(filtered.granges)
  filtered.granges <- sort(filtered.granges)
  
  return(filtered.granges)
  
}

# For extending start and end coordinates of Granges object:
#https://support.bioconductor.org/p/78652/
extend <- function(x, upstream=0, downstream=0)     
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}




##Input regions expected in GRanges-format
fisher_test_of_TF_binding = function(tf.of.interest, tf.regions, background.regions, 
                                     regions.of.interest, tf.col.in.tf.regions="tfTitle"){
  
  
  tf.specific.regions = tf.regions[subsetter(tf.regions,tf.col.in.tf.regions) == tf.of.interest]
  
  dm.tf.bound = length(background.regions[background.regions %over% tf.specific.regions
                                          & background.regions %over% regions.of.interest])
  dm.tf.unbound = length(background.regions[!background.regions %over% tf.specific.regions
                                            & background.regions %over% regions.of.interest])
  normal.tf.bound = length(background.regions[background.regions %over% tf.specific.regions
                                              & ! background.regions %over% regions.of.interest])
  normal.tf.unbound = length(background.regions[!background.regions %over% tf.specific.regions 
                                                & ! background.regions %over% regions.of.interest])
  
  test.matrix <- matrix(c(dm.tf.bound, dm.tf.unbound, normal.tf.bound, normal.tf.unbound), 2, 2, byrow=T)
  p = fisher.test(test.matrix, alternative = 'g')$p.value
  
  result <- list(fisher.p = p, fisher.matrix = test.matrix)
  return(result)
  
}


##Run test for each TF using regions of interest. Extracts hypermethylated and hypomethylated regions. Runs the tests for
# 1) all regions, 2) hypermethylated regions, 3) hypomethylated regions. 
fisher_tests_for_regions_of_interest = function(regions.of.interest, binding.site.db.granges, tf.column, background.regions) {
  
  
  #methylkit data
  if( is.element("mcols.meth.diff", names(mcols(regions.of.interest)))){
    #Extracting also hyper and hypomethylated regions:
    hypo.regions.of.interest = regions.of.interest[regions.of.interest$mcols.meth.diff < 0] #only hypos!
    hyper.regions.of.interest = regions.of.interest[regions.of.interest$mcols.meth.diff > 0] #only hypers!
  }
  
  #DMRcate data
  if( is.element("meanbetafc", names(mcols(regions.of.interest)))){
    #Extracting also hyper and hypomethylated regions:
    hypo.regions.of.interest = regions.of.interest[regions.of.interest$meanbetafc < 0] #only hypos!
    hyper.regions.of.interest = regions.of.interest[regions.of.interest$meanbetafc > 0] #only hypers!
  }
  
  # Extracting the list of TFs to be tested 
  tf.list = unique(subsetter(binding.site.db.granges, tf.column))
  
  # Performing the test for each TF separately and extracting p-value vector and matrix lists:
  TF.pvals= list()
  TF.matrices = list()
  
  
  ### both hypers and hypos
  result.all = sapply(tf.list, fisher_test_of_TF_binding,
                      tf.regions=binding.site.db.granges, background.regions=background.regions,
                      regions.of.interest=regions.of.interest)
  
  TF.pvals$pvals.both = unlist(result.all["fisher.p", ]) 
  TF.matrices$both = result.all["fisher.matrix", ]
  
  
  
  if( is.element("mcols.meth.diff", names(mcols(regions.of.interest))) | is.element("meanbetafc", names(mcols(regions.of.interest)))){
    ### Only hypos
    result.hypo = sapply(tf.list, fisher_test_of_TF_binding,
                         tf.regions=binding.site.db.granges, background.regions=background.regions,
                         regions.of.interest=hypo.regions.of.interest)
    
    TF.pvals$pvals.hypo = unlist(result.hypo["fisher.p", ])
    TF.matrices$hypo = result.hypo["fisher.matrix", ]
    
    
    ### Only hypers
    result.hyper = sapply(tf.list, fisher_test_of_TF_binding,
                          tf.regions=binding.site.db.granges, background.regions=background.regions,
                          regions.of.interest=hyper.regions.of.interest)
    
    TF.pvals$pvals.hyper = unlist(result.hyper["fisher.p", ])
    TF.matrices$hyper = result.hyper["fisher.matrix", ]
    
  }
  
  return(list(pvals = TF.pvals, matrices = TF.matrices))
  
}



# Function to run The tests for each comparison (members in DMRs.for.groups list). 
# Not parallelized, thus requires lot of time and is also memory-greed. Can also test annotatr annotation subsets, 
# but this will take even longer, since then the amount of tests to calculate is multiplied.
# All region inputs (DMRs.for.groups, background.DMRs, binding.site.db)	as GRanges!
# Output is a nested list.
fisher_tests_for_groups = function(DMRs.for.groups, background.DMRs, binding.site.db, tf.column, analyze.annotatr.subsets=FALSE){
  
  
  if(analyze.annotatr.subsets){
    # #Building annotations 
     assembly.name = "hg38"
     annotatr.annotations = c(paste0(assembly.name,"_cpgs"), paste0(assembly.name, "_basicgenes"),
                              paste0(assembly.name, "_genes_intergenic"), paste0(assembly.name,"_enhancers_fantom"))
     annos = build_annotations(genome=assembly.name, annotations = annotatr.annotations)
     annos = keepStandardChromosomes(annos, pruning.mode = "coarse")
     
    
  }
  
  TF.pvals = list()
  TF.matrices = list()
  
  for(i in 1:length(DMRs.for.groups)){
    
    TF.pvals[[i]] = list()
    TF.matrices[[i]] = list()
    
    
    
    print("Analyzing all DMRs")
    
    pvals.and.matrices = fisher_tests_for_regions_of_interest(DMRs.for.groups[[i]], binding.site.db,
                                                              tf.column, background.DMRs)
    
    
    TF.pvals[[i]]$all.DMRs = pvals.and.matrices$pvals
    TF.matrices[[i]]$all.DMRs = pvals.and.matrices$matrices
    
    
    #    ************************************************************************************************
    if(analyze.annotatr.subsets){  
      
      ## NOTE that for each subsequent tests, the background is recalculated (based on annotations!!)
      
      print("Analyzing CpG islands in DMRs")#CpG-islands
      CpG.islands = annos[subsetter(annos, "type") == "hg38_cpg_islands"]
      CpG.islands = keepStandardChromosomes(CpG.islands, pruning.mode = "coarse")
      
      background = annotate_regions(background.DMRs, annotations = CpG.islands)
      
      CpG.islands = annotate_regions(DMRs.for.groups[[i]], annotations = CpG.islands) #, minoverlap = anno.min.overlap)
      
      pvals.and.matrices = fisher_tests_for_regions_of_interest(CpG.islands, binding.site.db, tf.column, background)
      TF.pvals[[i]]$CpG.islands = pvals.and.matrices$pvals
      TF.matrices[[i]]$CpG.islands = pvals.and.matrices$matrices
      # #*************************************************************************************************
      print("Analyzing enhancers in DMRs") #Enhancers
      enhancers = annos[subsetter(annos, "type") == "hg38_enhancers_fantom"]
      enhancers = keepStandardChromosomes(enhancers, pruning.mode = "coarse")
      
      background = annotate_regions(background.DMRs, annotations = enhancers)
      
      enhancers = annotate_regions(DMRs.for.groups[[i]], annotations = enhancers)
      
      pvals.and.matrices = fisher_tests_for_regions_of_interest(enhancers, binding.site.db, tf.column, background)
      TF.pvals[[i]]$enhancers = pvals.and.matrices$pvals
      TF.matrices[[i]]$enhancers = pvals.and.matrices$matrices
      
      #***************************************************************************************************
      ## For promoters and genomic neighbourhoods, the default promoters of annotatr (1kb from TSS)
      ## needs to be extended:
      promoters1kbTSS = annos[subsetter(annos, "type") == "hg38_genes_promoters"]
      promoters1kbTSS = keepStandardChromosomes(promoters1kbTSS, pruning.mode = "coarse")
      promoters1kbTSS = subset(promoters1kbTSS, !is.na(promoters1kbTSS$symbol))
      
      #***************************************************************************************************
      print("Analyzing promoters (2kb upstream and 500b downstream from TSS) in DMRs")
      promoters2kbTSS = extend(promoters1kbTSS, upstream=1000, downstream=500)
      elementMetadata(promoters2kbTSS)[[ "dist.to.TSS" ]] = "2kb_upstr_500b_downstr"
      names(elementMetadata(promoters2kbTSS)) = c("id", "tx_id", "gene_id", "symbol", "annot.type","dist.to.TSS")
      
      background = annotate_regions(background.DMRs, annotations = promoters2kbTSS)
      
      promoters2kbTSS = annotate_regions(DMRs.for.groups[[i]], annotations = promoters2kbTSS)
      
      pvals.and.matrices = fisher_tests_for_regions_of_interest(promoters2kbTSS, binding.site.db, tf.column, background)
      TF.pvals[[i]]$promoters = pvals.and.matrices$pvals
      TF.matrices[[i]]$promoters = pvals.and.matrices$matrices
      #*****************************************************************************************************
      print("Analyzing genomic neighbourhoods (200kb upstream and downstream from TSS) in DMRs")#Genomic neighbourhoods 200 kb
      neighbourhoods200kb = extend(promoters1kbTSS, upstream=199000, downstream=200000)
      neighbourhoods200kb = keepStandardChromosomes(neighbourhoods200kb, pruning.mode = "coarse") #trimming extras
      background = annotate_regions(background.DMRs, annotations = neighbourhoods200kb)
      neighbourhoods200kb = annotate_regions(DMRs.for.groups[[i]], annotations = neighbourhoods200kb)
      
      pvals.and.matrices = fisher_tests_for_regions_of_interest(neighbourhoods200kb,
                                                                binding.site.db,
                                                                tf.column,background)
      TF.pvals[[i]]$genomic.neighbourhoods200kb = pvals.and.matrices$pvals
      TF.matrices[[i]]$genomic.neighbourhoods200kb = pvals.and.matrices$matrices
      
    }
  }
  
  
  names(TF.pvals) = names(DMRs.for.groups)
  names(TF.matrices) = names(DMRs.for.groups)
  
  result <- list(TFBS.pvals = TF.pvals,
                 TFBS.fisher.matrices = TF.matrices)
  
  return(result)
  
}




mark_cancer_specific_concordant_areas_for_DMRcate_probes = function(comparison1,comparison2){
  olaps <- findOverlaps(comparison1,comparison2, minoverlap=600) 
  
  res <- comparison2[subjectHits(olaps)] # using tiles as coordinates
  mcols(res) <- mcols(comparison1[queryHits(olaps)])
  
  res2 <- comparison1[queryHits(olaps)] # using probes as coordinates
  #mcols(res2) <- cbind(mcols(res2),data.frame(comparison2[subjectHits(olaps)]))
  mcols(res2) <- data.frame(comparison2[subjectHits(olaps)])
  #concordant=res[sign(res$meanbetafc) == sign(res2$meanbetafc)]
  
  elementMetadata(res)[[ "DM.status" ]] = 
    ifelse( res$meanbetafc > 0, "hyper", "hypo")
  elementMetadata(res)[[ "Behavior" ]] = 
    ifelse( sign(res$meanbetafc) == sign(res2$meanbetafc), "concordant", "non-corcondant")
  
  res = unique(res) 
  return(res)
}


###
output.dir = "PATH/FOR/ANANLYSIS/RESULTS"

all.gtrd.tfs = read.table("/YOUR/PATH/Homo_sapiens_meta_clusters.intervals.GTRD19_TFs_w_names.bed", sep='\t', header=F, stringsAsFactors = F)
#The input "Homo_sapiens_meta_clusters.intervals.GTRD19_TFs_w_names.bed" generated from the full GTRD 19.04 metacluster file (Homo_sapiens_meta_clusters.interval.gz) 
#by taking the genomic coordinate columns and tfTitle
colnames(all.gtrd.tfs) = c("chrom", "start", "end", "tfTitle")

all.gtrd.tfs = makeGRangesFromDataFrame(all.gtrd.tfs, keep.extra.columns=TRUE,
                                        ignore.strand=TRUE, seqinfo=NULL,
                                        seqnames.field="chrom",
                                        start.field="start",
                                        end.field="end")

all.gtrd.tfs = remove_unwanted_chromosomes_inGranges(all.gtrd.tfs)



#Producing background tiles
## Data for 450k
heatmap.hg38.probes.Granges = readRDS("YOUR/PATH/hg38.probes_for_heatmaps.rds")
probe.regions = remove_unwanted_chromosomes_inGranges(heatmap.hg38.probes.Granges) #This also sorts the ranges
probe.regions = extend(heatmap.hg38.probes.Granges, 500, 500)


#Reading DMRs
DMR.dir = "PATH/TO/DMRs"
DMR.names = list.files(DMR.dir)[grep("_located_DMRcate_GRanges_filtered_by_controls.rds", list.files(DMR.dir))]
dmrcate.DMRs = lapply(file.path(DMR.dir, DMR.names), readRDS)
names(dmrcate.DMRs) = gsub("_located_DMRcate_GRanges_filtered_by_controls.rds", "", DMR.names)

dmrcate.DMRs = lapply(dmrcate.DMRs, remove_unwanted_chromosomes_inGranges)
dmrcate.DMRs = lapply(dmrcate.DMRs, make_expanded_ranges_for_450k_DMR_probes,  probe.Granges=heatmap.hg38.probes.Granges)




###Preparing Cancer-specific regions for tests

ATRT = mark_cancer_specific_concordant_areas_for_DMRcate_probes(dmrcate.DMRs$ATRT.vs.MB, dmrcate.DMRs$ATRT.vs.PLEX)


MB.vs.ATRT.tmp = dmrcate.DMRs$ATRT.vs.MB
MB.vs.ATRT.tmp$meanbetafc = MB.vs.ATRT.tmp$meanbetafc*(-1)

MB.vs.PLEX.tmp = dmrcate.DMRs$PLEX.vs.MB
MB.vs.PLEX.tmp$meanbetafc = MB.vs.PLEX.tmp$meanbetafc*(-1)

MB = mark_cancer_specific_concordant_areas_for_DMRcate_probes(MB.vs.ATRT.tmp, MB.vs.PLEX.tmp)

PLEX.vs.ATRT.tmp = dmrcate.DMRs$ATRT.vs.PLEX
PLEX.vs.ATRT.tmp$meanbetafc = PLEX.vs.ATRT.tmp$meanbetafc*(-1)

PLEX = mark_cancer_specific_concordant_areas_for_DMRcate_probes(dmrcate.DMRs$PLEX.vs.MB, PLEX.vs.ATRT.tmp)

DMRs.for.groups.granges = list(ATRT = subset(ATRT, ATRT$Behavior == "concordant"),
                               MEDULLO = subset(MB, MB$Behavior == "concordant"),
                               PLEXUS = subset(PLEX, PLEX$Behavior == "concordant"))


saveRDS(DMRs.for.groups.granges, file.path(DMR.dir, "Cancer_DMRs.RDS"))


#### running tests and storing the results:
group.analysis.with.GTRD.full = fisher_tests_for_groups(DMRs.for.groups=DMRs.for.groups.granges,
                                                        binding.site.db=all.gtrd.tfs,
                                                        tf.column="tfTitle", background.DMRs = probe.regions)
saveRDS(group.analysis.with.GTRD.full, file.path(output.dir, "Fisher_results_with_GTRD.full_for.cancers_450k.RDS"))


###############################################

#Examples for further usage:

correct_nested_list_of_pvals = function(TF.pvals, pval.thr, method = "BH"){
  adjusted.pvals = TF.pvals
  for(roi in 1:length(names(TF.pvals[[1]]))){
    
    for (binding in 1:length(names(TF.pvals[[1]][[1]]))){
      
      
      for(comparison in 1:length(TF.pvals)){
        #names = names(subset.TF.pvals)[comparison]
        pvalues = TF.pvals[[comparison]][[roi]][[binding]]
        print(paste0(names(TF.pvals)[comparison], ":", names(TF.pvals[[comparison]])[roi], ":", names(TF.pvals[[comparison]][[roi]])[binding]))
        #cat(paste0(sum(pvalues < pval.thr), "significant in original \n"))
        pvals.corrected = p.adjust(pvalues, method = method)
        #cat(paste0(sum(pvals.corrected < pval.thr), "significant after ", method, " correction \n"))
        adjusted.pvals[[comparison]][[roi]][[binding]] =  pvals.corrected
        
      }
      
    }
  }	
  return(adjusted.pvals)
  
}


#read and adjust P-vals
group.TF.pvals = group.analysis.with.GTRD.full$TFBS.pvals 
group.TF.pvals = correct_nested_list_of_pvals(TF.pvals=group.TF.pvals, pval.thr=0.001, method = "BH") #NOTE! Using adjusted pvals

#Read Test matrices
group.TF.test.matrices =group.analysis.with.GTRD.full$TFBS.fisher.matrices













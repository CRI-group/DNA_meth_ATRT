library(GenomicRanges)
library(annotatr)


#Helper functions:

## mark DMRs as "hyper" or "hypo" methylated based on meth.diff column sign
add_DM.status_column = function(dmr.as.GRanges, meth.diff.column){ 
  elementMetadata(dmr.as.GRanges)[[ "DM.status" ]] = 
    ifelse( subsetter(dmr.as.GRanges, meth.diff.column) > 0, "hyper", "hypo")
  return(dmr.as.GRanges)
}


#Retrieve given data from GRanges based on a variable								
subsetter = function(gr, cname) {
  return(mcols(gr)[[cname]])
}


#To get cleaner output
prepare_DMRcate_Granges_for_DE_intergration = function(dmr.as.GRanges  = NULL){
  
  elementMetadata(dmr.as.GRanges)[["overlapping.promoters"]] = NULL 
  dmr.as.GRanges = add_DM.status_column(dmr.as.GRanges, "meanbetafc")
  
  return(dmr.as.GRanges) 
  
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



#Get annotatr annotations as data frame format
produce_full_annotatr_annotations = function(dmr.as.GRanges = NULL, annotations = NULL, 
                                             anno.min.overlap = 1L){
  dm.annotated = annotate_regions(dmr.as.GRanges, annotations = annotations,
                                  minoverlap = anno.min.overlap, quiet = TRUE)
  df.annotated = data.frame(dm.annotated)
  
  return(df.annotated)
  
}


# Integrate annotated DMRs (given as data frame format from produce_full_annotatr_annotations-function ).
# If DE-comparison was opposite to DM-comparison (numerator and denominator switched place), DE-fold change can be inverted
integrate_annotated_DMRs_with_DE = function(annotatr.df, de.data, invert.foldchange.sign = TRUE, 
                                            logFC.col = "log2FoldChange", symbol.col = "Gene.symbol",
                                            DM.status.col = "DM.status", region.descr="promoter"){
  
  if(invert.foldchange.sign){
    de.data[[logFC.col]] = de.data[[logFC.col]] * (-1)
  }
  
  DE.DMR = annotatr.df[which(annotatr.df$annot.symbol %in% de.data[[symbol.col]]),]
  colnames(DE.DMR)[grep("annot.symbol", colnames(DE.DMR))] = symbol.col
  merged.DE.DMR = merge(x = DE.DMR, y = de.data, by = symbol.col, all.x = TRUE)
  
  hypers.down = merged.DE.DMR[(merged.DE.DMR[[DM.status.col]] == "hyper" & 
                                 merged.DE.DMR[[logFC.col]] < 0), ]
  
  hypos.up = merged.DE.DMR[(merged.DE.DMR[[DM.status.col]] == "hypo" & 
                              merged.DE.DMR[[logFC.col]] > 0), ]
  
  hypers.up = merged.DE.DMR[(merged.DE.DMR[[DM.status.col]] == "hyper" & 
                               merged.DE.DMR[[logFC.col]] > 0), ]
  
  hypos.down = merged.DE.DMR[(merged.DE.DMR[[DM.status.col]] == "hypo" & 
                                merged.DE.DMR[[logFC.col]] < 0), ]
  
  cat("Genes that seem to be regulated by methylation:\n")
  
  cat(paste0(length(unique(hypers.down[[symbol.col]])), " genes with hypermethylated ", region.descr, " and expression down \n"))
  
  cat(paste0(length(unique(hypos.up[[symbol.col]])), " genes with hypomethylated ", region.descr, " and expression up \n"))
  
  cat(paste0("***", "\n"))
  
  cat("Conflicts (currently not filtered out) :\n")
  
  cat(paste0(length(unique(hypers.up[[symbol.col]])), " genes with hypermethylated ", region.descr, " and expression up \n"))
  
  cat(paste0(length(unique(hypos.down[[symbol.col]])), " genes with hypomethylated ", region.descr, " and expression down \n"))
  
  return(merged.DE.DMR)
  
}


#NOTE!!! expects to have annotatr default promoters as an input!
broaden_annoations = function(input.annotatr.annotation.set =promoters1kbTSS, annotation.type = "promoter", 
                              bases.upstream = 1000, bases.downstream = 0){
  
  kb.upstream = bases.upstream/1000 +1
  kb.downstream = bases.downstream/1000
  dist.to.TSS.value = paste0(kb.upstream, "kb_upstr_",  kb.downstream, "kb_downstr" )
  
  new.annotation.set = extend(input.annotatr.annotation.set, upstream=bases.upstream, downstream=bases.downstream)
  elementMetadata(new.annotation.set)[[ "dist.to.TSS" ]] = dist.to.TSS.value
  names(elementMetadata(new.annotation.set)) = c("id", "tx_id", "gene_id", "symbol", "annot.type","dist.to.TSS")
  
  elementMetadata(new.annotation.set)[[ "annot.type" ]] = annotation.type
  
  return(new.annotation.set)
  
}



### Function for annotating several comparisons (members of DMRs.granges.list) at once
annotate_and_integrate_DMR_objects_with_DE_tables = function(DMRs.granges.list, de.tables, annotatr.annotation.set.granges, region.description, 
                                                             invert.foldchange.sign = TRUE, symbol.col = "gene_symbol", logFC.col="log2FoldChange", qvalue.cutoff, meth.change.cutoff, output.dir=NA, used.DMR.method="methylKit" ){
  
  merged.DE.DMRs = list()
  for(i in 1:length(DMRs.granges.list)){
    
    print(paste0("Integrating ", names(de.tables)[i], " DEGs with ", names(DMRs.granges.list)[i], " DMRs"))
    
    
    annotated.df = produce_full_annotatr_annotations(dmr.as.GRanges = DMRs.granges.list[[i]], 
                                                     annotations = annotatr.annotation.set.granges, 
                                                     anno.min.overlap = 1L)
    
    #For some reason, the word "annot" is duplicated when using customized. annotations. Let's ensure the data is always similarly named:
    colnames(annotated.df) = gsub("annot.annot.", "annot.", colnames(annotated.df))
    
    merged.DE.DMR = integrate_annotated_DMRs_with_DE(annotatr.df = annotated.df, 
                                                     de.data = de.tables[[i]], invert.foldchange.sign = invert.foldchange.sign,
                                                     symbol.col = symbol.col, logFC.col = logFC.col, region.descr = region.description)
    
    if(!is.na(output.dir)){
      
      fname=paste0(names(DMRs.granges.list[i]), "_", region.description,"_DMR_DE_table_", used.DMR.method, "_q", qvalue.cutoff, "_methDiff",  
                   meth.change.cutoff,".tsv")
      
      fname = file.path(output.dir, fname)
      write.table(merged.DE.DMR, fname, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    merged.DE.DMRs[[i]] =  merged.DE.DMR
  }
  
  if(is.na(output.dir)){
    names(merged.DE.DMRs) = names(DMRs.granges.list)
    return(merged.DE.DMRs)
  }
}

#############################
#Producing annotation inputs:

library(annotatr)
 assembly.name = "hg38"
 annotatr.annotations = c(paste0(assembly.name,"_cpgs"), paste0(assembly.name, "_basicgenes"), paste0(assembly.name, "_genes_intergenic"), paste0(assembly.name,"_enhancers_fantom"))
 annos = build_annotations(genome=assembly.name, annotations = annotatr.annotations)
 annos = keepStandardChromosomes(annos, pruning.mode = "coarse")

enhancers = annos[subsetter(annos, "type") == "hg38_enhancers_fantom"]
CpGislands= annos[subsetter(annos, "type") == "hg38_cpg_islands"]

#annos[subsetter(annos, "type") == "hg38_genes_promoters"]
genomic.annos = annos[grep("genes", subsetter(annos, "type"))]
genomic.annos = subset(genomic.annos, !is.na(genomic.annos$symbol))

promoters1kbTSS = genomic.annos[subsetter(genomic.annos, "type") == "hg38_genes_promoters"]
promoters2kbTSS = broaden_annoations(input.annotatr.annotation.set = promoters1kbTSS,annotation.type = "promoter", bases.upstream = 1000, bases.downstream = 500)

TSSes = broaden_annoations(input.annotatr.annotation.set = promoters1kbTSS,annotation.type = "TSS", bases.upstream = -999, bases.downstream = 0)

neighbourhoods200kb = broaden_annoations(input.annotatr.annotation.set = promoters1kbTSS,annotation.type = "genomic_neighbourhood", bases.upstream = 199000, bases.downstream = 200000)
neighbourhoods200kb = keepStandardChromosomes(neighbourhoods200kb, pruning.mode = "coarse")
#Adding TSS coordinates:
elementMetadata(neighbourhoods200kb)[[ "TSS" ]] = start(TSSes)



####FANTOM5 annotations:

## Downloading and processing fantom5 data for promoter-enhancer analysis
#following protocol presented in 
#https://github.com/enricoferrero/bioconductor-regulatory-genomics-workflow/blob/master/biocondutor-regulatory-genomics-workflow.md


#download.file("http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed", destfile = "enhancer_tss_associations.bed")
fantom <- read.delim("UPDATE/PATH/enhancer_tss_associations.bed", skip = 1, stringsAsFactors = FALSE)
head(fantom, 3)


#Parsing relevant fields:
library(splitstackshape)
fantom <- as.data.frame(cSplit(fantom, splitCols = "name", sep = ";", direction = "wide"))

locs <- strsplit(as.character(fantom$name_1), "[:-]")
fantom$chr <- sapply(locs, "[", 1)
fantom$start <- as.numeric(sapply(locs, "[", 2))
fantom$end <- as.numeric(sapply(locs, "[", 3))
fantom$gene_symbol <- fantom$name_3
fantom$corr <- sub("R:", "", fantom$name_4)
fantom$fdr <- sub("FDR:", "", fantom$name_5)

#Filtering and tidying the data:
fantom <- unique(subset(fantom, corr >= 0.25 & fdr < 0.01, select = c("chr", "start", "end", "gene_symbol")))


#Liftover
library(liftOver)
#Chain file for liftover
hg19ToHg38 ="UPDATE/PATH/hg19ToHg38.over.chain"
ch = import.chain(hg19ToHg38, exclude = "_")


fantom <- makeGRangesFromDataFrame(fantom, keep.extra.columns = TRUE)
fantom <- unlist(liftOver(fantom, ch))


########
#Annotation examples, using 450K data:
genehancer.confidence.score.min = 5 ##As in https://dash.harvard.edu/bitstream/handle/1/34491893/5606860.pdf?sequence=1
qvalue.cutoff = 0.05    #To filter DMRs
meth.change.cutoff = 5  #To filter DMRs


#Reading DMRs
DMR.dir = "PATH/TO/DMR_files"
DMR.names = list.files(DMR.dir)[grep("_located_DMRcate_GRanges_filtered_by_controls.rds", list.files(DMR.dir))]
dmrcate.DMRs = lapply(file.path(DMR.dir, DMR.names), readRDS)
names(dmrcate.DMRs) = gsub("_located_DMRcate_GRanges_filtered_by_controls.rds", "", DMR.names)

dmrcate.DMRs = lapply(dmrcate.DMRs, prepare_DMRcate_Granges_for_DE_intergration)
output.dir.450k ="/YOUR/PATH"

###reading DE-data

#Illumina expression array results
public_DEGs = list()
public_DEGs$ATRT.vs.MB = read.table(file.path(DMR.dir, "GSE42658_ATRT.vs.MB_filtered_DEs_fdr0.05_logFC1.tsv"), sep = "\t")
public_DEGs$ATRT.vs.PLEX = read.table(file.path(DMR.dir, "GSE42658_ATRT.vs.PLEX_filtered_DEs_fdr0.05_logFC1.tsv"), sep = "\t")
public_DEGs$PLEX.vs.MB = read.table(file.path(DMR.dir, "GSE42658_PLEX.vs.MB_filtered_DEs_fdr0.05_logFC1.tsv"), sep = "\t")

#Loading affy data and preparing the union with illumina data:
affy.DE.input.paths = list()
affy.DE.input.paths[[1]] = "YOUR/PATH/GSE35493_ATRT.vs.MB_filtered_DEs_fdr0.05_logFC1.tsv"
affy.DE.input.paths[[2]] = "YOUR/PATH/GSE73038_ATRT.vs.MB_filtered_DEs_fdr0.05_logFC1.tsv"

#Collecting genes in both tables
affy.ATRT.MB.genes = list()
for( i in 1:length(affy.DE.input.paths)){
  affy.DE.input.path = affy.DE.input.paths[[i]]
  filtered.affy.ATRT.vs.MB = read.table(affy.DE.input.path, sep = "\t")
  genes = rownames(filtered.affy.ATRT.vs.MB)
  affy.ATRT.MB.genes[[i]] = genes
}

affy1 = read.table(affy.DE.input.paths[[1]], sep = "\t")
affy2 = read.table(affy.DE.input.paths[[2]], sep = "\t")

#affy.union.table = rbind(affy2, affy1[!(rownames(affy1) %in% rownames(affy2)),])

affy.intersect = intersect(affy.ATRT.MB.genes[[1]], affy.ATRT.MB.genes[[2]])
affy.intersect.table =  affy2[affy.intersect,]

#
new.ATRT.vs.MB.DEs = rbind(affy.intersect.table, public_DEGs$ATRT.vs.MB[!(rownames(public_DEGs$ATRT.vs.MB) %in% rownames(affy.intersect.table)),])

public_DEGs$ATRT.vs.MB = new.ATRT.vs.MB.DEs

#Adding gene symbol column
for(i in 1:length(public_DEGs)){
  public_DEGs[[i]]$gene_symbol = rownames(public_DEGs[[i]])
}


### A) Promoters & Genomic neighborhoods
annotate_and_integrate_DMR_objects_with_DE_tables(DMRs.granges.list=dmrcate.DMRs, 
                                                  de.tables=public_DEGs, 
                                                  annotatr.annotation.set.granges=promoters2kbTSS, 
                                                  region.description="promoters2kbTSS", 
                                                  invert.foldchange.sign = FALSE, 
                                                  symbol.col = "gene_symbol", 
                                                  logFC.col="logFC",
                                                  qvalue.cutoff=qvalue.cutoff, 
                                                  meth.change.cutoff=meth.change.cutoff, 
                                                  output.dir = output.dir.450k, 
                                                  used.DMR.method = "DMRcate")


annotate_and_integrate_DMR_objects_with_DE_tables(DMRs.granges.list=dmrcate.DMRs, 
                                                  de.tables=public_DEGs, 
                                                  annotatr.annotation.set.granges=neighbourhoods200kb, 
                                                  region.description="neighbourhoods200kb", 
                                                  invert.foldchange.sign = FALSE, 
                                                  symbol.col = "gene_symbol", 
                                                  logFC.col="logFC",
                                                  qvalue.cutoff=qvalue.cutoff, 
                                                  meth.change.cutoff=meth.change.cutoff, 
                                                  output.dir = output.dir.450k, 
                                                  used.DMR.method = "DMRcate")





# B) Fantom5 enhancers,
# see also https://github.com/enricoferrero/bioconductor-regulatory-genomics-workflow/blob/master/biocondutor-regulatory-genomics-workflow.md
for(i in 1:length(dmrcate.DMRs)){
  print(paste0("Comparing ", names(public_DEGs)[i], " with ", names(dmrcate.DMRs)[i], " DMRs"))
  region.description = "fantom5enhancer_mapped_to_TSS"
  
  hits <- findOverlaps(dmrcate.DMRs[[i]], fantom)
  DMRs.in.fantom = dmrcate.DMRs[[i]][queryHits(hits)]
  fantom_with_DMRs = fantom[subjectHits(hits)]
  mcols(DMRs.in.fantom) <- cbind(mcols(DMRs.in.fantom), mcols(fantom_with_DMRs))
  names(DMRs.in.fantom) = NULL
  DMRs.in.fantom <- as.data.frame(DMRs.in.fantom)
  
  fantom.DMRs.and.DE = merge(DMRs.in.fantom, public_DEGs[[i]], by = "gene_symbol", all = FALSE)
  fantom.DMRs.and.DE$annot.type = region.description #Adding annot.type
  
  fname=paste0(names(dmrcate.DMRs[i]), "_", region.description,"_DMR_DE_table_DMRcate_q", qvalue.cutoff, "_methDiff",
               meth.change.cutoff,".tsv")
  fname = file.path(output.dir.450k, fname)
  
  write.table(fantom.DMRs.and.DE, fname, sep = "\t", quote = FALSE, row.names = FALSE)
  
}

# C) Genehancer 
genehancer.annos=read.xlsx("/UPDATE/PATH/genehancer.xlsx")  #From  [https://genecards.weizmann.ac.il/geneloc_prev/genehancer.xlsx],
genehancer.annos=as(genehancer.annos, "GRanges")

if(filter.genhancer){
  genehancer.annos=genehancer.annos[grep("Promoter/Enhancer", subsetter(genehancer.annos, "feature.name"))]
}

for(i in 1:length(dmrcate.DMRs)){
  print(paste0("Comparing ", names(public_DEGs)[i], " with ", names(dmrcate.DMRs)[i], " DMRs"))
  region.description = "genehancer_mapped_to_TSS"
  
  hits <- findOverlaps(dmrcate.DMRs[[i]], genehancer.annos)
  DMRs.in.genehancer = dmrcate.DMRs[[i]][queryHits(hits)]
  genehancer_with_DMRs = genehancer.annos[subjectHits(hits)]
  mcols(DMRs.in.genehancer) <- cbind(mcols(DMRs.in.genehancer), mcols(genehancer_with_DMRs))
  names(DMRs.in.genehancer) = NULL
  DMRs.in.genehancer <- as.data.frame(DMRs.in.genehancer)
  
  genes.of.interest = public_DEGs[[i]]$gene_symbol
  
  relevant.row=c()
  gene.symbols=c()
  for(j in 1:nrow(DMRs.in.genehancer)){
    records.of.genehancer = DMRs.in.genehancer$attributes[j]
    records.of.genehancer = unlist(strsplit(records.of.genehancer, ";"))
    genes.of.genehancer = records.of.genehancer[grep("connected_gene", records.of.genehancer)]
    scores.for.genehancer=records.of.genehancer[grep("score", records.of.genehancer)]
    
    
    genes.of.genehancer = gsub("connected_gene=", "", genes.of.genehancer)
    scores.for.genehancer=as.numeric(gsub("score=", "", scores.for.genehancer))
    
    #Filtering by score:
    genes.of.genehancer= genes.of.genehancer[scores.for.genehancer > genehancer.confidence.score.min]
    
    no.match = sum(genes.of.genehancer %in% genes.of.interest) == 0
    if(no.match){
      relevant.row =c(relevant.row,FALSE)
    } else {
      relevant.row = c(relevant.row, TRUE)
      gene.symbols = c(gene.symbols, paste(genes.of.genehancer[genes.of.genehancer %in% genes.of.interest], collapse = ","))
      
      #print(length(genes.of.genehancer[genes.of.genehancer %in% genes.of.interest]))
    }
    #genes.of.genehancer[genes.of.genehancer %in% genes.of.interest]
  }
  
  #filtering
  DMRs.in.genehancer= DMRs.in.genehancer[relevant.row,]
  #adding column 
  DMRs.in.genehancer$gene_symbol=gene.symbols
  #DMRs.in.genehancer$attributes = NULL
  DMRs.in.genehancer = DMRs.in.genehancer %>% tidyr::separate_rows(1:length(colnames(DMRs.in.genehancer)), sep=",")
  
  
  genehancer.DMRs.and.DE = merge(DMRs.in.genehancer, public_DEGs[[i]], by = "gene_symbol", all = FALSE)
  genehancer.DMRs.and.DE$annot.type = region.description #Adding annot.type
  
  fname=paste0(names(dmrcate.DMRs[i]), "_", region.description,"_DMR_DE_table_DMRcate_q", qvalue.cutoff, "_methDiff",
               meth.change.cutoff,".tsv")
  fname = file.path(output.dir.450k, fname)
  
  if(filter.genhancer){
    fname=gsub(".tsv", "_most_confident.tsv", fname)
  }
  
  write.table(genehancer.DMRs.and.DE, fname, sep = "\t", quote = FALSE, row.names = FALSE)
  
}


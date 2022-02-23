library(GenomicRanges)
library(methylKit)
library(scales)
library(ggplot2)
library(reshape2)


add_DM.status_column = function(dmr.as.GRanges, meth.diff.column){ 
  elementMetadata(dmr.as.GRanges)[[ "DM.status" ]] = 
    ifelse( subsetter(dmr.as.GRanges, meth.diff.column) > 0, "hyper", "hypo")
  return(dmr.as.GRanges)
}


## prepare_methylkit_object()
#
#INPUTS: location {character}, path to methylKit-compatible txt.files
#        pattern.to.remove.from.names {character}, to extract sample IDs
#        assembly.name {character}, reference genome name in alignment 
#        sample.sheet {character/dataframe} has to contain col "Sample_id"
#                      matching with filenames after pattern is removed
#                      AND col "Group_number". Groups are good to be 
#                      numbered with the growing severity, but this only
#                      affects into comparison formation.
#        sample.sheet.separator, if sample sheet was given as a file, 
#                      give the mark how to split data into  columns  
#
#OUTPUT: methylKit::methRead()-return value, i.e. list of raw-GpC-data
#///////////////////////////////////////////////////////////////////  
prepare_methylkit_object <- function(location=getwd(), 
                                     pattern.to.remove.from.names="_CpG.txt$",
                                     assembly.name = "hg38", 
                                     sample.sheet=NULL, sample.sheet.separator=";"){
  
  if(typeof(sample.sheet) == "character") {
    sample.info = read.csv(sample.sheet, sep=sample.sheet.separator)
  }
  else if(typeof(sample.sheet) == "data.frame"){
    sample.info = sample.sheet
  }
  else {
    stop("Provide sample sheet as a dataframe or as path to csv-file")
    
  }
  
  id.column.exists = sum(grepl("Sample_id$", colnames(sample.info))) == 1
  if(!id.column.exists){
    stop("Provide sample sheet with the column 'Sample_id'")
  }
  
  group.column.exists = sum(grepl("Group_number$", colnames(sample.info))) == 1
  if(!group.column.exists){
    stop("Provide sample sheet with the column 'Group_number'")
  }
  
  
  pattern.to.search = pattern.to.remove.from.names
  
  
  files = list.files(path = location, 
                     pattern = pattern.to.search, 
                     full.names=TRUE)
  
  basenames = list.files(path = location, 
                         pattern = pattern.to.search,
                         full.names=FALSE)
  
  sample.ids =  sub(pattern.to.remove.from.names, '', basenames) 
  
  if( all(sample.ids %in% sample.info$Sample_id) ){
    
    #order the sample sheet based on filename to get groups in appropriate order
    rownames(sample.info) = sample.info$Sample_id
    sample.info = sample.info[sample.ids,]
    treatments = sample.info$Group_number
    
    
    methObj = methRead(location = as.list(files),
                       sample.id = as.list(sample.ids),
                       assembly = assembly.name, 
                       treatment = treatments)
    
    
    
    return(methObj)
    
    
  }
  else {
    stop("Tried to process samples not found in sample.info$Sample_id. 
                 Check the naming politics and pattern.to.remove.from.names")
  }
  
}





#preprocess_methylkit_object()
#
#INPUTS: methObj {methylKit::}
#
#OUTPUT: filtered and possibly normalised methylKit object
#//////////////////////////////////////////////////////////////
preprocess_methylkit_object = function(methObj=NULL, normalize=TRUE, filter=TRUE, min.cov.as.reads=10, max.cov.as.perc=99.9){
  
  if(do.filtering==TRUE){
    processed.methObj = filterByCoverage(methObj,lo.count=min.cov.as.reads,
                                         lo.perc=NULL,
                                         hi.count=NULL,
                                         hi.perc=max.cov.as.perc)
  } 
  else {
    processed.methObj = methObj 
  }
  
  if(do.normalisation){
    processed.methObj = normalizeCoverage(processed.methObj)
  }
  
  return( processed.methObj )
}



# prepare_pairwise_group_comparisons_from_united_methObj()
#
#///////////////////////////////////////////////////////////////
prepare_pairwise_group_comparisons_from_united_methObj = function(united.methObj, group.key){
  
  # 2-member combinations of the treatments,
  # order is not important (returns a matrix):
  comparisons = combn(unique(getTreatment(united.methObj)), 2)
  
  # Number of possible comparisons is ncol(comparisons)
  comparison.objects = list()
  comparison.names = c()
  for(i in 1:ncol(comparisons)){
    
    # Collecting information needed by methylKit::reorganize()
    
    #group.A = comparisons[1,i]
    group.A = max(comparisons[,i]) # Taking the severity into account...
    # ...(Setting the group with higher number as the numerator )
    group.A.samples = getSampleID(united.methObj)[getTreatment(united.methObj) == group.A]
    
    #...(Setting the group with smaller number as denominator )
    group.B = min(comparisons[,i])
    group.B.samples = getSampleID(united.methObj)[getTreatment(united.methObj) == group.B]
    
    sample.selector = c(group.A.samples, group.B.samples)
    treatment.selector = c(rep(group.A, length(group.A.samples)),
                           rep(group.B, length(group.B.samples)) )
    
    #Taking samles of each group
    #new.comp.object = reorganize(united.methObj, sample.ids = sample.selector,
    #                             treatment = treatment.selector)
    
    #print(typeof(new.comp.object))
    #comparison.objects[[i]] = new.comp.object
    comparison.objects[[i]] = reorganize(united.methObj, sample.ids = sample.selector,treatment = treatment.selector)
    print(typeof(comparison.objects[[i]]))
    
    
    group.A.name = as.character(
      group.key$Sample_group[group.key$Group_number == group.A])
    
    group.B.name = as.character(
      group.key$Sample_group[group.key$Group_number == group.B])
    
    comparison.names[i] = paste0(group.A.name, ".vs.", group.B.name)
    #print(comparison.names[i])
  }
  
  #print(length(comparison.objects))
  #print(length(comparison.names))
  names(comparison.objects) = comparison.names 
  #sapply(comparison.objects, typeof)
  
  return(comparison.objects)
  
}

methylDiff_object_to_GRanges = function(methylDiff_object, assembly.name = "hg38"){
  
  dmr.as.GRanges = as(methylDiff_object, "GRanges")
  
  #Removing pseudochromosomes and sex chromosomes
  dmr.as.GRanges = keepSeqlevels(dmr.as.GRanges, paste0("chr", c(1:22)), pruning.mode="coarse")
  
  dmr.as.GRanges = add_DM.status_column(dmr.as.GRanges, "meth.diff")
  
  genome(dmr.as.GRanges) = assembly.name
  
  return(dmr.as.GRanges)
}



#########################
data.path="/PATH/TO/methylKit_CpGs/"
sample.sheet="/PATH/TO/Sample_info_RRBS.csv" 
sample.info=as.data.frame(read.csv(sample.sheet, sep=";"))

output.dir="/YOUR/PATH/results"
#dir.create(output.dir)


#ANALYSIS VARIABLES:

available.cores = 40
assembly.name="hg38"

do.filtering=TRUE
do.normalisation=TRUE
save.DMRs = FALSE #Whether to save results into disk 

#Filtering:
min.coverage.as.reads=10 #DEFAULT: 10
max.coverage.as.perc=99.8 #DEFAULT: 99.9 (percentile)

#Comparative analysis:
remove.samples.from.tests = FALSE #Whether there is a known sample(s) that should be removed
samples.to.remove = "" #Value(s) from Sample_Id column!!



#Note that the comparisons are generated such that the group with 
#the higher number is set as A in comparison  A.vs.B (or A/B)
group.name.number.key = sample.info[, c("Sample_group", "Group_number")]
group.name.number.key = subset(group.name.number.key,
                               !duplicated(group.name.number.key$Sample_group))
rownames(group.name.number.key) = as.integer(group.name.number.key$Group_number)
group.name.number.key = group.name.number.key[order(group.name.number.key$Group_number),]


group.name.colour.key = data.frame("Sample_group" =
                                     group.name.number.key$Sample_group,
                                   "Color" = c("#E69F00", "#009E73", "#0072B2"))  #as.character(hue_pal()(3) __< ggplot default colors
rownames(group.name.colour.key) = as.integer(group.name.number.key$Group_number)



#DMR-detection
tiling.window.size=1000 #DEFAULT: 1000
tiling.step.size=1000    #DEFAULT: 1000
min.covered.bases= 0     #DEFAULT: 0


meth.change.cutoff = 25 #DEFAULT: 25. Abs. methylation % change between test and control
qvalue.cutoff = 0.05 #DEFAULT: 0.01 cutoff for qvalue of differential methylation statistic 
overdispersion.method = "MN" #DEFAULT: "none"
diff.meth.test = "Chisq" #DEFAULT: Chisq"
adjust.method = "fdr" #DEFAULT: SLIM


###################
#Reading data

methObj = prepare_methylkit_object(location=data.path, pattern.to.remove.from.names = "_CpG.txt$", assembly.name = assembly.name, sample.sheet = sample.sheet)


## Ordering the object based on sample sheet
methObj = reorganize(methObj, sample.ids=as.character(sample.info$Sample_id), treatment=sample.info$Group_number, chunk.size = 1e+06,save.db = FALSE)

filtered.methObj = preprocess_methylkit_object(methObj, normalize=do.normalisation,
                                               filter=do.filtering, 
                                               min.cov.as.reads = min.coverage.as.reads,
                                               max.cov.as.perc = max.coverage.as.perc)



##################
#DMR-detection

tiles=tileMethylCounts(filtered.methObj, win.size=tiling.window.size, step.size=tiling.step.size, cov.bases=min.covered.bases, mc.cores = available.cores)
tiles.united = methylKit::unite(tiles, destrand=FALSE, mc.cores = available.cores)

tiles.united = methylKit::select(tiles.united, !grepl("_|X|Y", tiles.united$chr)) 

meth.matrix = percMethylation(tiles.united)


#Forming unique name for each tile from tiles.united.object
regions = paste0(tiles.united$chr, ":", tiles.united$start,"-", tiles.united$end)
rownames(meth.matrix) = regions
rownames(tiles.united) = regions

meth.matrix = meth.matrix[, sample.info$Sample_id] #Ensuring we are still having the data in the same order



## Storing the full methylation matrix
meth.matrix=percMethylation(tiles.united)

if(save.DMRs){
  
  fname = paste0("FULL_methylation_matrix_with_window_size_", tiling.window.size, "_read.thr",
                 min.coverage.as.reads, ".tsv")
  fname = file.path(output.dir, fname)
  write.table(meth.matrix, fname, sep = "\t", quote = FALSE, row.names = TRUE)
  
  
}

if(remove.samples.from.tests){
  
  #removing sample from the sample(s) by using sample sheet
  ##### !!!!!!!!!!!!
  sample.info=sample.info[! sample.info$Sample_id %in% samples.to.remove,] #Removing metastasis sample
  ##### !!!!!!!!!!!!
  
  tiles.united = reorganize(tiles.united, sample.ids=as.character(sample.info$Sample_id), treatment=sample.info$Group_number,
                            chunk.size = 1e+06,save.db = FALSE)
  
}



comparisons.tiles = prepare_pairwise_group_comparisons_from_united_methObj(tiles.united, group.key=group.name.number.key)


### Storing DMRs and histograms
comparison.wise.DMRs = list()
meth.diff.histograms = list()

par(mfrow = c(length(comparisons.tiles),2))
for (i in 1:length(comparisons.tiles)){
  tiles.diff = calculateDiffMeth(comparisons.tiles[[i]], mc.cores=available.cores, overdispersion = overdispersion.method, test= diff.meth.test, adjust = adjust.method )
  av.meth.diff = round(median(tiles.diff$meth.diff), digits = 2)
  hist(tiles.diff$meth.diff, main = paste0(names(comparisons.tiles)[i]," methylation difference (%), m=",av.meth.diff ), plot = TRUE)
  hist(tiles.diff$qvalue, main = paste0(names(comparisons.tiles)[i]," q-values"), xlim = c(0, 1.0), plot = TRUE)
  
  comparison.wise.DMRs[[i]] = getMethylDiff(tiles.diff,difference = meth.change.cutoff,qvalue=qvalue.cutoff)
  
}
names(comparison.wise.DMRs)= names(comparisons.tiles)  


if(save.DMRs){
  
  for(comparison in 1:length(comparison.wise.DMRs)){
    fname = paste0(names(comparison.wise.DMRs)[comparison], "_methylKit_DMRs_q",qvalue.cutoff, "_methChange", meth.change.cutoff, ".tsv")
    fname = file.path(output.dir, fname)
    
    write.table(getData(comparison.wise.DMRs[[comparison]]), fname, sep = "\t", quote = FALSE, row.names = FALSE)
    
  }
}



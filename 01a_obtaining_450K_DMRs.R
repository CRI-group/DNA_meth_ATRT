######## Subsetting preprocessed Capper et al 2018 data:
#The original dataset contains 2801 samples of several brain cancer types. Only subset of them is relevant for this analysis. 

library(readxl)
library(minfi)

resource.dir = "your/path/"
result.dir = "your/path2/"
load("PATH/betas.ba.RData") # --> Capper et al result, contains two dataframes: "anno" for experinental info and "betas" for beta-values
clinical.anno = read_excel(file.path(resource.dir, "nature26000-s3.xlsx"), skip = 1)  # Capper Supplementary data

#merging clinical annotation with GEO-annotations:

#Extraction of array-ids from rownames,format  GSM2402854_5684819014_R03C02. Ensuring, that data is in same order
geo.accessions = unlist(lapply(rownames(betas), substr, start=1, stop=10))
#identical(geo.accessions, anno$geo_accession) #Should be true!
anno[["Sentrix ID (.idat)"]] = unlist(lapply(rownames(betas), substr, start=12, stop=28))

#-->taking only uniques (2795), and those found in our GEO annos
anno = anno[!(duplicated(anno$`Sentrix ID (.idat)`)),]
anno = anno[ anno$`Sentrix ID (.idat)` %in% clinical.anno$`Sentrix ID (.idat)`,]

full.anno = merge(x = anno, y = clinical.anno, by = "Sentrix ID (.idat)", all.x = TRUE)

targets = data.frame(full.anno[grep("^ATRT|^MB|^PLEX|CONTR", full.anno$title),])

#Taking only relevant columns
targets = targets[,c("title", "geo_accession", "material.ch1", "methylation.class.ch1", "WHO.Grade", "Age..years.at.operation.", "Sex", "Tumour.Location", "Array.processing.date...batch")]

#Removing INFLAM- samples
targets = targets[!grepl("INFLAM", targets$title),]


#extracting primary cancer groups and subclass in separate columns
targets$primary.group = gsub("\\,.*", "", targets$methylation.class.ch1)
targets$sample.details= gsub(".*, ", "", targets$methylation.class.ch1)
# <--- primary group (ATRT, MB PLEX or CONTR)

#Ordering the data based on primary cancer type:
targets$primary.group = factor(targets$primary.group, levels = c("ATRT","MB", "PLEX", "CONTR"))
targets = targets[order(targets$primary.group),]
rownames(targets) = targets$geo_accession

#Choosing the samples from the matrix (naming data with original geo_accessions first)
rownames(betas) = geo.accessions
bVals = t(as.matrix(betas[rownames(targets),])) #Subsequent analysis expects to have probeids as rows

###Storing the matrix for further use
grset <- makeGenomicRatioSetFromMatrix(bVals, pData=targets)
fname = paste("450k_bVals_as_GenomicRatioSet.RDS")
fname = file.path(result.dir, fname)
saveRDS(grset, file = fname)
rm(grset)



###########******************************************************************
#steps to illumina 450k liftover:
#Check  http://www.imsbio.co.jp/RGM/R_rdfile?f=IlluminaHumanMethylation450kprobe/man/IlluminaHumanMethylation450kprobe.Rd&d=R_BC
library(liftOver)
library(IlluminaHumanMethylation450kprobe
        
data(IlluminaHumanMethylation450kprobe) 

## 1) importing the chain-file downloaded from UCSC:
hg19ToHg38 =file.path(resource.dir, "hg19ToHg38.over.chain")
chain = import.chain(hg19ToHg38, exclude = "_") #Alignments for chromosomes matching the exclude pattern are not imported. 


# 2)converting hg19 probes (data-frame) into GRanges, chromosome-notation needs an update

probes = IlluminaHumanMethylation450kprobe
probes$chr = paste0("chr",probes$chr)
hg19.probes.Granges = makeGRangesFromDataFrame(probes, keep.extra.columns=TRUE,
                                               ignore.strand=FALSE,seqinfo=NULL,
                                               seqnames.field="chr",
                                               start.field="start",
                                               end.field="end",
                                               strand.field="strand")

# Sorting the probes by coordinates (for distance estimation)
hg19.probes.Granges = sortSeqlevels(hg19.probes.Granges)
hg19.probes.Granges = sort(hg19.probes.Granges)


## 3)Liftover operation. This originally outputs GRangesList,
##  which is turned into GRranges by unlisting.

hg38.probes.Granges <- liftOver(hg19.probes.Granges, chain)
hg38.probes.Granges = unlist(hg38.probes.Granges)  # 485748 ranges
genome(hg38.probes.Granges) = "hg38"

# Ensuring the probes are ordered by coordinates (for distance estimation)
hg38.probes.Granges = sortSeqlevels(hg38.probes.Granges)
hg38.probes.Granges = sort(hg38.probes.Granges)

# removing those that are not uniquely mapped
hg38.probes.Granges = hg38.probes.Granges[!duplicated(hg38.probes.Granges$Probe_ID),] #485383 ranges 
names(hg38.probes.Granges) = hg38.probes.Granges$Probe_ID



############******************************************************************
 #DMR-calling:

library(DMRcate)
library(GenomicRanges)


#Functions
sort_GRanges_by_coordinates = function(granges){
  granges <- sortSeqlevels(granges)
  granges <- sort(granges)
  
  return(granges)	
} 


make_DMR_GRanges = function(meth.matrix=NULL, matrix.type = c("M", "beta"), 
                            design.matrix = NULL, contrast.matrix = NULL,
                            contrast.name = NULL, liftover.probes = NULL,
                            available.cores = 10) {
  
  #preprocesssing for DMRcate, getting coordinates and statistics:
  myAnnotation <- cpg.annotate(object = meth.matrix, datatype = "array", what = matrix.type,
                               analysis.type = "differential", design = design.matrix,
                               contrasts = TRUE, cont.matrix = contrast.matrix,
                               coef = contrast.name, arraytype = "450K", fdr = 0.05)
  
  #Updating genomic coordinates for the probes:
  sign.probes = as.character(myAnnotation$ID)
  sign.probe.idx = chmatch(sign.probes, liftover.probes$Probe_ID)
  hg38.chr.for.sign.probes =  liftover.probes[sign.probe.idx, "seqnames" ] 
  hg38.pos.for.sign.probes =  liftover.probes[sign.probe.idx, "start" ]
  
  myAnnotation$CHR = hg38.chr.for.sign.probes
  myAnnotation$pos = hg38.pos.for.sign.probes
  
  #full result is large data object
  res <- dmrcate(myAnnotation, lambda=1000, C=2, mc.cores = available.cores, betacutoff = 0.05, min.cpgs = 2 )
  
  #Extracting only genomic ranges (with hg38 promoter annotations9)
  res.granges <- extractRanges(res, genome = "hg38")
  
  return(res.granges)
  
}


##DMRcate filtering
bVals <- rmSNPandCH(bVals, dist=2, mafcut=0.05) 
remaining.probes = nrow(bVals)

mVals =  log2(bVals/(1-bVals))


## inputs for DMR-analysis
cellType <- factor(targets$primary.group)

design <- model.matrix(~0+cellType, data=targets)
colnames(design) <- c(levels(cellType))


# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(ATRT-PLEX,
                            PLEX-MB,
                            ATRT-MB,
                            levels=design)


contrast.names=c("ATRT - PLEX",
                 "PLEX - MB",
                 "ATRT - MB")

hg38.probes = data.frame(hg38.probes.Granges)
hg38.probes$seqnames = droplevels(hg38.probes$seqnames)



##### A) SIMPLE DESIGN

contrast.DMRs.Granges = list()
for(i in 1:length(contrast.names)){
  
  contrast.DMRs.Granges[[i]] = make_DMR_GRanges(meth.matrix=mVals, matrix.type = "M", 
                                                design.matrix = design, contrast.matrix = contMatrix,
                                                contrast.name = contrast.names[i], liftover.probes = hg38.probes,
                                                available.cores = available.cores)
  
}
names(contrast.DMRs.Granges) = contrast.names

### B) INCLUDING SIMPLIFIED LOCATION DATA INTO MODEL:

located.targets = targets[!grepl("available", targets$Tumour.Location),]
located.targets$simplified.location = sub('^([^,;]+([,;][^,;]+)?).*', '\\1', located.targets$Tumour.Location)
located.targets$simplified.location= gsub(", ", "...", located.targets$simplified.location)
located.targets$simplified.location= gsub(" ", ".", located.targets$simplified.location)
located.targets$simplified.location= gsub("-", "_", located.targets$simplified.location)
print(table(located.targets$simplified.location))

# Adding location into the model
cellType <- factor(located.targets$primary.group)
location <- factor(located.targets$simplified.location)
design3 <- model.matrix(~0+cellType+location, data=located.targets)
colnames(design3) <- c(levels(cellType),levels(location)[-1]) 
rownames(design3) <- located.targets$geo_accession

# create a contrast matrix for specific comparisons
contMatrix3 <- makeContrasts(ATRT-PLEX,
                             PLEX-MB,
                             ATRT-MB,
                             ATRT-CONTR,
                             PLEX-CONTR,
                             MB-CONTR,
                             levels=design3)


located.mVals = mVals[,rownames(design3)]


#1) Calculating DMRs for cancer comparisons (="contrasts")
contrast.location.DMRs.Granges = list()
for(i in 1:length(contrast.names)){
  
  contrast.location.DMRs.Granges[[i]] = make_DMR_GRanges(meth.matrix=located.mVals, matrix.type = "M", 
                                                         design.matrix = design3, contrast.matrix = contMatrix3,
                                                         contrast.name = contrast.names[i], liftover.probes = hg38.probes,
                                                         available.cores = available.cores)
  
}
names(contrast.location.DMRs.Granges) = contrast.names



#2) Calculating DMRs for control comparisons separately

control.comparison.names=c("ATRT - CONTR",
                           "PLEX - CONTR",
                           "MB - CONTR")

control.location.DMRs.Granges = list()
for(i in 1:length(control.comparison.names)){
  
  control.location.DMRs.Granges[[i]] = make_DMR_GRanges(meth.matrix=located.mVals, matrix.type = "M", 
                                                        design.matrix = design3, contrast.matrix = contMatrix3,
                                                        contrast.name = control.comparison.names[i], liftover.probes = hg38.probes,
                                                        available.cores = available.cores)
  
}
names(control.location.DMRs.Granges) = control.comparison.names


#Filtering  contrast.location.DMRs.Granges with control.location.DMRs.Granges. 
# store results into contrast.location.DMRs.Granges.
intersected.Granges = list()
for(i in 1:length(contrast.names)) {
  
  contrast.name =  contrast.names[i]
  print(contrast.name)
  control.ids = names(control.location.DMRs.Granges)[grep(gsub(" - ", "|", contrast.name), names(control.location.DMRs.Granges))]
  print(control.ids)
  
  sort_GRanges_by_coordinates                                                    
  controlA = sort_GRanges_by_coordinates(control.location.DMRs.Granges[[control.ids[1]]])
  controlB = sort_GRanges_by_coordinates(control.location.DMRs.Granges[[control.ids[2]]])
  A.vs.B = sort_GRanges_by_coordinates(contrast.location.DMRs.Granges[[i]])
  
  #print(controlA)
  #print(controlB)
  #print(A.vs.B)
  intersected.Granges[[i]] = Reduce(subsetByOverlaps, list(A.vs.B, controlB, controlA ))
  intersected.Granges[[i]] = intersected.Granges[[i]][width(intersected.Granges[[i]]) >= 50 ]
  
}

contrast.location.DMRs.Granges = intersected.Granges
names(contrast.location.DMRs.Granges) = contrast.names

for( i in 1:length(contrast.location.DMRs.Granges)){
  
  fname = gsub(" - ", ".vs.", names(contrast.location.DMRs.Granges)[i])
  fname=paste0(fname,
               "_located_DMRcate_GRanges_filtered_by_controls.rds")
  fname = file.path(result.dir, fname)
  saveRDS(contrast.location.DMRs.Granges[[i]], fname)
}

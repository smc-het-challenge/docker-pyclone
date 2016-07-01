#!/usr/bin/env Rscript

rm(list = ls())
library(dplyr)
library(ccube)
library(colorspace)
library(mcclust)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
myColors <- gg_color_hue(10)

compute_mpear_label <- function(label_traces){
  ltmat <- as.matrix(label_traces)
  ltmat <- ltmat + 1
  psm <- comp.psm(ltmat)
  mpear <- maxpear(psm)
  mpear_label <- mpear$cl
  return(mpear_label)
}

args <- commandArgs(trailingOnly = TRUE)
#vcfFile <- as.character(args[1])
#batternbergFile <- as.character(args[2])
vcfFile <- "Tumour4/Tumour4.mutect.vcf"
batternbergFile <- "Tumour4/Tumour4.battenberg.txt"

ssm_file <- "ssm_data.txt"
cnv_file <- "cnv_data.txt"

shellCommandMutectSmcHet <- paste(
  "python create_ccfclust_inputs.py -v mutect_smchet",
  " -b ", batternbergFile,
  " -c ", 1,
  " --output-cnvs ", cnv_file,
  " --output-variants ", ssm_file,
  " ", vcfFile, sep = ""
)
system(shellCommandMutectSmcHet, intern = TRUE)
ssm <- read.delim(ssm_file,
                  stringsAsFactors = F)
cnv <- read.delim(cnv_file,
                  stringsAsFactors = F)

ssm$major_cn = 1
ssm$minor_cn = 1
ssm$cn_frac = 1
ssm$mu_r <- NULL
ssm$mu_v <- NULL
ssm$cn_ref <- NULL
ssm$cn_tot <- NULL

cnv <- na.omit(cnv)
cnv <- dplyr::filter(cnv, ssms!="" )

if (nrow(cnv)>0) {
  cnvtmp1 <- strsplit(as.character(cnv$ssms), ";")
  for (j in seq_len(nrow(cnv))) {
    if (length(cnvtmp1[[j]])==0) { next }
    cnvtmp1[[j]] = paste(cnvtmp1[[j]], cnv[j,]$frac, sep="," )
  }
  cnvtmp1 <- unlist(cnvtmp1)
  cnvtmp2 <- Reduce(
    rbind, strsplit(cnvtmp1, ",")
  )
  
  if (is.null(dim(cnvtmp2) )) {
    cnvtmp2 = as.data.frame(t(cnvtmp2), stringsAsFactors=F)
  } else {
    cnvtmp2 = as.data.frame(cnvtmp2, stringsAsFactors=F)
  }
  
  for (j in 2:ncol(cnvtmp2)) {
    cnvtmp2[,j] = as.numeric(cnvtmp2[,j])
  }
  
  ssm <- dplyr::left_join(ssm, cnvtmp2, by=c("id"="V1"))
  ssm$major_cn <- ssm$V3
  ssm$minor_cn <- ssm$V2
  ssm$cn_frac <- ssm$V4
  
  ssm$V2 <- NULL
  ssm$V3 <- NULL
  ssm$V4 <- NULL
  
  ssm[is.na(ssm[,5]), 5] = 1
  ssm[is.na(ssm[,6]), 6] = 1
  ssm[is.na(ssm[,7]), 7] = 1
}

allSsm <- dplyr::mutate(rowwise(ssm), chr =  strsplit(gene, split = "_")[[1]][1],
                        pos = strsplit(gene, split = "_")[[1]][2]) 

clonalCnFrac <- sum(allSsm$cn_frac==1)/nrow(allSsm)
ssm <- dplyr::filter(allSsm, cn_frac==1 & !chr %in% c("x", "y") )
problemSsm <- dplyr::filter(allSsm, cn_frac!=1 | chr %in% c("x", "y") )

# maxSnv <- 30000
# if (nrow(ssm) > maxSnv) {
#   ssm <- dplyr::sample_n(ssm, maxSnv)
# }

ssm$normal_cn = 2
ssm <- dplyr::rename(ssm, ref_counts=a, total_counts=d)
ssm <- dplyr::mutate(ssm, var_counts=total_counts-ref_counts, mutation_id = gene)
ssm$purity <- GetPurity(ssm)
cellularity <- unique(ssm$purity)
write.table(cellularity, file = "1A.txt", sep = "\t", row.names = F, 
            col.names = F, quote = F)


pycloneFolder <- paste0("purity_ccube_pyclone_v0.3/")


unlink(pycloneFolder, recursive = T, force = T)
dir.create(pycloneFolder, recursive = T)

pycloneData <- data.frame(mutation_id = ssm$id, ref_counts = ssm$ref_counts, 
                          var_counts = ssm$var_counts,
                          normal_cn = ssm$normal_cn, minor_cn =ssm$minor_cn, major_cn = ssm$major_cn )
write.table(pycloneData, file=paste0(pycloneFolder, '/pyclone_data.tsv'), quote=FALSE, sep='\t', row.names = F)

shellCommand <- paste0("/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/PyClone build_mutations_file ", 
                       pycloneFolder, '/pyclone_data.tsv ',
                       pycloneFolder, '/pyclone_mutations.yaml ', 
                       "--var_prior ", "parental_copy_number")

system(shellCommand)

yamlScript <- cat("num_iters: 100", "\n",
                  "base_measure_params:", "\n",
                  "  alpha: 1", "\n",
                  "  beta: 1", "\n",
                  "concentration:", "\n",
                  "  value: 1.0", "\n",
                  "  prior:", "\n",
                  "    shape: 1.0", "\n",
                  "    rate: 0.001", "\n",
                "density: pyclone_beta_binomial", "\n",
                "beta_binomial_precision_params:", "\n",
                "  value: 1000", "\n",
                "  prior:", "\n",
                "    shape: 1.0", "\n",
                "    rate: 0.0001", "\n",
                "  proposal:", "\n",
                "    precision: 0.01", "\n",
                paste0('working_dir: ', pycloneFolder), "\n",
                "trace_dir: trace", "\n",
                "samples:", "\n",
                paste0("  ", "sampleName", ":"), "\n",
                paste0("    mutations_file: ", "pyclone_mutations.yaml"), "\n",
                "    tumour_content:", "\n",
                paste0("      value: ", cellularity), "\n",
                "    error_rate: 0.001", 
                sep = "",
                file = paste0(pycloneFolder, "/pyclone_configure.yaml")
)

dir.create(pycloneFolder,"/trace")

shellCommand <- paste0("/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/PyClone analyse ", 
                       pycloneFolder, '/pyclone_configure.yaml ',
                       "--seed 1234")

system(shellCommand)

labelTrace <- read.delim(paste0(pycloneFolder, "/trace/labels.tsv.bz2"), stringsAsFactors = F)
burnIn <- 50
mpearLabels <- compute_mpear_label(labelTrace[-1:-burnIn,])
uniqLabels <- unique(mpearLabels)
#write.table(length(uniqLabels), file = "1B.txt", sep = "\t", row.names = F, col.names=F, quote = F)




## write results
codeName <- paste0("results/", "purity_ccube_pyclone_v0.3/")
resultsFolder <- paste0(icgc, codeName, sampleName)
dir.create(resultsFolder, recursive = T)

# load tsv
sampleTsvs <- paste0(pycloneFolder, "/pyclone_data.tsv")
allData <- read.delim(sampleTsvs, stringsAsFactors = F)
allData <- mutate(allData, vaf = var_counts/(ref_counts+var_counts))

# load trace and mpear label
traceFile <- dir(paste0(pycloneFolder, "/trace"), 
                 pattern = "cellular_frequencies", full.names = T)
paramsTrace <- read.delim(traceFile, stringsAsFactors = F, header = F)
id <- as.character(paramsTrace[1,])
burnIn <- 2000
tmp <- as.matrix( paramsTrace[-1:-(burnIn+1), ])
class(tmp) <- "numeric"
paramsEst <- colMeans(tmp)
names(paramsEst) <- NULL
traceData <- data.frame(mutation_id = id, ccf = paramsEst)

mpearFile <- paste0(pycloneFolder, "/pyclone_mpear.tsv")
if (file.exists(mpearFile)) {
  mpear <- read.delim(paste0(pycloneFolder, "/pyclone_mpear.tsv"), stringsAsFactors = F)
} else {
  library(mcclust)
  compute_mpear_label <- function(label_traces){
    ltmat <- as.matrix(label_traces)
    ltmat <- ltmat + 1
    psm <- comp.psm(ltmat)
    mpear <- maxpear(psm)
    mpear_label <- mpear$cl
    return(mpear_label)
  }
  
  labelTrace <- read.delim(paste0(pycloneFolder, "/trace/labels.tsv.bz2"), stringsAsFactors = F)
  burnIn <- 2000
  mpearLabels <- compute_mpear_label(labelTrace[-1:-burnIn,])
  mpear <- data.frame(mutation_id = colnames(labelTrace), cluster_id = mpearLabels)
  write.table(mpear, file = paste0(pycloneFolder, '/pyclone_mpear.tsv'), 
              row.names = F, sep = "\t", quote = F)
} 

allData <- left_join(allData, traceData, by="mutation_id")
allData <- left_join(allData, mpear, by="mutation_id")

tt <- table(allData$cluster_id)
clusterMean <- vector(mode = "numeric", length = length(tt))
clusterSd <- clusterMean
for (ii in seq_along(tt)) {
  clusterMean[ii] <- if (tt[ii]/sum(tt) > 0.01) {
    mean(filter(allData, cluster_id == as.integer(names(tt[ii])))$ccf)
  } else {NA}
  
  clusterSd[ii] <- if (tt[ii]/sum(tt) > 0.01) {
    sd(filter(allData, cluster_id == as.integer(names(tt[ii])))$ccf)
  } else {NA}
}

clusterDf <- data.frame(cluster_id=as.integer(names(tt)), 
                        average_ccf = clusterMean,
                        lower_95_ci = clusterMean - 2*clusterSd,
                        upper_95_ci = clusterMean + 2*clusterSd)
allData <- left_join(allData, clusterDf, by="cluster_id" )


id <- Reduce(rbind, strsplit(as.character(ssm$gene), "_", fixed = T), c())
mutAssign <- data.frame(mutation_id = ssm$id, chr = id[,1], pos = id[,2])

allData <- left_join(allData, mutAssign, by = "mutation_id") 
allData <- dplyr::filter(allData, !is.na(average_ccf))

# co-assignments matrix file
labelFile <- dir(paste0(pycloneFolder, "/trace/"), pattern = "labels", full.names = T)
labelTrace <- read.delim(labelFile, stringsAsFactors = F)
ltmat <- as.matrix(labelTrace[-1:-burnIn, allData$mutation_id])
ltmat <- ltmat + 1
psm <- comp.psm(ltmat)
fn <- paste0(resultsFolder, "/",
             sampleName, "_coassignment_probabilities.txt")
write.table(psm, file = fn, sep = "\t", row.names = F, quote = F)
shellCommand <- paste0("gzip -f ", fn)
system(shellCommand, intern = TRUE)

# index file
index <- allData[, c("chr", "pos")]
index$col <- seq_along(allData$mutation_id)
fn <- paste0(resultsFolder, "/",
             sampleName, "_index.txt")
write.table(index, file = fn, sep = "\t", row.names = F, quote = F)
shellCommand <- paste0("gzip -f ", fn)
system(shellCommand, intern = TRUE)

# cluster certainty file 
clusterCertainty <- subset(allData, 
                           select = c("chr", "pos", "cluster_id", 
                                      "average_ccf", "lower_95_ci", "upper_95_ci"))
clusterCertainty <- rename(clusterCertainty, most_likely_assignment = cluster_id)
clusterCertainty$most_likely_assignment <- 
  match(clusterCertainty$most_likely_assignment, 
        sort(unique(clusterCertainty$most_likely_assignment)))

tmp11 <- clusterCertainty[, c("chr", "pos", "most_likely_assignment")] 
tmp11 <- rename(tmp11, cluster = most_likely_assignment)
fn <- paste0(resultsFolder, "/", 
             sampleName, "_mutation_assignments.txt")
write.table(tmp11, file = fn, sep = "\t", row.names = F, quote = F)
shellCommand <- paste0("gzip -f ", fn)
system(shellCommand, intern = TRUE)

# multiplicity
mult <- allData[, c("chr", "pos")]
mult$tumour_copynumber <- allData$major_cn+allData$minor_cn
mult$multiplicity <- 1
fn <- paste0(resultsFolder, "/", 
             sampleName, "_multiplicity.txt")
write.table(mult, file = fn, sep = "\t", row.names = F, quote = F)
shellCommand <- paste0("gzip -f ", fn)
system(shellCommand, intern = TRUE)

# subclonal_structure file
tmp1 <- as.data.frame(table(clusterCertainty$most_likely_assignment), stringsAsFactors = F)
tmp2 <- as.data.frame(table(clusterCertainty$average_ccf), stringsAsFactors = F)
tmp <- left_join(tmp1, tmp2, by ="Freq")
tmp <- rename(tmp, cluster = Var1.x, n_ssms = Freq, proportion = Var1.y)
tmp <- mutate(tmp, proportion = as.numeric(proportion) * cellularity)
fn <- paste0(resultsFolder, "/", 
             sampleName, "_subclonal_structure.txt")
write.table(tmp, file = fn, sep = "\t", row.names = F, quote = F)
shellCommand <- paste0("gzip -f ", fn)
system(shellCommand, intern = TRUE)

# save sample summary
allData$cluster_id <- clusterCertainty$most_likely_assignment 
fn <- paste0(pycloneFolder, "/", sampleName, "_pyclone_results_table.csv") 
write.csv(allData, file = fn , row.names = F)

# graph summary
fn = paste0(resultsFolder, "/", 
            sampleName, "_pyclone_results_summary.pdf")
pdf(fn, width=8, height=8)
myColors <- gg_color_hue(n_distinct(allData$cluster_id))
par(mfrow=c(2,2))
plot(allData$ccf, allData$vaf, col = myColors[allData$cluster_id], 
     xlab = "cancer cell fraction", ylab = "variant allele frequecy", 
     main = "ccf vs vaf (colored by cluster memebership)")
hist(allData$ccf, density=20, breaks=20, prob=TRUE, 
     main = "ccf histogram",
     xlab = "cancer cell fraction")
clusterSize <- table(allData$average_ccf)/nrow(allData)
names(clusterSize) <- as.character(format(round(as.numeric(names(clusterSize)), 2), nsmall = 2))
tmp1 <- as.data.frame(table(allData$average_ccf))
tmp2 <- as.data.frame(table(allData$cluster_id))
tmp3 <- left_join(tmp1, tmp2, by ="Freq")
barplot(clusterSize, las = 2, col = myColors[tmp3$Var1.y], xlab = "cluster mean", ylab="mutation proportions", 
        main = "cluster sizes")
dev.off()
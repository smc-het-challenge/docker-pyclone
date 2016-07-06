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

Assign <- function(x, centers, s) {
  
  n <- length(x)
  
  k <- length(centers)
  
  logRho <- array(0, dim= c(n ,k))
  
  for (ii in 1:k) {
    logRho[,ii] = ccube::bsxfun.se("-", -(x-centers[ii])^2/(2*s[ii]), log(s[ii]))
  }
  
  if (n==k) {
    logR <- ccube::bsxfun.se("-", logRho, ccube:::logsumexp(logRho, 1), expandByRow = F)	# 10.49
  } else {
    logR <- ccube::bsxfun.se("-", logRho, ccube:::logsumexp(logRho, 1))	# 10.49
  }
  
  R <- exp(logR)
  
}  
args <- commandArgs(trailingOnly = TRUE)
vcfFile <- as.character(args[1])
batternbergFile <- as.character(args[2])
purityFile <- as.character(args[3])
numMCMC <- as.integer(args[4])
burnIn <- as.integer(args[5])
maxSnv <- as.integer(args[6]) 

ssm_file <- "ssm_data.txt"
cnv_file <- "cnv_data.txt"

shellCommandMutectSmcHet <- paste(
  "create_ccfclust_inputs.py -v mutect_smchet",
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
if (nrow(problemSsm)>0) {
  problemSsmFlag <- T
} else {
  problemSsmFlag <- F
}


if (nrow(ssm) > maxSnv) {
  tmp_ssm <- dplyr::sample_n(ssm, maxSnv)
  holdOutData <- dplyr::anti_join(ssm, tmp_ssm)
  ssm <- tmp_ssm
  holdOutDataFlag <- T
} else {
  holdOutDataFlag <- F
}

ssm$normal_cn = 2
ssm <- dplyr::rename(ssm, ref_counts=a, total_counts=d)
ssm <- dplyr::mutate(ssm, var_counts=total_counts-ref_counts, mutation_id = gene)
cellularity <-read.delim(purityFile, stringsAsFactors=FALSE)$cellularity
ssm$purity <- cellularity
write.table(cellularity, file = "1A.txt", sep = "\t", row.names = F, 
            col.names = F, quote = F)

pycloneFolder <- paste0("pyclone_v0.3/")
unlink(pycloneFolder, recursive = T, force = T)
dir.create(pycloneFolder, recursive = T)
pycloneData <- data.frame(mutation_id = ssm$id, ref_counts = ssm$ref_counts, 
                          var_counts = ssm$var_counts,
                          normal_cn = ssm$normal_cn, minor_cn =ssm$minor_cn, major_cn = ssm$major_cn )
write.table(pycloneData, file=paste0(pycloneFolder, '/pyclone_data.tsv'), quote=FALSE, sep='\t', row.names = F)
shellCommand <- paste0("PyClone build_mutations_file ", 
                       pycloneFolder, '/pyclone_data.tsv ',
                       pycloneFolder, '/pyclone_mutations.yaml ', 
                       "--var_prior ", "parental_copy_number")

system(shellCommand)

yamlScript <- cat(paste0("num_iters: ", numMCMC), "\n",
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

shellCommand <- paste0("PyClone analyse ", 
                       pycloneFolder, '/pyclone_configure.yaml ',
                       "--seed 1234")

system(shellCommand)

labelTrace <- read.delim(paste0(pycloneFolder, "/trace/labels.tsv.bz2"), stringsAsFactors = F)
mpearLabels <- compute_mpear_label(labelTrace[-1:-burnIn,])
mpear <- data.frame(mutation_id = colnames(labelTrace), cluster_id = mpearLabels)

# load tsv
sampleTsvs <- paste0(pycloneFolder, "/pyclone_data.tsv")
allData <- read.delim(sampleTsvs, stringsAsFactors = F)
allData <- mutate(allData, vaf = var_counts/(ref_counts+var_counts))

# load trace and mpear label
traceFile <- dir(paste0(pycloneFolder, "/trace"), 
                 pattern = "cellular_frequencies", full.names = T)
paramsTrace <- read.delim(traceFile, stringsAsFactors = F, header = F)
id <- as.character(paramsTrace[1,])
tmp <- as.matrix( paramsTrace[-1:-(burnIn+1), ])
class(tmp) <- "numeric"
paramsEst <- colMeans(tmp)
names(paramsEst) <- NULL
traceData <- data.frame(mutation_id = id, ccf = paramsEst)
 
allData <- left_join(allData, traceData, by="mutation_id")
allData <- left_join(allData, mpear, by="mutation_id")

tt <- table(allData$cluster_id)
clusterMean <- vector(mode = "numeric", length = length(tt))
clusterSd <- clusterMean
for (ii in seq_along(tt)) {
  clusterMean[ii] <- if (tt[ii]/sum(tt) > 0.01 & 
                         nrow(filter(allData, cluster_id == as.integer(names(tt[ii])))) > 1
                         ) {
    mean(filter(allData, cluster_id == as.integer(names(tt[ii])))$ccf)
  } else {NA}
  
  clusterSd[ii] <- if (tt[ii]/sum(tt) > 0.01 & 
                       nrow(filter(allData, cluster_id == as.integer(names(tt[ii])))) > 1
                       ) {
   
    if (sd(filter(allData, cluster_id == as.integer(names(tt[ii])))$ccf) ==0) {
      1e-20
    } else { sd(filter(allData, cluster_id == as.integer(names(tt[ii])))$ccf)}
  } else {NA}
}

clusterDf <- data.frame(cluster_id=as.integer(names(tt)), 
                        average_ccf = clusterMean,
                        lower_95_ci = clusterMean - 2*clusterSd,
                        upper_95_ci = clusterMean + 2*clusterSd)
allData <- left_join(allData, clusterDf, by="cluster_id" )
clusterDf <- dplyr::filter(clusterDf, !is.na(average_ccf))
write.table(nrow(clusterDf), file = "1B.txt", sep = "\t", row.names = F, col.names=F, quote = F)

id <- Reduce(rbind, strsplit(as.character(ssm$gene), "_", fixed = T), c())
mutAssign <- data.frame(mutation_id = ssm$id, chr = id[,1], pos = id[,2], 
                        stringsAsFactors = F)
allData <- left_join(allData, mutAssign, by = "mutation_id") 

## Reassign low support data
lowSupportData <- dplyr::filter(allData, is.na(average_ccf))
allData <- dplyr::filter(allData, !is.na(average_ccf))
ssm <- allData
lowSupportDataFlag <- F
if (nrow(lowSupportData) > 0) {
  lowSupportDataFlag <- T
  lowSupportR <- Assign(lowSupportData$ccf, clusterDf$average_ccf, 
                        ((clusterDf$upper_95_ci - clusterDf$average_ccf)/2)^2)
  lowSupportLabel <- apply(lowSupportR, 1, which.max)
  lowSupportData$cluster_id <- clusterDf$cluster_id[lowSupportLabel]
  lowSupportData$average_ccf <- clusterDf$average_ccf[lowSupportLabel]
  lowSupportData$lower_95_ci <- clusterDf$lower_95_ci[lowSupportLabel]
  lowSupportData$upper_95_ci <- clusterDf$upper_95_ci[lowSupportLabel]
  allData <- rbind(allData, lowSupportData)
  rm(lowSupportR)
}

## Assign hold-out data
if(holdOutDataFlag) {
  holdOutData$normal_cn <- 2
  holdOutData$purity <- cellularity
  holdOutData <- rename(holdOutData, mutation_id = id)
  holdOutData <- rename(holdOutData, ref_counts = a, total_counts = d)
  holdOutData <- mutate(rowwise(holdOutData), 
                        var_counts = total_counts - ref_counts, 
                        vaf = var_counts/total_counts,
                        ccf1 = ccube::MapVaf2CcfPyClone(vaf, 
                                                        purity, 
                                                        normal_cn, 
                                                        major_cn+minor_cn, 
                                                        major_cn+minor_cn,
                                                        major_cn, 
                                                        lower = 0,
                                                        upper = 2), 
                        ccf2 = ccube::MapVaf2CcfPyClone(vaf, 
                                                        purity, 
                                                        normal_cn, 
                                                        major_cn+minor_cn, 
                                                        major_cn+minor_cn,
                                                        minor_cn,
                                                        lower = 0,
                                                        upper = 2),
                        ccf3 = ccube::MapVaf2CcfPyClone(vaf, 
                                                        purity, 
                                                        normal_cn, 
                                                        major_cn+minor_cn, 
                                                        major_cn+minor_cn,
                                                        1,
                                                        lower = 0,
                                                        upper = 2), 
                        ccf = mean(unique(c(ccf1, ccf2, ccf3)), na.rm=TRUE))
  
  holdOutData <- mutate(holdOutData, ccf = if (is.na(ccf)) {vaf} else {ccf} )
  holdOutR <- Assign(holdOutData$ccf, clusterDf$average_ccf, 
                     ((clusterDf$upper_95_ci - clusterDf$average_ccf)/2)^2) 
  holdOutLabel <- apply(holdOutR, 1, which.max)
  holdOutData$cluster_id <- clusterDf$cluster_id[holdOutLabel]
  holdOutData$average_ccf <- clusterDf$average_ccf[holdOutLabel]
  holdOutData$lower_95_ci <- clusterDf$lower_95_ci[holdOutLabel]
  holdOutData$upper_95_ci <- clusterDf$upper_95_ci[holdOutLabel]
  holdOutData$cn_frac <- NULL
  holdOutData$purity <- NULL
  holdOutData$ccf1 <- NULL
  holdOutData$ccf2 <- NULL
  holdOutData$ccf3 <- NULL
  holdOutData$gene <- NULL
  holdOutData$total_counts <- NULL
  allData <- rbind(allData, holdOutData)
  rm(holdOutR)
}

## Post assign problem SSMs-- temporal solution
MapVaf2CcfTest <- function(x, t, cv, bv, frac,
                           epi = 1e-3, constraint = T,
                           lower = -Inf, upper = Inf) {
  if(bv==0) {
    return(0)
  }
  
  cn2 = 2
  zz = (1-t)*2 + t*(frac*cv + (1-frac)*2 )
  
  
  ccf <- ( x * zz - t*(1-frac)*2*epi - (1-t) *2*epi - t*frac*cv*epi  ) /
    ( t*frac*( bv *(1-epi) - cv*epi ) )
  
  if (constraint) {
    if (is.na(ccf)) {
      return(as.numeric(NA))
    } else if (ccf < 0.9 && bv > 1) {
      return(as.numeric(NA))
    } else if (ccf < upper && ccf > lower) {
      return(ccf)
    } else {
      return(as.numeric(NA))
    }
  } else {
    return(ccf)
  }
}

if (problemSsmFlag) {
  problemSsm$normal_cn <- 2
  problemSsm$purity <- cellularity
  problemSsm <- rename(problemSsm, mutation_id = id)
  problemSsm <- rename(problemSsm, ref_counts = a, total_counts = d)
  problemSsm <- mutate(rowwise(problemSsm), 
                       var_counts = total_counts - ref_counts, 
                       vaf = var_counts/total_counts,
                       ccf1 = MapVaf2CcfTest(var_counts/total_counts, 
                                             purity, 
                                             major_cn+minor_cn, 
                                             major_cn, cn_frac, constraint = F), 
                       ccf2 = MapVaf2CcfTest(var_counts/total_counts, 
                                             purity, 
                                             major_cn+minor_cn, 
                                             minor_cn, cn_frac, constraint = F),
                       ccf3 = MapVaf2CcfTest(var_counts/total_counts, 
                                             purity, 
                                             major_cn+minor_cn, 
                                             1, cn_frac, constraint = F), 
                       ccf = mean( unique( c(ccf1, ccf2, ccf3)) ) )
  
  problemSsm$purity <- NULL
  problemSsm$ccf1 <- NULL
  problemSsm$ccf2 <- NULL
  problemSsm$ccf3 <- NULL
  problemSsm$gene <- NULL
  problemSsm$total_counts <- NULL
  problemSsm$cn_frac <- NULL
  problemSsmR <- Assign(problemSsm$ccf, clusterDf$average_ccf, 
                        ((clusterDf$upper_95_ci - clusterDf$average_ccf)/2)^2)
  problemSsmLabel <- apply(problemSsmR, 1, which.max)
  problemSsm$cluster_id <- clusterDf$cluster_id[problemSsmLabel]
  problemSsm$average_ccf <- clusterDf$average_ccf[problemSsmLabel]
  problemSsm$lower_95_ci <- clusterDf$lower_95_ci[problemSsmLabel]
  problemSsm$upper_95_ci <- clusterDf$upper_95_ci[problemSsmLabel]
  
  allData <- rbind(allData, problemSsm)
  
  rm(problemSsmR)
}


## co-assignments matrix file
allData <- allData[order(as.numeric(gsub("[^\\d]+", "", allData$mutation_id, perl=TRUE))), ]
allR <- Assign(allData$ccf, clusterDf$average_ccf, 
               ((clusterDf$upper_95_ci - clusterDf$average_ccf)/2)^2)
allRR <- Matrix::tcrossprod(allR)
colnames(allRR) <- allData$mutation_id
rownames(allRR) <- allData$mutation_id
rm(allR)

ltmat <- as.matrix(labelTrace[-1:-burnIn, ssm$mutation_id])
ltmat <- ltmat + 1
psm <- comp.psm(ltmat)
colnames(psm) <- ssm$mutation_id
rownames(psm) <- ssm$mutation_id
allRR[ssm$mutation_id, ssm$mutation_id] <- psm
diag(allRR) <- 1
write.table(allRR, file = "2B.txt", sep = "\t", row.names = F, col.names = F, quote = F)


# cluster certainty file 
clusterCertainty <- subset(allData, 
                           select = c("chr", "pos", "cluster_id", 
                                      "average_ccf", "lower_95_ci", "upper_95_ci"))
clusterCertainty <- rename(clusterCertainty, most_likely_assignment = cluster_id)

# subclonal_structure file
tmp1 <- as.data.frame(table(clusterCertainty$most_likely_assignment), stringsAsFactors = F)
tmp1 <- mutate(tmp1, Var1 = as.integer(Var1))
tmp <- left_join(tmp1, clusterDf, by = c("Var1"="cluster_id"))
tmp <- rename(tmp, cluster = Var1, n_ssms = Freq, proportion = average_ccf)
tmp <- mutate(tmp, proportion = as.numeric(proportion) * cellularity)
tmp <- mutate(tmp, cluster = seq_along(cluster))
tmp$lower_95_ci <- NULL
tmp$upper_95_ci <- NULL
write.table(tmp, file = "1C.txt", sep = "\t", row.names = F, col.names = F, quote = F)


clusterCertainty$most_likely_assignment <- 
  match(clusterCertainty$most_likely_assignment, 
        sort(unique(clusterCertainty$most_likely_assignment)))
tmp11 <- clusterCertainty[, c("chr", "pos", "most_likely_assignment")] 
tmp11 <- rename(tmp11, cluster = most_likely_assignment)
write.table(tmp11$cluster, file = "2A.txt", sep = "\t", row.names = F, col.names = F, quote = F)


# graph summary
fn = "clonal_results_summary.pdf"
pdf(fn, width=8, height=8)
myColors <- gg_color_hue(max(unique(ssm$cluster_id)))
par(mfrow=c(2,2))
plot(ssm$ccf, ssm$vaf, col = myColors[ssm$cluster_id], 
     xlab = "cancer cell fraction", ylab = "variant allele frequecy", 
     main = "ccf vs vaf (colored by cluster memebership)")
hist(ssm$ccf, density=20, breaks=20, prob=TRUE, 
     main = "ccf histogram",
     xlab = "cancer cell fraction")
clusterSize <- table(ssm$average_ccf)/nrow(ssm)
names(clusterSize) <- as.character(format(round(as.numeric(names(clusterSize)), 2), nsmall = 2))
tmp1 <- as.data.frame(table(ssm$average_ccf), stringsAsFactors = F)
tmp2 <- as.data.frame(table(ssm$cluster_id), stringsAsFactors = F)
tmp3 <- left_join(tmp1, tmp2, by ="Freq")
barplot(clusterSize, las = 2, col = myColors[as.integer(tmp3$Var1.y)], xlab = "cluster mean", ylab="mutation proportions", 
        main = "cluster sizes")
dev.off()
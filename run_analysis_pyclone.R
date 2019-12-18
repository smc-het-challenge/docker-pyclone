#!/usr/bin/env Rscript

rm(list = ls())
library(dplyr)
library(ccube)
library(mcclust)
options(stringsAsFactors = F)

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
# Parse vcf file
vcfParserPath <- "create_ccfclust_inputs.py"
shellCommandMutectSmcHet <- paste(
  vcfParserPath,
  " -v mutect_smchet",
  " -c ", 1,
  " --output-variants ", ssm_file,
  " ", vcfFile, sep = ""
)

system(shellCommandMutectSmcHet, intern = TRUE)

ssm <- read.delim(ssm_file,
                  stringsAsFactors = F)


# Parse Battenberg CNA data
cna <- read.delim(batternbergFile, stringsAsFactors = F)
ssm <- ccube::ParseSnvCnaBattenberg(ssm, cna) 

cellularity <-read.delim(purityFile, stringsAsFactors=FALSE)$cellularity
ssm$purity <- cellularity
ssm$vaf = ssm$var_counts/ssm$total_counts
ssm <- ccube:::CheckAndPrepareCcubeInupts(ssm)

ssm <- dplyr::mutate(rowwise(ssm), 
                        chr =  strsplit(mutation_id, split = "_")[[1]][1],
                        pos = strsplit(mutation_id, split = "_")[[1]][2]) 

# 1st filter: major_cn == 0
falsePositiveSsmIDs <- filter(ssm, major_cn == 0 )$id
ssm <- filter(ssm, major_cn > 0)


# 2nd filter: fp_qval, rough_ccf1, rough_ccf0, total_counts, var_counts 
ssm<- mutate(rowwise(ssm), 
             fp_pval = binom.test(var_counts, total_counts, 5e-2, alternative = "great")$p.value, 
             fp_qval = p.adjust(fp_pval, method = "hochberg")
)

ssm<- mutate(rowwise(ssm), 
             rough_ccf1 =  
               MapVaf2CcfLinear( vaf, purity, normal_cn, total_cn, 1) 
)

ssm<- mutate(rowwise(ssm), 
             rough_ccf0 =  
               MapVaf2CcfLinear( vaf, purity, normal_cn, total_cn, 0) 
)


falsePositiveSsmIDs <- c(falsePositiveSsmIDs, filter(ssm, fp_qval > 0.05  &
                                           rough_ccf1 < 0.2 &
                                           rough_ccf0 > -90
                                         )$id)



if (length(falsePositiveSsmIDs)>0) {
  ssm = filter(ssm, !id %in% falsePositiveSsmIDs)
  FalsePositiveSsmFlag <- T
} else {
  FalsePositiveSsmFlag <- F
}

subClonalSsmIDs <- filter(ssm, subclonal_cn)$id
if (length(subClonalSsmIDs)>0) {
  problemSsm = filter(ssm, subclonal_cn)
  ssm = filter(ssm, !id %in% subClonalSsmIDs)
  subClonalSsmFlag <- T
} else {
  subClonalSsmFlag <- F
}


if (nrow(ssm) > maxSnv) {
  downSampledIDs <- dplyr::sample_n(ssm, maxSnv)$id
  holdOutData = filter(ssm, !id %in% downSampledIDs)
  ssm = filter(ssm, id %in% downSampledIDs) 
  holdOutDataFlag <- T
} else {
  holdOutDataFlag <- F
}



pycloneFolder <- paste0("pyclone_v0.3/")
unlink(pycloneFolder, recursive = T, force = T)
dir.create(pycloneFolder, recursive = T)
pycloneData <- data.frame(mutation_id = ssm$id, ref_counts = ssm$ref_counts, 
                          var_counts = ssm$var_counts,
                          normal_cn = ssm$normal_cn, minor_cn =ssm$minor_cn, major_cn = ssm$major_cn )
rm(ssm)
write.table(pycloneData, file=paste0(pycloneFolder, '/pyclone_data.tsv'), quote=FALSE, sep='\t', row.names = F)
pycloneDataIDs <- pycloneData$mutation_id
rm(pycloneData)
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

dir.create(paste0(pycloneFolder,"trace"))

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
rm(paramsTrace)
class(tmp) <- "numeric"
paramsEst <- colMeans(tmp)
rm(tmp)

names(paramsEst) <- NULL
traceData <- data.frame(mutation_id = id, ccf = paramsEst)
 
allData <- left_join(allData, traceData, by="mutation_id")
rm(traceData)
allData <- left_join(allData, mpear, by="mutation_id")
rm(mpear)

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
allData <- dplyr::left_join(allData, clusterDf, by="cluster_id" )
clusterDf <- dplyr::filter(clusterDf, !is.na(average_ccf))
write.table(nrow(clusterDf), file = "1B.txt", sep = "\t", row.names = F, col.names=F, quote = F)

allData <- dplyr::left_join(allData, 
                            data.frame(mutation_id = pycloneDataIDs), 
                            by = "mutation_id") 

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
  holdOutData <- GetCcf(holdOutData, use = "use_base")
  holdOutData <- dplyr::rename(holdOutData, ccf = rough_ccf)
  holdOutR <- Assign(holdOutData$ccf, clusterDf$average_ccf, 
                     ((clusterDf$upper_95_ci - clusterDf$average_ccf)/2)^2) 
  holdOutLabel <- apply(holdOutR, 1, which.max)
  holdOutData$cluster_id <- clusterDf$cluster_id[holdOutLabel]
  holdOutData$average_ccf <- clusterDf$average_ccf[holdOutLabel]
  holdOutData$lower_95_ci <- clusterDf$lower_95_ci[holdOutLabel]
  holdOutData$upper_95_ci <- clusterDf$upper_95_ci[holdOutLabel]
  holdOutData$mutation_id <- holdOutData$id
  holdOutData <- dplyr::select(holdOutData, colnames(allData))
  allData <- rbind(allData, holdOutData)
  rm(holdOutR, holdOutData)
}

## Post assign subclonal SSMs-- temporal solution
if (subClonalSsmFlag) {
  problemSsm <- mutate(rowwise(problemSsm), 
                       ccf1 = MapVaf2CcfLinear(vaf, 
                                             purity, 
                                             normal_cn,
                                             total_cn, 
                                             mult = major_cn), 
                       ccf1 = if (ccf1 <0) {0} else {ccf1}, 
                       ccf2 = MapVaf2CcfLinear(vaf, 
                                             purity, 
                                             normal_cn,
                                             total_cn, 
                                             minor_cn),
                       ccf2 = if (ccf2 <0) {0} else {ccf2}, 
                       ccf3 = MapVaf2CcfLinear(vaf, 
                                             purity, 
                                             normal_cn,
                                             total_cn, 
                                             frac_cn_sub1), 
                       ccf3 = if (ccf3 <0) {0} else {ccf3}, 
                       ccf4 = MapVaf2CcfLinear(vaf, 
                                             purity, 
                                             normal_cn,
                                             total_cn, 
                                             frac_cn_sub2), 
                       ccf4 = if (ccf4 <0) {0} else {ccf4}, 
                       ccf5 = MapVaf2CcfLinear(vaf, 
                                             purity, 
                                             normal_cn,
                                             total_cn, 
                                             frac_cn_sub1+frac_cn_sub2), 
                       ccf5 = if (ccf5 <0) {0} else {ccf5}, 
                       ccf = mean( unique( c(ccf1, ccf2, ccf3, ccf4, ccf5)) ) )
  
  problemSsmR <- Assign(problemSsm$ccf, clusterDf$average_ccf, 
                        ((clusterDf$upper_95_ci - clusterDf$average_ccf)/2)^2)
  problemSsmLabel <- apply(problemSsmR, 1, which.max)
  problemSsm$cluster_id <- clusterDf$cluster_id[problemSsmLabel]
  problemSsm$average_ccf <- clusterDf$average_ccf[problemSsmLabel]
  problemSsm$lower_95_ci <- clusterDf$lower_95_ci[problemSsmLabel]
  problemSsm$upper_95_ci <- clusterDf$upper_95_ci[problemSsmLabel]
  problemSsm <- dplyr::select(problemSsm, colnames(allData))
  allData <- rbind(allData, problemSsm)
  
  rm(problemSsmR, problemSsm)
}


## co-assignments matrix file
allData <- allData[order(as.numeric(gsub("[^\\d]+", "", allData$mutation_id, perl=TRUE))), ]
allR <- Assign(allData$ccf, clusterDf$average_ccf, 
               ((clusterDf$upper_95_ci - clusterDf$average_ccf)/2)^2)

if (FalsePositiveSsmFlag) {
  allR <- cbind(allR, matrix(0, nrow = nrow(allR), 1 ))
  fpR <- cbind(matrix(0, nrow = length(falsePositiveSsmIDs), ncol = nrow(clusterDf) ), 
               matrix(1, nrow = length(falsePositiveSsmIDs), ncol = 1 ) )
  allR <- rbind(allR, fpR)
  
  ssmPyClone = data.frame(id = c(allData$mutation_id, falsePositiveSsmIDs) )
  ssmPyClone = cbind(ssmPyClone,  allR)
} else {
  ssmPyClone = data.frame(id = allData$mutation_id)
  ssmPyClone = cbind(ssmPyClone,  allR)
}
rm(allR, fpR)

ssmPyClone <- ssmPyClone[order(as.numeric(gsub("[^\\d]+", "", ssmPyClone$id, perl=TRUE))), ]
rownames(ssmPyClone)<-NULL
allIDs <- ssmPyClone$id

ssmPyClone$id <- NULL
allRR <- Matrix::tcrossprod(as.matrix(ssmPyClone))
rm(ssmPyClone)
colnames(allRR) <- allIDs
rownames(allRR) <- allIDs


ltmat <- as.matrix(labelTrace[-1:-burnIn, ssm$mutation_id])
rm(labelTrace)
ltmat <- ltmat + 1
psm <- comp.psm(ltmat)
rm(ltmat)
colnames(psm) <- ssm$mutation_id
rownames(psm) <- ssm$mutation_id
allRR[ssm$mutation_id, ssm$mutation_id] <- psm
rm(psm)
diag(allRR) <- 1

write.table(allRR, file = "2B.txt", sep = "\t", row.names = F, col.names = F, quote = F)
rm(allRR)

# cluster certainty file 
clusterCertainty <- subset(allData, 
                           select = c("mutation_id", "cluster_id", 
                                      "average_ccf"))
rm(allData)

if (FalsePositiveSsmFlag) {
  clusterCertainty <- rbind(clusterCertainty, data.frame(mutation_id = falsePositiveSsmIDs, 
                                                         cluster_id = max(clusterCertainty$cluster_id)+1,
                                                         average_ccf = 0)  )
  clusterDf <- rbind(clusterDf, data.frame(cluster_id = max(clusterDf$cluster_id)+1,
                                            average_ccf = 0, lower_95_ci = 0, upper_95_ci = 0 ) )
}

# subclonal_structure file
tmp1 <- as.data.frame(table(clusterCertainty$cluster_id), stringsAsFactors = F)
tmp1 <- mutate(tmp1, Var1 = as.integer(Var1))
tmp <- left_join(tmp1, clusterDf, by = c("Var1"="cluster_id"))
tmp <- rename(tmp, cluster = Var1, n_ssms = Freq, proportion = average_ccf)
tmp <- mutate(tmp, proportion = as.numeric(proportion) * cellularity)
tmp <- mutate(tmp, cluster = seq_along(cluster))
tmp$lower_95_ci<-NULL
tmp$upper_95_ci<-NULL
write.table(tmp, file = "1C.txt", sep = "\t", row.names = F, col.names = F, quote = F)


clusterCertainty$cluster_id <- match(clusterCertainty$cluster_id, 
                                     sort(unique(clusterCertainty$cluster_id)))
write.table(clusterCertainty$cluster_id, file = "2A.txt", sep = "\t", row.names = F, col.names = F, quote = F)

#' bsxfun {pracma} with single expansion (real Matlab style)
#' @param func the function used by bsxfun
#' @param x a matrix
#' @param y a vector need to be expanded
#' @param  expandByRow applies only when x is a square matrix
#' @return value of func
#' @export
bsxfun.se <- function(func, x, y, expandByRow=TRUE) {

  if(length(y) == 1) return(pracma::arrayfun(func, x, y)) else
    stopifnot(nrow(x) == length(y) || ncol(x) == length(y))

  expandCol <- nrow(x) == length(y)
  expandRow <- ncol(x) == length(y)

  if(expandCol & expandRow & expandByRow) expandCol <- FALSE
  if(expandCol & expandRow & !expandByRow) expandRow <- FALSE

  # repeat row (if dim2expand = 1, then length(y) = ncol(x))
  if(expandRow) y.repmat <- matrix(rep(as.numeric(y), each=nrow(x)), nrow=nrow(x))

  # repeat col (if dim2expand = 2, then length(y) = nrow(x))
  if(expandCol) y.repmat <- matrix(rep(as.numeric(y), ncol(x)), ncol=ncol(x))

  pracma::bsxfun(func, x, y.repmat)
}

#' Compute log(sum(exp(x),dim)) while avoiding numerical underflow
#' @param x a matrix
#' @param margin used for apply
#' @return log(sum(exp(x),dim)) a matrix of the sample size as x
logsumexp <- function(x, margin=1) {

  if ( ! is.matrix(x) ) {
    x <- as.matrix(x)
  }

  # subtract the largest in each column
  y <- apply(x, margin, max)

  if (nrow(x) == ncol(x)) {
    x <- bsxfun.se("-", x, y, expandByRow = F)
  } else {
    x <- bsxfun.se("-", x, y)
  }

  s <- y + log(apply(exp(x), margin, sum))

  i <- which(!is.finite(s))

  if(length(i) > 0) s[i] <- y[i]

  s
}


#' Convert ccf to vaf
#' @param x ccf
#' @param t purity
#' @param cn normal total copy number
#' @param cr total copy number in reference population
#' @param cv total copy number in variant population
#' @param bv mutation multplicity
#' @param epi sequencing error
#' @return vaf
#' @export
cp2ap <- function(x, t, cn, cr, cv, bv, epi = 1e-3){
  p1 <- (1-t)*cn
  p2 <- t*(1-x)*cr
  p3 <- t*x*cv
  p <- p1+p2+p3
  u1 <- epi
  u2 <- epi
  u3 <- (bv/cv)*(1-epi)
  y <- p1*u1+p2*u2+p3*u3
  return(y/p)
}

#' Convert vaf to ccf using pyclone style mapping
#' @param x vaf
#' @param t purity
#' @param cn normal total copy number
#' @param cr total copy number in reference population
#' @param cv total copy number in variant population
#' @param bv mutation multplicity
#' @param epi sequencing error
#' @param constraint whether constriant the transformation
#' @param lower lower bound of the constriant
#' @param upper upper bound of the constriant
#' @return ccf
#' @export
MapVaf2CcfPyClone <- function(x, t, cn, cr, cv, bv,
                              epi = 1e-3, constraint = T,
                              lower = -Inf, upper = Inf) {
  un <- epi
  ur <- epi
  if (bv == cv) {
    uv <-  1-epi
  } else if ( bv == 0) {
    uv <- epi
  } else {
    uv <- (bv/cv)*(1-epi)
  }

  ccf <- ( (1-t) * cn * (un - x) + t * cr *(ur - x) ) /
    ( t * cr * (ur - x) - t * cv * (uv - x) )

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

#' Convert vaf to ccf using ccube linear mapping
#' @param vaf vaf
#' @param purity purity
#' @param normal_cn normal total copy number
#' @param total_cn (average) total copy number in cancer population
#' @param mult (average) mutation multplicity
#' @param epi sequencing error
#' @return ccf
#' @export
MapVaf2CcfLinear <- function ( vaf, purity, normal_cn, total_cn, mult, epi = 1e-3 ) {
  w = ( purity * (mult * (1-epi) - total_cn * epi ) ) / ( (1-purity) * normal_cn + purity * total_cn )
  return( (vaf-epi)/w )
}

#' Convert ccf to vaf using ccube linear mapping
#' @param ccf ccf
#' @param purity purity
#' @param normal_cn normal total copy number
#' @param total_cn (average) total copy number in cancer population
#' @param mult (average) mutation multplicity
#' @param epi sequencing error
#' @return vaf
#' @export
MapCcf2VafLinear <- function ( ccf, purity, normal_cn, total_cn, mult, epi = 1e-3 ) {
  w = ( purity * (mult * (1-epi) - total_cn * epi ) ) / ( (1-purity) * normal_cn + purity * total_cn )
  return( w*ccf + epi )
}

############ Vector dot product ############
# handle single row matrix by multiplying each value
# but not sum them up
dot.ext <- function(x,y,mydim) {

  if(missing(mydim)) pracma::dot(x,y) else {

    if(1 %in% pracma::size(x) & mydim == 1) x * y else pracma::dot(x,y)
  }
}

############ Logarithmic Multivariate Gamma function ############
# Compute logarithm multivariate Gamma function.
# Gamma_p(x) = pi^(p(p-1)/4) prod_(j=1)^p Gamma(x+(1-j)/2)
# log Gamma_p(x) = p(p-1)/4 log pi + sum_(j=1)^p log Gamma(x+(1-j)/2)
logmvgamma <- function(x, d) {

  s <- pracma::size(x)

  x <- matrix(as.numeric(x), nrow=1)

  x <- bsxfun.se("+", kronecker(matrix(1,d,1), x), (1 - matrix(1:d))/2)

  y <- d*(d-1)/4*log(pi) + colSums(lgamma(x))

  y <- matrix(as.numeric(y), nrow=s[1], ncol=s[2])

  y
}

logChoose <- function(n, k) {
  return(
    lgamma(n + 1) - lgamma(k+1) - lgamma(n-k+1)
  )
}


#' Parse old cna_data.txt files
#' @param ssm ssms
#' @param cna copy number
#' @return ccube data frame
#' @export
ParseSnvCnaOld <- function (ssm, cna) {
  ssm$major_cn = 1
  ssm$minor_cn = 1
  ssm$cn_frac = 1
  ssm$mu_r <- NULL
  ssm$mu_v <- NULL
  ssm$cn_ref <- NULL
  ssm$cn_tot <- NULL

  cna <- na.omit(cna)
  cna <- dplyr::filter(cna, ssms!="" )

  if (nrow(cna)>0) {
    cnatmp1 <- strsplit(as.character(cna$ssms), ";")
    for (j in seq_len(nrow(cna))) {
      if (length(cnatmp1[[j]])==0) { next }
      cnatmp1[[j]] = paste(cnatmp1[[j]], cna[j,]$frac, sep="," )
    }
    cnatmp1 <- unlist(cnatmp1)
    cnatmp2 <- do.call(
      rbind, strsplit(cnatmp1, ",")
    )

    if (is.null(dim(cnatmp2) )) {
      cnatmp2 = as.data.frame(t(cnatmp2), stringsAsFactors=F)
    } else {
      cnatmp2 = as.data.frame(cnatmp2, stringsAsFactors=F)
    }

    for (j in 2:ncol(cnatmp2)) {
      cnatmp2[,j] = as.numeric(cnatmp2[,j])
    }

    ssm <- dplyr::left_join(ssm, cnatmp2, by=c("id"="V1"))
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
  return(ssm)
}

#' Parse PCAWG-11 cna format
#' @param ssm ssms
#' @param cna copy number profile
#' @return ccube data frame
#' @export
ParseSnvCnaPcawg11Format <- function (ssm, cna) {

  id <- do.call(rbind, strsplit(as.character(ssm$gene), "_", fixed = T))
  ssm$chr = id[,1]
  ssm$pos = as.integer(id[,2])
  ssm$mu_r <- NULL
  ssm$mu_v <- NULL
  ssm$cn_ref <- NULL
  ssm$cn_tot <- NULL
  ssm$cn_frac = NA
  ssm$major_cn = NA
  ssm$minor_cn = NA

  for (jj in seq_len(nrow(cna)) ) {
    cc = cna[jj,]
    idx = which(ssm$chr == cc$chromosome &  (ssm$pos >= cc$start & ssm$pos <= cc$end) )
    if (length(idx) > 0) {
      ssm[idx, ]$major_cn <- cc$major_cn
      ssm[idx, ]$minor_cn <- cc$minor_cn
      ssm[idx, ]$cn_frac <- cc$ccf
    }
  }
  ssm$chr <-NULL
  ssm$pos <-NULL
  ssm <- dplyr::filter(ssm, !is.na(major_cn) & !is.na(minor_cn) & !is.na(cn_frac) & major_cn > 0)
  return(ssm)
  return(ssm)
}

#' Parse new consensus files
#' @param ssm ssms
#' @param cna copy number
#' @return ccube data frame
#' @export
ParseSnvCnaPreConsensus <- function(ssm, cna) {

  id <- do.call(rbind, strsplit(as.character(ssm$gene), "_", fixed = T))
  ssm$chr = as.integer(id[,1])
  ssm$pos = as.integer(id[,2])
  ssm$mu_r <- NULL
  ssm$mu_v <- NULL
  ssm$cn_ref <- NULL
  ssm$cn_tot <- NULL
  ssm$cn_frac = 1
  ssm$major_cn = NA
  ssm$minor_cn = NA
  ssm$star = NA
  for (jj in seq_len(nrow(cna)) ) {
    cc = cna[jj,]
    idx = which(ssm$chr == cc$chromosome &  (ssm$pos >= cc$start & ssm$pos <= cc$end) )
    if (length(idx) > 0) {
      ssm[idx, ]$major_cn <- cc$major_cn
      ssm[idx, ]$minor_cn <- cc$minor_cn
      ssm[idx, ]$star <- cc$star
    }
  }
  ssm$chr <-NULL
  ssm$pos <-NULL
  ssm <- dplyr::filter(ssm, !is.na(major_cn) & !is.na(minor_cn) & !is.na(star) & major_cn > 0)
  return(ssm)
}

#' Parse final consensus files
#' @param ssm ssms
#' @param cna copy number
#' @return ccube data frame
#' @export
ParseSnvCnaConsensus <- function(ssm, cna) {

  id <- do.call(rbind, strsplit(as.character(ssm$gene), "_", fixed = T))
  ssm$chr = id[,1]
  ssm$pos = as.integer(id[,2])
  ssm$mu_r <- NULL
  ssm$mu_v <- NULL
  ssm$cn_ref <- NULL
  ssm$cn_tot <- NULL
  ssm$cn_frac = 1
  ssm$major_cn = NA
  ssm$minor_cn = NA
  ssm$star = NA
  ssm$level = NA
  for (jj in seq_len(nrow(cna)) ) {
    cc = cna[jj,]
    idx = which(ssm$chr == cc$chromosome &  (ssm$pos >= cc$start & ssm$pos <= cc$end) )
    if (length(idx) > 0) {
      ssm[idx, ]$major_cn <- cc$major_cn
      ssm[idx, ]$minor_cn <- cc$minor_cn
      ssm[idx, ]$star <- cc$star
      ssm[idx, ]$level <- cc$level
    }
  }
  ssm$chr <-NULL
  ssm$pos <-NULL
  ssm <- dplyr::filter(ssm, !is.na(major_cn) & !is.na(minor_cn) & major_cn > 0)
  return(ssm)
}

#' Parse Battenberg copy number files
#' @param ssm ssms
#' @param cna copy number
#' @return ccube data frame
#' @export
ParseSnvCnaBattenberg <- function(ssm, cna) {

  id <- do.call(rbind, strsplit(as.character(ssm$gene), "_", fixed = T))
  ssm$chr = id[,1]
  ssm$pos = as.integer(id[,2])
  ssm$major_cn_sub1 = NA
  ssm$minor_cn_sub1 = NA
  ssm$frac_cn_sub1 = NA
  ssm$major_cn_sub2 = -100
  ssm$minor_cn_sub2 = -100
  ssm$frac_cn_sub2 = 0
  ssm$mu_r <- NULL
  ssm$mu_v <- NULL

  for (jj in seq_len(nrow(cna)) ) {
    cc = cna[jj,]
    idx = which(ssm$chr == cc$chr &  (ssm$pos >= cc$startpos & ssm$pos <= cc$endpos) )
    if (length(idx) > 0) {
      ssm[idx, ]$major_cn_sub1 <- cc$nMaj1_A
      ssm[idx, ]$minor_cn_sub1 <- cc$nMin1_A
      ssm[idx, ]$frac_cn_sub1 <- cc$frac1_A
      ssm[idx, ]$frac_cn_sub2 <- 1 - cc$frac1_A

      if (  !is.na(  cc$nMaj2_A ) ) {
        ssm[idx, ]$major_cn_sub2 <- cc$nMaj2_A
        ssm[idx, ]$minor_cn_sub2 <- cc$nMin2_A
      } else {
        ssm[idx, ]$major_cn_sub2 <- -100
        ssm[idx, ]$minor_cn_sub2 <- -100
      }
    }
  }



  ssm$normal_cn = 2
  ssm <- dplyr::rename(ssm, ref_counts=a, total_counts=d, mutation_id = gene)
  ssm <- dplyr::mutate(ssm, var_counts=total_counts-ref_counts)

  HasNonOverLappingSsm <-  sum( is.na(ssm$major_cn_sub1) ) > 0
  if ( HasNonOverLappingSsm ) {

    nonOverLappingSsm <- dplyr::filter(ssm, is.na(major_cn_sub1) )

    for ( ii in 1:nrow(nonOverLappingSsm)  ) {
      ref_range_ssm = which(ssm$chr == nonOverLappingSsm[ii, ]$chr & !is.na(ssm$major_cn_sub1) )
      ref_idx = which.min(  abs( ssm[ref_range_ssm,]$pos - nonOverLappingSsm[ii, ]$pos) )
      idx = ref_range_ssm[ref_idx]
      nonOverLappingSsm[ii, ]$major_cn_sub1 <- ssm[idx, ]$major_cn_sub1
      nonOverLappingSsm[ii, ]$minor_cn_sub1 <- ssm[idx, ]$minor_cn_sub1
      nonOverLappingSsm[ii, ]$frac_cn_sub1 <- ssm[idx, ]$frac_cn_sub1
      nonOverLappingSsm[ii, ]$frac_cn_sub2 <- ssm[idx, ]$frac_cn_sub2
      nonOverLappingSsm[ii, ]$major_cn_sub2 <- ssm[idx, ]$major_cn_sub2
      nonOverLappingSsm[ii, ]$minor_cn_sub2 <- ssm[idx, ]$minor_cn_sub2
    }

    ssm[which(ssm$id %in% nonOverLappingSsm$id), ] = nonOverLappingSsm

  }

  # check gender
  maleCna <- dplyr::filter(cna, chr %in% c("Y", "y") & !is.na(nMaj1_A) )
  isMale <- nrow(maleCna) > 0

  if (isMale) {
    ssm <- dplyr::mutate(dplyr::rowwise(ssm),
                         normal_cn = if (chr %in% c("X", "Y", "x", "y") ) {1} else {2}  )
  }

  ssm$chr <-NULL
  ssm$pos <-NULL
  ssm
}


#' Assign data to centers
#' @param x datas
#' @param centers cluster centers
#' @param s cluster variances
#' @return assignment probabilites and log probabilites
#' @export
Assign <- function(x, centers, s) {
  n <- length(x)
  k <- length(centers)
  logRho <- array(0, dim= c(n ,k))

  for (ii in 1:k) {
    logRho[,ii] = bsxfun.se("-", -(x-centers[ii])^2/(2*s[ii]), log(s[ii]))
  }

  if (n==k) {
    logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)  # 10.49
  } else {
    logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1)) # 10.49
  }

  R <- exp(logR)
  return(list(R=R, logR=logR))
}

#' Ccube color
#' @param n of colors
#' @return colors
#' @export
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


#' Simulate a clonal copy number profile
#' @param cnPoolMaj a pool of possible major copy numbers
#' @param cnPoolMin a pool of possible minor copy numbers
#' @param cnPoolMajFractions prevalence of possible major copy numbers
#' @param cnPoolMinFractions prevalence of possible minor copy numbers
#' @param numVariants number of variants
#' @return cnProfile, a matrix with the minor copy number, major copy number, and total copy number as its 1st, 2nd, and 3rd column
#' @export
GenerateCopyNumberProfile <- function(cnPoolMaj, cnPoolMin,cnPoolMajFractions, cnPoolMinFractions, numVariants){
  tmp1 <- sample(cnPoolMaj, numVariants, cnPoolMajFractions, replace =T)
  tmp2 <- sample(cnPoolMin, numVariants, cnPoolMinFractions, replace =T)
  cnProfileTot <- tmp1 + tmp2
  cnProfile <- cbind(tmp1, tmp2, cnProfileTot )
  cnProfile <- t( apply(cnProfile, 1, sort) )
  return(cnProfile)
}

#' Simulate a 2-subclone subclonal copy number profile.
#' @param cnPoolMaj a pool of possible major copy numbers
#' @param cnPoolMin a pool of possible minor copy numbers
#' @param cnPoolMajFractions prevalence of possible major copy numbers
#' @param cnPoolMinFractions prevalence of possible minor copy numbers
#' @param numVariants number of variants
#' @param subclonal a binary indicator of subclonal status, 1 for subclonal, 0 for clonal.
#' @param ccfCN prevalence of the subclonal copy numbers
#' @export
GenerateSubClonalCNProfile <- function(cnPoolMaj, cnPoolMin,
                                       cnPoolMajFractions, cnPoolMinFractions,
                                       numVariants, subclonal, ccfCN) {
  cnProfile <- matrix(-100, nrow = numVariants, ncol = 8)
  for (ii in 1:numVariants) {
    if (subclonal[ii]) {
      tmp1 = tmp2 = NA
      while (  identical(tmp1, tmp2) ) {
        tmp1 = GenerateCopyNumberProfile(cnPoolMaj, cnPoolMin,cnPoolMajFractions, cnPoolMinFractions, 1)
        tmp2 = GenerateCopyNumberProfile(cnPoolMaj, cnPoolMin,cnPoolMajFractions, cnPoolMinFractions, 1)
      }
      cnProfile[ii, ] = c(tmp1, ccfCN[1], tmp2, ccfCN[2])
    } else {
      cnProfile[ii, 1:4] = c( GenerateCopyNumberProfile(cnPoolMaj, cnPoolMin,
                                                        cnPoolMajFractions, cnPoolMinFractions,
                                                        1), 1)
      cnProfile[ii, 8] = 0
    }
  }

  return(cnProfile)
}

python_false_true_converter <- function(x) {
  if (x == "False") {
    x = FALSE
  }

  if (x == "True") {
    x = TRUE
  }
  x
}



#' Combine a SNV and a SV result list.
#' @param resSNV SNV result list
#' @param resSV SV result list
#' @return a list containing combined results
#' @export
CombineSNVandSVResults <- function(resSNV, resSV) {
  fm1 = resSNV$full.model
  fm2 = resSV$full.model

  fm =
    list(
      ccfMean = c(fm1$ccfMean, fm2$ccfMean),
      ccfCov = c(fm1$ccfCov, fm2$ccfCov),
      dirichletConcentration = c(),
      Epi = c(),
      responsibility = c(),
      logResponsibility = c(),
      dirichletConcentration0 = mean(fm1$dirichletConcentration/sum(fm1$dirichletConcentration) ),
      normalMean = fm1$normalMean,
      invWhishartScale = fm1$invWhishartScale,
      bv = fm1$bv,
      bv_sub1 = fm1$bv_sub1,
      bv_sub2 = fm1$bv_sub2,
      bv1 = fm2$bv1,
      bv2 = fm2$bv2,
      bv1_sub1 = fm2$bv1_sub1,
      bv1_sub2 = fm2$bv1_sub2,
      bv2_sub1 = fm2$bv2_sub1,
      bv2_sub2 = fm2$bv2_sub1
      )

  return( list(full.model = fm) )
}

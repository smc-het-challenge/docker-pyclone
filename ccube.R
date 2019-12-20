############ Sort the model  ############
# Sort model paramters in increasing order of averaged means
# of d variables
SortClusters <- function(model) {

  if ( length(model$ccfMean) > 1 ) {

    if (is.matrix(model$ccfMean)) {
      idx <- order(apply(model$ccfMean, 2, mean))

      model$ccfMean <- model$ccfMean[, idx]

      names(model$ccfMean) <- NULL

      model$ccfCov <- model$ccfCov[,idx]

      names(model$ccfCov) <- NULL

    } else {

      idx <- order(model$ccfMean)

      model$ccfMean <- model$ccfMean[idx]

      names(model$ccfMean) <- NULL

      model$ccfCov <- model$ccfCov[idx]

      names(model$ccfCov) <- NULL
    }


    model$responsibility <- model$responsibility[, idx]

    model$logResponsibility <- model$logResponsibility[, idx]

    if (!is.matrix(model$responsibility) ) {
      model$responsibility <- t ( as.matrix(model$responsibility)  )
      model$logResponsibility <- t ( as.matrix(model$responsibility)  )
    }

    model$dirichletConcentration <- model$dirichletConcentration[idx]

  } else if (! is.matrix(model$responsibility) & length(model$responsibility) == 1 ) {

    model$responsibility <- as.matrix(model$responsibility)
    model$logResponsibility <- as.matrix(model$logResponsibility)

  }

  model
}

#' Run Ccube with model 6: Normal-Binomial. TODO: to be removed in the final release
#' @param mydata mutation data frame
#' @param epi sequencing error
#' @param init scalar number of clusters, or vector all possible cluster numbers
#' @param prior prior parameters
#' @param tol VBEM convergence threshold
#' @param maxiter VBEM maximum iteration
#' @param fit_mult whether or not to estimate multiplicities
#' @param fit_hyper whether or not to estimat hyperparameters
#' @param use methods to get rough estimates of ccf
#' @param verbose show progress
#' @export
#' @return a list containing model parameters
ccube_m6 <- function(mydata, epi=1e-3, init=2, prior=NULL, tol=1e-20, maxiter=1e3, fit_mult = F, fit_hyper = T, use = c("use_base", "use_one"),verbose=FALSE) {

  stopifnot(
    all(c("var_counts","ref_counts","normal_cn",
          "major_cn","minor_cn","purity") %in% names(mydata)))
  mydata <- GetCcf(mydata, use=use)

  dn <- mydata$ref_counts + mydata$var_counts
  bn <- mydata$var_counts
  cn <- mydata$normal_cn
  cr <- mydata$major_cn + mydata$minor_cn
  major_cn <- mydata$major_cn
  purity <- unique(mydata$purity)
  bv <- mydata$rough_mult
  rawCcf <- mydata$rough_ccf
  rawCcf <- as.matrix(rawCcf)
	n <- nrow(rawCcf)
	d <- ncol(rawCcf)

	X <- t(rawCcf) # Work with D by N for convenience

	message(sprintf("Running VB-Normal-Binomial on a %d-by-%d data with %d clusters ...\n", n, d, init))

	if(is.null(prior)) {

		if( length(init) == 1 ) {
			# General prior with equal dirichletConcentration
		  k <-  init
		  prior <- list(
		    dirichletConcentration = 1e-2,
		    normalMean = 1,
		    invWhishartScale = var(rawCcf)*(d+1)# M = inv(W)
			)
		}
	} else {

		stopifnot(
			all(names(prior) %in% c("dirichletConcentration","normalMean","invWhishartScale")) &
			all(sapply(prior, is.numeric)) & nrow(prior$normalMean) == d &
			ncol(prior$normalMean) == 1 &
			nrow(prior$normalMean) == d & ncol(prior$invWhishartScale) == d)
	}

	# variational lower bound (objective function)
	L <- rep(-Inf, maxiter)
	converged <- FALSE
	degenerated <- FALSE
	vbiter <- 1

	model <- list()

	initParams <- initialization(X, init, prior) # initialize responsibility and hidden scale
	model$responsibility <- initParams$R
	model$ccfMean <- initParams$ccfMean
	model$ccfCov <- initParams$ccfCov
	model$bv <- bv
	model$dirichletConcentration0 <- prior$dirichletConcentration
	model$normalMean <- prior$normalMean
	model$invWhishartScale <- prior$invWhishartScale

	while(!converged & vbiter < maxiter & !degenerated) {
	  vbiter <- vbiter + 1
		model <- VariationalMaximimizationStep(bn, dn, cn, cr, major_cn, epi, purity, model,
		                                       fit_mult = fit_mult, fit_hyper = fit_hyper)
		model <- VarationalExpectationStep(bn, dn, cn, cr, epi, purity, model)
		L[vbiter] <- VariationalLowerBound(bn, dn, cn, cr, epi, purity, model)/n
		converged <- abs(L[vbiter] - L[vbiter-1]) < (tol) * abs(L[vbiter])
		degenerated <- (L[vbiter] - L[vbiter-1]) < 0
		#degenerated = F
		if(verbose) cat(sprintf("\rVB-EM-%d: L = %.8f \r", vbiter, L[vbiter]))
	}

	L <- L[2:vbiter]

	model <- SortClusters(model)

	if (init > 1) {
	  label <- apply(model$responsibility, 1, which.max)
	  nk <- colSums(model$responsibility)
	  Epi <- (model$dirichletConcentration + nk) / (k*model$dirichletConcentration0 + n)
	  model$Epi <- Epi/sum(Epi)
	} else {
	  label <- rep(1, n)
	  model$Epi <- 1
	}


	if(converged) cat(sprintf("\nConverged in %d steps.\n", vbiter-1)) else
	  cat(sprintf("Not converged in %d steps.\n", maxiter))
	if(degenerated)
	  cat(sprintf("Degenerated at %d steps.\n", vbiter-1))

	list(label=label, R=model$responsibility, mu=model$ccfMean, full.model=model, L=L)
}

#' Run Ccube Core function
#' @param mydata mutation data frame
#' @param epi sequencing error
#' @param init scalar number of clusters, or vector all possible cluster numbers
#' @param prior prior parameters
#' @param tol VBEM convergence threshold
#' @param maxiter VBEM maximum iteration
#' @param fit_mult whether or not to estimate multiplicities
#' @param fit_hyper whether or not to estimat hyperparameters
#' @param use methods to get rough estimates of ccf
#' @param verbose show progress
#' @export
#' @return a list containing model parameters
CcubeCore <- function(mydata, epi=1e-3, init=2, prior, tol=1e-20, maxiter=1e3, fit_mult = F, fit_hyper = T, use = c("use_base", "use_one"),verbose=FALSE) {

  stopifnot(
    all(c("var_counts","ref_counts","normal_cn",
          "purity", "subclonal_cn") %in% names(mydata)))

  stopifnot(
    all(c("major_cn_sub1","major_cn_sub2","minor_cn_sub1", "minor_cn_sub2",
          "frac_cn_sub1", "frac_cn_sub2")
        %in% names(mydata)))


  mydata <- GetCcf(mydata, use=use)

  dn <- mydata$ref_counts + mydata$var_counts
  bn <- mydata$var_counts
  cn <- mydata$normal_cn
  cr <- mydata$total_cn

  max_mult_cn_sub1 <- mydata$major_cn_sub1
  max_mult_cn_sub2 <- mydata$major_cn_sub2
  frac_cn_sub1 <- mydata$frac_cn_sub1
  frac_cn_sub2 <- mydata$frac_cn_sub2
  bv <- mydata$rough_mult
  bv_sub1 <- rep(-100, length(max_mult_cn_sub1))
  bv_sub2 <- rep(-100, length(max_mult_cn_sub2))
  subclonal_cn <- mydata$subclonal_cn

  purity <- unique(mydata$purity)
  rawCcf <- mydata$rough_ccf
  rawCcf <- as.matrix(rawCcf)
  n <- nrow(rawCcf)
  d <- ncol(rawCcf)

  X <- t(rawCcf) # Work with D by N for convenience

  message(sprintf("Running VB-Normal-Binomial on a %d-by-%d data with %d clusters ...\n", n, d, init))

  if(missing(prior)) {

    if( length(init) == 1 ) {
      # General prior with equal dirichletConcentration
      k <-  init
      prior <- list(
        dirichletConcentration = 1e-2,
        normalMean = 1,
        invWhishartScale = var(rawCcf)*(d+1)# M = inv(W)
      )
    }
  } else {

    stopifnot(
      all(names(prior) %in% c("dirichletConcentration","normalMean","invWhishartScale")) &
        all(sapply(prior, is.numeric)) & nrow(prior$normalMean) == d &
        ncol(prior$normalMean) == 1 &
        nrow(prior$normalMean) == d & ncol(prior$invWhishartScale) == d)
  }

  # variational lower bound (objective function)
  L <- rep(-Inf, maxiter)
  converged <- FALSE
  degenerated <- FALSE
  vbiter <- 1

  model <- list()

  initParams <- initialization(X, init, prior) # initialize responsibility and hidden scale
  model$responsibility <- initParams$R
  model$ccfMean <- initParams$ccfMean
  model$ccfCov <- initParams$ccfCov
  model$bv <- bv
  model$bv_sub1 <- bv_sub1
  model$bv_sub2 <- bv_sub2
  model$dirichletConcentration0 <- prior$dirichletConcentration
  model$normalMean <- prior$normalMean
  model$invWhishartScale <- prior$invWhishartScale

  while(!converged & vbiter < maxiter & !degenerated) {
    vbiter <- vbiter + 1
    model <- VariationalMaximimizationStep(bn, dn, cn, cr, max_mult_cn_sub1, max_mult_cn_sub2,
                                           frac_cn_sub1, frac_cn_sub2, subclonal_cn, epi, purity, model,
                                           fit_mult = fit_mult, fit_hyper = fit_hyper)
    model <- VarationalExpectationStep(bn, dn, cn, cr, epi, purity, model)
    L[vbiter] <- VariationalLowerBound(bn, dn, cn, cr, epi, purity, model)/n
    converged <- abs(L[vbiter] - L[vbiter-1]) < (tol) * abs(L[vbiter])
    degenerated <- (L[vbiter] - L[vbiter-1]) < 0
    #degenerated = F
    if(verbose) cat(sprintf("\rVB-EM-%d: L = %.8f \r", vbiter, L[vbiter]))
  }

  L <- L[2:vbiter]

  model <- SortClusters(model)

  if (init > 1) {
    label <- apply(model$responsibility, 1, which.max)
    nk <- colSums(model$responsibility)
    Epi <- (model$dirichletConcentration + nk) / (k*model$dirichletConcentration0 + n)
    model$Epi <- Epi/sum(Epi)
  } else {
    label <- rep(1, n)
    model$Epi <- 1
  }


  if(converged) cat(sprintf("\nConverged in %d steps.\n", vbiter-1)) else
    cat(sprintf("Not converged in %d steps.\n", maxiter))
  if(degenerated)
    cat(sprintf("Degenerated at %d steps.\n", vbiter-1))

  list(label=label, full.model=model, L=L)
}



#' Get a rough estimate of ccf and multiplicities
#' @param mydata mutation data frame
#' @param use methods for get rough estimates of ccf
#' @return mydata mutation data frame
#' @export
GetCcf <- function(mydata, use = c("use_base", "use_one")) {

  if ("rough_ccf" %in% names(mydata) &
      "rough_mult" %in% names(mydata)){
    return(mydata)
  }

  GetMult <- function(x, y, z) {
    index <- which(x==z)
    mean(y[index])
  }

  if (!"total_counts" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, total_counts = ref_counts + var_counts)
  }



  if (!"major_cn" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, major_cn = frac_cn_sub1 * major_cn_sub1 +
                              frac_cn_sub2 * major_cn_sub2)
  }

  if (!"minor_cn" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, minor_cn = frac_cn_sub1 * minor_cn_sub1 +
                              frac_cn_sub2 * minor_cn_sub2 )
  }

  if (!"total_cn" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, total_cn = major_cn + minor_cn )
  }

  if (!"vaf" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, vaf = var_counts/total_counts )
  }


  if (use=="use_base") {
    mydata <- dplyr::mutate(dplyr::rowwise(mydata),
                            ccf1 = MapVaf2CcfPyClone(vaf,
                                                     purity,
                                                     normal_cn,
                                                     total_cn,
                                                     total_cn,
                                                     major_cn, lower = 0, upper = 2),
                            ccf2 = MapVaf2CcfPyClone(vaf,
                                                     purity,
                                                     normal_cn,
                                                     total_cn,
                                                     total_cn,
                                                     minor_cn, lower = 0, upper = 2),
                            ccf3 = MapVaf2CcfPyClone(vaf,
                                                     purity,
                                                     normal_cn,
                                                     total_cn,
                                                     total_cn,
                                                     1, lower = 0, upper = 2))

    mydata <- dplyr::mutate(dplyr::rowwise(mydata), rough_ccf = mean(c(ccf1, ccf2, ccf3), na.rm = T))
    dd <- dplyr::filter(mydata, ccf1 != ccf2 | ccf1 != ccf3 | ccf2 != ccf3)

    if (nrow(dd) > 0) {
      dd <- dplyr::mutate(dplyr::rowwise(dd), rough_ccf = min(c(ccf1, ccf2, ccf3), na.rm = T ))
      mydata[mydata$mutation_id %in% dd$mutation_id,]$rough_ccf = dd$rough_ccf
    }

    mydata <- dplyr::mutate(dplyr::rowwise(mydata),
                            rough_mult = GetMult(c(ccf1, ccf2, ccf3),
                                           c(major_cn, minor_cn, 1), rough_ccf))

    dd1 <- dplyr::filter(mydata, is.na(rough_ccf))
    if (nrow(dd1) > 0) {
      dd1 <- dplyr::mutate(dplyr::rowwise(dd1),
                           rough_ccf = MapVaf2CcfPyClone(vaf,
                                                   purity,
                                                   normal_cn,
                                                   total_cn,
                                                   total_cn,
                                                   1, constraint = F))


      mydata[mydata$mutation_id %in% dd1$mutation_id,]$rough_ccf = dd1$rough_ccf
      mydata[mydata$mutation_id %in% dd1$mutation_id,]$rough_mult = 1
    }
  }

  if (use == "use_one") {
    mydata$rough_mult = 1
    mydata <- dplyr::mutate(dplyr::rowwise(mydata),
                            rough_ccf = MapVaf2CcfPyClone(vaf,
                                                 purity,
                                                 normal_cn,
                                                 total_cn,
                                                 total_cn,
                                                 1, constraint = F))
  }

  mydata$ccf1 <- NULL
  mydata$ccf2 <- NULL
  mydata$ccf3 <- NULL
  return(mydata)
}


############ helper function to initialize responsibility and other parameters ############
initialization <- function(X, init, prior) {

  d <- nrow(X)

  n <- ncol(X)

  stopifnot(length(init) %in% c(1, n) ||
              (nrow(init) == d  & ncol(init) == k))

  k <- init
  numOfDistinctDataPoints <- dplyr::n_distinct(X)

  if ( k > numOfDistinctDataPoints | k == numOfDistinctDataPoints  ) {

    normalMean <- t(rep(1, k))
    R <- as.matrix(Matrix::sparseMatrix(1:n, sample(1:k, n, replace = T), x=1, dims = c(n, k)))

  } else {

    res <- kmeans(t(X), k)
    label <- res$cluster
    normalMean <- t(res$centers)
    if (length(which(normalMean > 1) ) >0 ){
      normalMean[which(normalMean>1)] = 1
    }
    R <- as.matrix(Matrix::sparseMatrix(1:n, label, x=1, dims = c(n, k)))

  }

  invWhishartScale <- prior$invWhishartScale
  ccfMean = unname(normalMean)
  ccfCov = array(invWhishartScale, dim = c(d,k)) # TODO temporal fix

  return(list(R = R,
              ccfMean = ccfMean, ccfCov = ccfCov))
}


############ Variational-Maximimization ############
VariationalMaximimizationStep <- function(bn, dn, cn, cr, max_mult_cn_sub1, max_mult_cn_sub2,
                                          frac_cn_sub1, frac_cn_sub2, subclonal_cn,
                                          epi, purity, model, fit_mult = T, fit_hyper = T, fit_params_ccf = T,
                                          fit_params_pi = T) {

  bv <- model$bv
  bv_sub1 <- model$bv_sub1
  bv_sub2 <- model$bv_sub2

  dirichletConcentration0 <- model$dirichletConcentration0
  responsibility <- model$responsibility
  normalMean <- model$normalMean
  invWhishartScale <- model$invWhishartScale



  # Gaussian approximation for each clone
  ccfMean = model$ccfMean
  ccfCov = model$ccfCov

  if (fit_params_ccf){
    k <- length(ccfMean)
    Bn = (1-purity)*cn + purity*cr
    Cn = purity*(bv*(1-epi) - cr*epi)

    ccfMeanOld <- ccfMean

    for (i in 1:k){
      term1 = 1/invWhishartScale*normalMean
      term2 = 1/invWhishartScale
      upper <- 1
      lower <- 0
      tmp <- NULL
      jj <- 0
      while (!is.numeric(tmp)) {
        if (jj >=  1000) {
          tmp <- ccfMeanOld[i]
          break
        }
        jj <- jj + 1
        if (jj > 1) {upper <- upper + 0.1}
        tmp <- try(suppressWarnings( uniroot(
          function(x) {
            An = (1-purity)*cn*epi + purity*(1-x)*cr*epi + purity*x*bv*(1-epi)
            term4 = (bn/An - (dn-bn)/(Bn-An)) * Cn
            term5 = term2*x
            return(sum(responsibility[,i]*term4) - term5 + term1)
          },
          c(lower, upper), extendInt = "yes")$root), T)
      }
      ccfMean[i] <- tmp

      An = (1-purity)*cn*epi + purity*(1-ccfMean[i])*cr*epi + purity*ccfMean[i]*bv*(1-epi)
      ccfCov[i] = solve(term2 + sum(responsibility[,i]*(bn/(An^2) + (dn-bn)/((Bn-An)^2) * Cn^2)))
    }


    if (length(which(ccfMean > 1) ) >0 ){
      ccfMean[which(ccfMean>1)] = 1
    }
    if (length(which(ccfMean < 1e-20)) > 0){
      ccfMean[which(ccfMean<1e-20)] = 1e-20
    }

    model$ccfMean <- ccfMean
    model$ccfCov <- ccfCov

  }


  if (fit_params_pi) {
    dirichletConcentration <- model$dirichletConcentration0 + colSums(responsibility) # 10.58
    model$dirichletConcentration <- dirichletConcentration
  }




  numberOfDataPoints <- nrow(responsibility)
  if (fit_mult) {
    for (ii in 1:numberOfDataPoints) {
      # TODO: General implementation of more than two subclonal copy number populations

      if (subclonal_cn[ii]) {

        sub_cn_mults = pracma::meshgrid(0:max_mult_cn_sub1[ii], 0:max_mult_cn_sub2[ii])
        bvPool <- frac_cn_sub1[ii] * sub_cn_mults$X + frac_cn_sub2[ii] * sub_cn_mults$Y
        bvPool[ bvPool < 0.6 ] = NA
        bvPoolMat <- t(my_repmat(as.vector(bvPool), length(ccfMean)))
        ccfMeanMat <- my_repmat(ccfMean, length(bvPool))
        ccfCovMat <- my_repmat(ccfCov, length(bvPool))
        respMat <- my_repmat(responsibility[ii, ], length(bvPool))
        w <- purity * (bvPoolMat*(1-epi) -cr[ii]*epi) / ((1-purity)*cn[ii] + purity * cr[ii])
        ww <- w^2
        ef <- w * ccfMeanMat +epi
        term1 <- bn[ii] * (log (ef) - ww*ccfCovMat/(2 * ef^2 ) )
        term2 <- (dn[ii] - bn[ii]) * (log (1 - ef) - ww*ccfCovMat/(2 * (1 - ef)^2)  )
        term3 <- logChoose(dn[ii], bn[ii])
        qq <- rowSums ( respMat *  (term1 + term2 + term3))
        qq[1] = NA # remove both multiplicities being zero
        maxQq <- which.max(qq)


        if ( length(maxQq) == 0 ) {

          if (frac_cn_sub1[ii] >= frac_cn_sub2[ii]) {

            if (max_mult_cn_sub1[ii] > 0) {
              bv[ii] <- frac_cn_sub1[ii]
              bv_sub1[ii] = 1
              bv_sub2[ii] = 0
            } else {
              bv[ii] <- frac_cn_sub2[ii]
              bv_sub1[ii] = 0
              bv_sub2[ii] = 1
            }

          } else {
            if (max_mult_cn_sub2[ii] > 0) {
              bv[ii] <- frac_cn_sub2[ii]
              bv_sub1[ii] = 0
              bv_sub2[ii] = 1
            } else {
              bv[ii] <- frac_cn_sub1[ii]
              bv_sub1[ii] = 1
              bv_sub2[ii] = 0
            }

          }
        } else {
          bv[ii] <- bvPool[maxQq]
          bv_sub1[ii] <- sub_cn_mults$X[maxQq]
          bv_sub2[ii] <- sub_cn_mults$Y[maxQq]
        }

        # Todo: remove
        # if ( length(maxQq) == 0 ) {
        #   if (frac_cn_sub1[ii] >= frac_cn_sub2[ii]) {
        #     bv[ii] <- frac_cn_sub1[ii]
        #     bv_sub1[ii] = 1
        #     bv_sub2[ii] = 0
        #   } else {
        #     bv[ii] <- frac_cn_sub2[ii]
        #     bv_sub1[ii] = 0
        #     bv_sub2[ii] = 1
        #   }
        # } else {
        #   bv[ii] <- bvPool[maxQq]
        #   bv_sub1[ii] <- sub_cn_mults$X[maxQq]
        #   bv_sub2[ii] <- sub_cn_mults$Y[maxQq]
        # }


      } else {

        if ( cr[ii] < 1 ) {

          bv[ii] <- max_mult_cn_sub1[ii]

        } else {

          bvPool <- 1:max_mult_cn_sub1[ii]
          qq <- rep(NA, length(bvPool))
          for (jj in seq_along(bvPool) ) {
            aa <- purity * (bvPool[jj] *(1-epi) -cr[ii]*epi) / ((1-purity)*cn[ii] + purity * cr[ii])
            aa2 <- aa^2
            bb <- epi
            term1 <- sum(responsibility[ii, ] * bn[ii] * (log (aa * ccfMean +bb) - aa2*ccfCov/(2 * (aa * ccfMean +bb)^2 ) ))
            term2 <- sum(responsibility[ii, ] * (dn[ii] - bn[ii]) * (log (1 - aa * ccfMean - bb) - aa2*ccfCov/(2 * (1 - aa * ccfMean -bb)^2)  ))
            term3 <- sum( responsibility[ii, ]*  logChoose(dn[ii], bn[ii]) )
            qq[jj] <- term1 + term2 + term3
          }
          maxQq <- which.max(qq)
          if ( length(maxQq) == 0 ) {
            bv[ii] <- 1
          } else {
            bv[ii] <- bvPool[maxQq]
          }

        }

      }

    }
    model$bv <- bv
    model$bv_sub1 <- bv_sub1
    model$bv_sub2 <- bv_sub2
  }

  # estimate hyper-parameters
  if (fit_hyper) {
    # Elogpi <- sum(digamma(dirichletConcentration) - digamma(sum(dirichletConcentration)))
    # tmpOld <- model$dirichletConcentration0
    # jj <- 0
    # upper <- 1e-2
    # lower <- 1e-99
    # while (!is.numeric(tmp)) {
    #   if (jj > 1000) {
    #     break
    #   }
    #   jj <- jj + 1
    #   if (jj > 1) {upper <- upper + 1e-3}
    #   tmp <- try(suppressWarnings( uniroot(
    #     function(x) {
    #       term1 <- digamma(k*x) - digamma(x)
    #       return(term1+Elogpi/k)
    #     },
    #     c(lower, upper), extendInt = "yes")$root), T)
    # }
    # if (jj > 1000) {
    #   model$dirichletConcentration0 <- tmpOld
    #   cat(jj)
    # } else {
    #   model$dirichletConcentration0 <- tmp
    # }
    model$normalMean <- mean(ccfMean)
    model$invWhishartScale <-mean((ccfMean - model$normalMean)^2 + ccfCov)
  }

  model
}


############ Variational-Expectation ############
VarationalExpectationStep <- function(bn, dn, cn, cr, epi, purity, model, no.weights = FALSE) {

  bv <- model$bv
  dirichletConcentration <- model$dirichletConcentration	# Dirichlet
  ccfMean <- model$ccfMean
  ccfCov <- model$ccfCov
  Epi <- model$Epi

  n <- length(bn)
  if (is.matrix(ccfMean)) {
    d <- nrow(ccfMean)
    k <- ncol(ccfMean)
  } else {
    d <- 1
    k <- length(ccfMean)
  }

  Epbk <- array(0, dim=c(n,k))

  aa <- purity * (bv *(1-epi) -cr*epi) / ((1-purity)*cn + purity * cr)
  aa2 <- aa^2
  bb <- epi

  # todo: only works for d = 1
  for(i in 1:k) {
    Epbk[,i] <- bn * (log (aa * ccfMean[i] +bb) - aa2*ccfCov[i]/(2 * (aa * ccfMean[i] +bb)^2 ) ) +
      (dn - bn) * (log (1 - aa * ccfMean[i] - bb) - aa2*ccfCov[i]/(2 * (1 - aa * ccfMean[i] -bb)^2)  )
  }


  if (no.weights) {
    logRho <-  Epbk  #eq (19)  10.46
  } else {
    Elogpi <- digamma(dirichletConcentration) - digamma(sum(dirichletConcentration))
    logRho <- bsxfun.se("+", Epbk, Elogpi)
  }

  if (n==k) {
    logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
  } else {
    logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
  }

  R <- exp(logR)

  model$logResponsibility <- logR
  model$responsibility<- R

  model
}


############ Variational-(lower)-Bound Evaluation ############
VariationalLowerBound <- function(bn, dn, cn, cr, epi, purity, model) {

  bv <- model$bv
  dirichletConcentration0 <- model$dirichletConcentration0
  m0 <- model$normalMean
  M0 <- model$invWhishartScale

  dirichletConcentration <- model$dirichletConcentration	# Dirichlet
  m <- model$ccfMean			# Gasusian
  M <- model$ccfCov
  R <- model$responsibility
  logR <- model$logResponsibility
  ccfMean <- model$ccfMean
  ccfCov <- model$ccfCov
  nk <- colSums(R)									# 10.51

  n <- length(bn)

  if (is.matrix(m)) {
    d <- nrow(m)
    k <- ncol(m)
  } else {
    d <- 1
    k <- length(m)
  }




  Elogpi <- digamma(dirichletConcentration) - digamma(sum(dirichletConcentration))		# 10.66

  Epz = pracma::dot(nk, Elogpi)								# 10.72
  Eqz = pracma::dot(as.numeric(R), as.numeric(logR))			# 10.75
  logCdirichletConcentration0 = lgamma(k*dirichletConcentration0) - k*lgamma(dirichletConcentration0)
  Eppi <- logCdirichletConcentration0+(dirichletConcentration0-1)*sum(Elogpi)
  logCdirichletConcentration <- lgamma(sum(dirichletConcentration))-sum(lgamma(dirichletConcentration))
  Eqpi = pracma::dot(dirichletConcentration-1, Elogpi) + logCdirichletConcentration				# 10.76

  # part of 10.70
  L <- Epz - Eqz + Eppi - Eqpi

  U0 <- chol(M0)
  trM0W <- array(0, dim = c(1, k))
  mm0Wmm0 <- array(0, dim = c(1, k))

  Epbk <- array(0, dim=c(1,k))

  aa <- purity * (bv *(1-epi) -cr*epi) / ((1-purity)*cn + purity * cr)
  aa2 <- aa^2
  bb <- epi

  for(i in 1:k) {

    Epbk[,i] <- sum(R[,i]*(bn * (log (aa * ccfMean[i] +bb) - aa2*ccfCov[,i]/(2 * (aa * ccfMean[i] +bb)^2 ) ) +
      (dn - bn) * (log (1 - aa * ccfMean[i] - bb) - aa2*ccfCov[,i]/(2 * (1 - aa * ccfMean[i] -bb)^2)))) +
      sum( R[,i]*logChoose(dn, bn) )
    q <- solve(t(U0), m[,i,drop=F]-m0)
    mm0Wmm0[i] <- pracma::dot(q, q)
    U <- chol(M[,i])
    Q <- solve(U0, U)
    trM0W[i] <- pracma::dot(as.numeric(Q), as.numeric(Q))
  }

  # first half of 10.74
  Epmu <- sum(-log(2*pi)-log(M)-mm0Wmm0-trM0W)/2

  Eqmu <- -0.5*sum(1+log(2*pi)+log(M))	# 10.77 (1/2)

  L <- L + Epmu - Eqmu + sum(Epbk)	# 10.70

  L
}

#' Estimate purity
#' @param mydata ccube data
#' @param wgd whole genome duplication status
#' @param K total number of cluster in student't mixture
#' @param th threshold on weights to be eligible for clonal cluster
#' @return purity
#' @export
GetPurity <- function(mydata, wgd=F, K = 6, th = 1.5e-2) {

  vtox<-function(v,nA,nB,tA,tB)
  {
    (nB - nB*v - nA*v) / (tA*v + tB*v + nB - tB -nA*v - nB*v)
  }

  mydata <- CheckAndPrepareCcubeInupts(mydata, estimatePurity = F)

  tmpdata <- dplyr::filter(mydata, frac_cn_sub1 == 1) # Use only clonal CN regions

  tmpdata <- dplyr::filter(mydata, major_cn == minor_cn & major_cn != 0)

  if (!wgd) {
    tmpdata <- dplyr::filter(mydata, major_cn == minor_cn & major_cn == 1)
  }

  if (nrow(tmpdata) == 0) {
    return(NA)
  }
  tmpdata <- dplyr::mutate(dplyr::rowwise(tmpdata), ploidy = major_cn + minor_cn)
  tmpdata <- dplyr::mutate(tmpdata, vaf = var_counts/(var_counts+ref_counts))
  tmpdata <- dplyr::mutate(dplyr::rowwise(tmpdata), cp = vtox(vaf, 2, 0, ploidy/2, ploidy/2))
  tmpdata <- dplyr::filter(tmpdata, !is.infinite(cp) & !is.na(cp))

  if (K >= nrow(tmpdata)) {
    K = nrow(tmpdata) - 1
    if (K == 0) {K = 1}
  }
  if (nrow(tmpdata) == 1) {
    purity <- tmpdata$cp
  } else {
    res <- vbsmm(tmpdata$cp, init = K, tol = 1e-5,  verbose = F)
    uniqLabels <- sort(unique(res$label))
    ww <- uniqLabels[which( (table(res$label)/length(res$label) ) > th)]
    pool <- res$mu[ww]
    maxCp <- max(pool[pool<=1])
    purity  <- if (maxCp > 1) 1 else maxCp
  }

  return(purity)
}

#' Make Ccube results plot
#' @param ssm data
#' @param res Ccube result list
#' @param myColors colors
#' @param printPlot output flag
#' @param fn output file name
#' @return NULL
#' @export
MakeCcubeStdPlot <- function(ssm, res, myColors=gg_color_hue(10), printPlot = F, fn = NULL) {

  if (printPlot) {
    pdf(fn, width=8, height=8)
  }

  par(mfrow=c(2,2))
  plot(ssm$ccube_ccf, ssm$vaf, col = myColors[res$label],
       xlab = "cancer cell fraction", ylab = "variant allele frequecy",
       main = "ccf vs vaf (colored by cluster membership)")
  cellularity <- unique(ssm$purity)
  ssm$total_cn =ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
    ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2)
  uniqueTotCn = unique(ssm$total_cn)
  xx = seq(0,2, length.out = 100)
  for (cn in uniqueTotCn) {
    for (i in 1:cn) {
      lines(MapVaf2CcfPyClone(xx, cellularity, 2, cn, cn, i, constraint = F), xx, lty = 6, col = 80)
    }
  }


  if (is.matrix(res$full.model$ccfMean)) {
    Emu <- res$full.model$ccfMean[,]
  } else {
    Emu <- res$full.model$ccfMean
  }

  if (is.matrix(res$full.model$ccfCov)) {
    Esigma <- res$full.model$ccfCov[,]
  } else {
    Esigma <- res$full.model$ccfCov
  }

  Epi <- res$full.model$Epi

  params <- data.frame(Emu, Esigma, Epi)
  xx <- seq(0,2,  length.out = 1000)
  ll <- 0

  for (j in seq_len(nrow(params))) {
    ll <- ll + params[j,]$Epi * dnorm(xx, mean = params[j,]$Emu, sd = sqrt(params[j,]$Esigma))
  }

  hist(ssm$ccube_ccf, density=20, breaks=20, prob=TRUE,
       main = "ccf histogram +
       cluster uncertainties",
       xlab = "cancer cell fraction")
  lines(xx,ll, lwd=2, col = "darkred")

  numSnv <- table(res$label)
  uniqLabels = unique(res$label)
  names(numSnv) <- as.character(format(round(Emu[sort(uniqLabels)], 2), nsmall = 2))
  barplot(numSnv, las = 2, col = myColors[sort(uniqLabels)],
          xlab = "cluster mean", ylab="number of variants",
          main = "cluster prevalence")
  if (printPlot) {
    dev.off()
  }

}

#' Check if the result has a clonal cluster
#' @param res Ccube result list
#' @return TRUE/FALSE
#' @export
HasClonalCluster <- function (res) {
  assignments <- as.data.frame(table(res$label), stringsAsFactors = F)
  assignments <- mutate(assignments, Var1 = as.integer(Var1))
  mu <- res$full.model$ccfMean[assignments$Var1]
  ci_95 <- 2*sqrt(res$full.model$ccfCov)
  ci <- ci_95[assignments$Var1]
  num_clonal <- sum(assignments$Freq[which(mu-ci >= 0.9 | mu+ci >=0.9)])
  return(num_clonal > 0)
}

#' Check if the result has empty cluster
#' @param res Ccube result list
#' @return TRUE/FALSE
#' @export
HasEmptyCluster <- function (res) {
  return(length(unique(res$label)) < length(res$full.model$ccfMean))
}

#' Remove empty clusters
#' @param res Ccube result list
#' @param ssm Ccub input data
#' @param useEstep flag to use Estep for reassignment
#' @param epi sequencing error
#' @return Ccube result list
#' @export
CullEmptyClusters <- function(res, ssm, useEstep = T, epi = 1e-3) {

  idx <- which(! seq_along(res$full.model$ccfMean) %in% unique(res$label) )

  if (useEstep) {
    return(RemoveClusterAndReassignVariantsWithEstep(res = res, removeIdx = idx, ssm = ssm, epi = epi))
  } else {
    return(RemoveClusterAndReassignVariants(res = res, removeIdx = idx, ssm = ssm))
  }
}

#' Remove small clusters
#' @param res Ccube result list
#' @param ssm Ccub input data
#' @param th threshold for small clusters
#' @param epi sequencing error
#' @param useEstep use Estep for reassignment
#' @return Ccube result list
#' @export
CullSmallClusters <- function(res, ssm, th = 1e-2, epi = 1e-3, useEstep = T) {

  tt <- table(res$label)
  idx <- which( tt/sum(tt) < th )

  if (useEstep) {
    return(RemoveClusterAndReassignVariantsWithEstep(res = res, removeIdx = idx, ssm = ssm, epi = epi))
  } else {
    return(RemoveClusterAndReassignVariants(res = res, removeIdx = idx, ssm = ssm))
  }
}

#' Remove a (or more) cluster and reassign its data if the cluster is nonempty
#' @param res Ccube result list
#' @param removeIdx clusters to remove
#' @param ssm data
#' @param label assigned labels if res doesn't have label variable
#' @return Ccube result list
#' @export
RemoveClusterAndReassignVariants <- function(res, removeIdx, ssm = NULL, label = NULL) {

  if (length(removeIdx) == 0) {
    return(res)
  }

  if (!is.null(res$label)) {
    uniqLabels <- sort(unique(res$label))
  } else {
    uniqLabels <- sort(unique(label))
  }

  remainedLabels <- uniqLabels[which(!uniqLabels %in% removeIdx)]
  uniqLabels[ which(uniqLabels %in% removeIdx) ] <- NA
  newLabels <- match(res$label, uniqLabels)
  reassignIdx <- which(is.na(newLabels))
  reassignSsm <- ssm[reassignIdx, ]

  if (! is.null(res$full.model) ) {
    res$full.model$ccfMean <- res$full.model$ccfMean[-removeIdx]
    res$full.model$ccfCov <- res$full.model$ccfCov[-removeIdx]

    logRho <- res$full.model$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$full.model$logResponsibility <-
      if (length(res$label) ==length(res$full.model$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }

    if(!is.null(ssm) & length(reassignIdx) >0) {
      reassignList <- Assign(reassignSsm$ccube_ccf, res$full.model$ccfMean, res$full.model$ccfCov)
      res$full.model$logResponsibility[reassignIdx, ] <- reassignList$logR
    }

    res$full.model$responsibility <- exp(res$full.model$logResponsibility)

    res$label <- apply(res$full.model$responsibility, 1, which.max)

    res$full.model$dirichletConcentration <- res$full.model$dirichletConcentration0 + colSums(res$full.model$responsibility)

    if (!is.null(res$R)) {
      res$R=res$full.model$responsibility
    }

    if (!is.null(res$mu)) {
      res$mu=res$full.model$ccfMean
    }
    if (! is.null(res$full.model$Epi)) {
      Epi <- (res$full.model$dirichletConcentration + colSums(res$full.model$responsibility)) /
        (length(res$full.model$ccfMean) * res$full.model$dirichletConcentration0 + length(res$label))
      res$full.model$Epi <- Epi/sum(Epi)
    }

  } else {

    res$ccfMean <- res$ccfMean[-removeIdx]
    res$ccfCov <- res$ccfCov[-removeIdx]

    logRho <- res$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$logResponsibility <-
      if (length(label) ==length(res$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }

    if(!is.null(ssm) & length(reassignIdx) >0) {
      reassignList <- Assign(reassignSsm$ccube_ccf, res$ccfMean, res$ccfCov)
      res$logResponsibility[reassignIdx, ] <- reassignList$logR
    }

    res$responsibility <- exp(res$logResponsibility)

    res$dirichletConcentration <- res$dirichletConcentration0 + colSums(res$responsibility)
  }

  res
}

#' Remove a (or more) cluster and reassign its data if the cluster is nonempty
#' @param res Ccube result list
#' @param removeIdx clusters to remove
#' @param ssm data
#' @param label assigned labels if res doesn't have label variable
#' @param epi sequencing error
#' @return Ccube result list
#' @export
RemoveClusterAndReassignVariantsWithEstep <- function(res, removeIdx, ssm = NULL, label = NULL, epi = 1e-3) {

  if (length(removeIdx) == 0) {
    return(res)
  }

  if (!is.null(res$label)) {
    uniqLabels <- sort(unique(res$label))
  } else {
    uniqLabels <- sort(unique(label))
  }

  remainedLabels <- uniqLabels[which(!uniqLabels %in% removeIdx)]
  uniqLabels[ which(uniqLabels %in% removeIdx) ] <- NA
  newLabels <- match(res$label, uniqLabels)
  reassignIdx <- which(is.na(newLabels))
  reassignSsm <- ssm[reassignIdx, ]

  if (! is.null(res$full.model) ) {
    res$full.model$ccfMean <- t(as.matrix(res$full.model$ccfMean[-removeIdx]))
    res$full.model$ccfCov <- t(as.matrix(res$full.model$ccfCov[-removeIdx]))

    logRho <- res$full.model$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$full.model$logResponsibility <-
      if (length(res$label) ==length(res$full.model$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }

    res$full.model$responsibility <- exp(res$full.model$logResponsibility)
    if(!is.null(ssm)) {
      res$full.model$dirichletConcentration <- res$full.model$dirichletConcentration0 + colSums(res$full.model$responsibility)
      res$full.model <- VarationalExpectationStep(bn = ssm$var_counts,
                                       dn = ssm$ref_counts + ssm$var_counts,
                                       cn = ssm$normal_cn,
                                       cr = ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
                                         ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2),
                                       epi = epi,
                                       purity = unique(ssm$purity),
                                       model = res$full.model)

      res$full.model <- SortClusters(res$full.model)
      res$label <- apply(res$full.model$responsibility, 1, which.max)
      res$full.model$dirichletConcentration <- res$full.model$dirichletConcentration0 + colSums(res$full.model$responsibility)

      # if (is.null(dim(res$full.model$responsibility)) ) {
      #   res$label <- rep(1, length(res$full.model$responsibility))
      # } else {
      #   res$label <- apply(res$full.model$responsibility, 1, which.max)
      # }
    }



    if (!is.null(res$R)) {
      res$R=res$full.model$responsibility
    }

    if (!is.null(res$mu)) {
      res$mu=res$full.model$ccfMean
    }
    if (! is.null(res$full.model$Epi)) {

      Epi <- (res$full.model$dirichletConcentration + colSums(res$full.model$responsibility)) /
        (length(res$full.model$ccfMean) * res$full.model$dirichletConcentration0 + length(res$label))
      res$full.model$Epi <- Epi/sum(Epi)

      # if (is.null(dim(res$full.model$responsibility)) ) {
      #   res$full.model$Epi <- 1
      # } else {
      #   Epi <- (res$full.model$dirichletConcentration + colSums(res$full.model$responsibility)) /
      #     (length(res$full.model$ccfMean) * res$full.model$dirichletConcentration0 + length(res$label))
      #   res$full.model$Epi <- Epi/sum(Epi)
      # }
    }

  } else {

    res$ccfMean <- res$ccfMean[-removeIdx]
    res$ccfCov <- res$ccfCov[-removeIdx]

    logRho <- res$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$logResponsibility <-
      if (length(label) ==length(res$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }

    if(!is.null(ssm)) {
      res$dirichletConcentration <- res$dirichletConcentration0 + colSums(res$responsibility)
      res <- VarationalExpectationStep(bn = ssm$var_counts,
                                       dn = ssm$ref_counts + ssm$var_counts,
                                       cn = ssm$normal_cn,
                                       cr = ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
                                         ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2),
                                       epi = 1e-3,
                                       purity = unique(ssm$purity),
                                       model = res)
      res$dirichletConcentration <- res$dirichletConcentration0 + colSums(res$responsibility)
    }

  }

  res
}

#' Remove a (or more) cluster and reassign its data if the cluster is nonempty
#' @param res Ccube result list
#' @param removeIdx clusters to remove
#' @param ssm data
#' @param label assigned labels if res doesn't have label variable
#' @param tol stopping condition
#' @param maxiter maximum iteration
#' @param epi sequencing error
#' @param verbose show progress
#' @param fit_mult flag to estimate multiplicities
#' @param fit_hyper flag to estimate hyperparameters
#' @return Ccube result list
#' @export
RemoveClusterAndReassignVariantsWithEMsteps <- function(res, removeIdx, ssm = NULL, label = NULL, tol = 1e-8, maxiter = 100, epi = 1e-3, verbose = F,
                                                        fit_mult = T, fit_hyper = T) {

  if (length(removeIdx) == 0) {
    return(res)
  }

  if (!is.null(res$label)) {
    uniqLabels <- sort(unique(res$label))
  } else {
    uniqLabels <- sort(unique(label))
  }

  remainedLabels <- uniqLabels[which(!uniqLabels %in% removeIdx)]
  uniqLabels[ which(uniqLabels %in% removeIdx) ] <- NA
  newLabels <- match(res$label, uniqLabels)
  reassignIdx <- which(is.na(newLabels))
  reassignSsm <- ssm[reassignIdx, ]

  if (! is.null(res$full.model) ) {
    res$full.model$ccfMean <- t(as.matrix(res$full.model$ccfMean[-removeIdx]))
    res$full.model$ccfCov <-  t(as.matrix(res$full.model$ccfCov[-removeIdx]))

    logRho <- res$full.model$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$full.model$logResponsibility <-
      if (length(res$label) ==length(res$full.model$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }

    res$full.model$responsibility <- exp(res$full.model$logResponsibility)
    res$full.model$dirichletConcentration <- res$full.model$dirichletConcentration0 + colSums(res$full.model$responsibility)

    if(!is.null(ssm)) {
      ll = rep(-Inf, maxiter)
      vbiter = 1
      converged = F
      degenerated = F
      while (!converged & vbiter < maxiter & !degenerated) {
        vbiter = vbiter + 1

        res$full.model <- VarationalExpectationStep(bn = ssm$var_counts,
                                                    dn = ssm$ref_counts + ssm$var_counts,
                                                    cn = ssm$normal_cn,
                                                    cr = ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
                                                      ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2),
                                                    epi = epi,
                                                    purity = unique(ssm$purity),
                                                    model = res$full.model)

        res$full.model <- VariationalMaximimizationStep(bn = ssm$var_counts,
                                                        dn = ssm$ref_counts + ssm$var_counts,
                                                        cn = ssm$normal_cn,
                                                        cr = ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
                                                          ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2),
                                                        max_mult_cn_sub1 = ssm$major_cn_sub1,
                                                        max_mult_cn_sub2 = ssm$major_cn_sub2,
                                                        frac_cn_sub1 = ssm$frac_cn_sub1,
                                                        frac_cn_sub2 = ssm$frac_cn_sub2,
                                                        subclonal_cn = ssm$subclonal_cn,
                                                        epi = epi,
                                                        purity = unique(ssm$purity),
                                                        model = res$full.model,
                                                        fit_mult = fit_mult, fit_hyper = fit_hyper)

        ll[vbiter] = VariationalLowerBound(bn = ssm$var_counts,
                                   dn = ssm$ref_counts + ssm$var_counts,
                                   cn = ssm$normal_cn,
                                   cr = ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
                                     ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2),
                                   epi = epi,
                                   purity = unique(ssm$purity),
                                   model = res$full.model)/length(res$label)

        converged <- abs(ll[vbiter] - ll[vbiter-1]) < (tol * abs(ll[vbiter]))
        degenerated <- (ll[vbiter] - ll[vbiter-1]) < 0
        if(verbose) cat(sprintf("\rVB-EM-%d: L = %.8f \r", vbiter, ll[vbiter]))
      }

      res$full.model <- SortClusters(res$full.model)
      res$label <- apply(res$full.model$responsibility, 1, which.max)

      # if (is.null(dim(res$full.model$responsibility)) ) {
      #   res$label <- rep(1, length(res$full.model$responsibility))
      # } else {
      #   res$label <- apply(res$full.model$responsibility, 1, which.max)
      # }

    }



    if (!is.null(res$R)) {
      res$R=res$full.model$responsibility
    }

    if (!is.null(res$mu)) {
      res$mu=res$full.model$ccfMean
    }

    if (! is.null(res$full.model$Epi)) {

      Epi <- (res$full.model$dirichletConcentration + colSums(res$full.model$responsibility)) /
        (length(res$full.model$ccfMean) * res$full.model$dirichletConcentration0 + length(res$label))
      res$full.model$Epi <- Epi/sum(Epi)

      # if (is.null(dim(res$full.model$responsibility)) ) {
      #   res$full.model$Epi <- 1
      # } else {
      #   Epi <- (res$full.model$dirichletConcentration + colSums(res$full.model$responsibility)) /
      #     (length(res$full.model$ccfMean) * res$full.model$dirichletConcentration0 + length(res$label))
      #   res$full.model$Epi <- Epi/sum(Epi)
      # }


    }

  } else {

    res$ccfMean <- res$ccfMean[-removeIdx]
    res$ccfCov <- res$ccfCov[-removeIdx]

    logRho <- res$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$logResponsibility <-
      if (length(label) ==length(res$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }
    res$dirichletConcentration <- res$dirichletConcentration0 + colSums(res$responsibility)

    if(!is.null(ssm)) {
      res <- VarationalExpectationStep(bn = ssm$var_counts,
                                       dn = ssm$ref_counts + ssm$var_counts,
                                       cn = ssm$normal_cn,
                                       cr = ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
                                         ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2),
                                       epi = 1e-3,
                                       purity = unique(ssm$purity),
                                       model = res)
      res$dirichletConcentration <- res$dirichletConcentration0 + colSums(res$responsibility)
    }

  }

  res
}



#' Merge clusters
#' @param res ccube results list
#' @param ssm ccube data frame
#' @param tol stopping condition in VBEM
#' @param maxiter maximum iteration in VBEM
#' @param epi sequencing error
#' @return res ccube results list
#' @export
MergeClusters <- function(res = res, ssm = ssm, tol = 1e-8, maxiter = 100, epi = 1e-3) {

  res <- CullEmptyClusters(res = res, ssm = ssm)

  if (is.matrix(res$full.model$ccfMean)) {
    ccfCentersMap = res$full.model$ccfMean[,]
  } else {
    ccfCentersMap = res$full.model$ccfMean
  }

  clusterWeights <- as.data.frame(table(res$label), stringsAsFactors = F)
  clusterWeights <- dplyr::mutate(clusterWeights, Var1 = as.integer(Var1))
  ccfDistMat <- dist(ccfCentersMap)
  idx <- as.data.frame(which(as.matrix( ccfDistMat )< 0.1, arr.ind = TRUE), stringsAsFactors =F)
  idx <- dplyr::filter(idx, row != col)
  idx <- dplyr::left_join(idx, clusterWeights, by = c("col"="Var1"))
  idx <- dplyr::rename(idx, col_weights = Freq)
  idx <- dplyr::left_join(idx, clusterWeights, by = c("row"="Var1"))
  idx <- dplyr::rename(idx, row_weights = Freq)
  idx <- dplyr::mutate(rowwise(idx), remove_idx = c(row, col)[which.min(c(row_weights, col_weights))] )
  removeIdx <- unique(idx$remove_idx)
  res$mergeCluster <- length(removeIdx)>0
  return(RemoveClusterAndReassignVariantsWithEMsteps(res = res,
                                                     removeIdx = removeIdx, ssm = ssm,
                                                     tol = tol, maxiter = maxiter, epi = epi))
}


#' Write files in PCAWG-11 formats, works for both CcubeCore and ccube_m6 output
#' @param ssm data
#' @param res Ccube result list
#' @param resultFolder path to file
#' @param sampleName sample name
#' @param allFormats all PCAWG output formats
#' @param basicFormats basic PCAWG output formats
#' @param outputMult format for multiplicities
#' @param outputAssignProb format for assignment probabilities
#' @param outputAssign format for assignment
#' @param outputSubStruct format for subclonal structure
#' @param outputCcm format for co-clustering matrix
#' @param outputCcmIdx format for variant indices in co-clustering matrix
#' @return NULL
#' @export
WritePcawgFormats <- function(ssm, res, resultFolder, sampleName,
                              allFormats = F, basicFormats = T,
                              outputMult = F, outputAssignProb = F, outputAssign = F,
                              outputSubStruct = F, outputCcm = F, outputCcmIdx = F) {
  if (basicFormats) {
    outputMult <- T
    outputAssign <- T
    outputSubStruct <- T
  }

  if (allFormats) {
    outputMult <- T
    outputAssignProb <- T
    outputAssign <- T
    outputSubStruct <- T
    outputCcm <- T
    outputCcmIdx <- T
  }

  ## output calibration format
  uniqLabels <- unique(res$label)
  dir.create(resultFolder, recursive = T)

  id <- do.call(rbind, strsplit(as.character(ssm$mutation_id), "_", fixed = T))

  # Multiplicity
  if (outputMult) {
    mult <- data.frame(chr = id[,1], pos = id[,2])
    mult$tumour_copynumber <- ssm$major_cn+ssm$minor_cn
    mult$multiplicity <- ssm$ccube_mult
    fn <- paste0(resultFolder, "/",
                 sampleName, "_multiplicity.txt")
    write.table(mult, file = fn, sep = "\t", row.names = F, quote = F)
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
    rm(mult)
  }


  # Assignment Prob
  if(outputAssignProb){
    mutAssign <- data.frame(chr = id[,1], pos = id[,2])

    if (length(uniqLabels) == 1) {
      mutR = data.frame(res$full.model$responsibility)
      colnames(mutR) <- "cluster_1"
    } else {
      mutR <- data.frame(res$full.model$responsibility[, sort(uniqLabels)])
      colnames(mutR) <- paste0("cluster_", seq_along(uniqLabels))
    }

    mutAssign <- data.frame(mutAssign, mutR)
    fn <- paste0(resultFolder, "/",
                 sampleName, "_assignment_probability_table.txt")
    write.table(mutAssign, file = fn, sep = "\t", row.names = F, quote = F)
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
    rm(mutR, mutAssign)
  }

  # Assignment
  if (outputAssign) {
    mutAssign <- data.frame(chr = id[,1], pos = id[,2])
    if (length(uniqLabels) == 1) {
      mutAssign$cluster = 1
    } else {
      mutAssign$cluster <- res$label
    }
    fn <- paste0(resultFolder, "/",
                 sampleName, "_mutation_assignments.txt")
    write.table(mutAssign, file = fn, sep = "\t", row.names = F, quote = F)
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
    rm(mutAssign)
  }

  # subclonal structure
  if (outputSubStruct) {
    cellularity <- unique(ssm$purity)
    clusterCertainty <- as.data.frame(table(res$label), stringsAsFactors = F)
    clusterCertainty <- dplyr::rename(clusterCertainty, cluster = Var1, n_ssms = Freq)
    clusterCertainty$proportion <- res$full.model$ccfMean[as.integer(clusterCertainty$cluster)] * cellularity
    clusterCertainty$cluster <- seq_along(uniqLabels)
    fn <- paste0(resultFolder, "/",
                 sampleName, "_subclonal_structure.txt")
    write.table(clusterCertainty, file = fn, sep = "\t", row.names = F, quote = F)
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
  }


  # ccm
  if (outputCcm) {
    coClustMat <- Matrix::tcrossprod(res$full.model$responsibility)
    diag(coClustMat) <- 1
    fn = paste0(resultFolder, "/", sampleName, "_coassignment_probabilities.txt")
    write.table(coClustMat, file = fn, sep = "\t", col.names = F, row.names = F, quote = F )
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
    rm(coClustMat)
  }


  # ccm index
  if (outputCcmIdx) {
    indexFile <- cbind(id[,1], id[, 2], seq_along(id[,1]))
    colnames(indexFile) <- c("chr", "pos", "col")
    fn = paste0(resultFolder, "/", sampleName, "_index.txt")
    write.table(indexFile, file = fn, sep = "\t", row.names = F, quote = F )
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
  }

}


CheckAndPrepareCcubeInupts <- function(mydata, estimatePurity = T) {

  stopifnot(
    all(c("var_counts","ref_counts") %in% names(mydata)))

  if ( "major_cn_sub1" %in% names(mydata) &
       "minor_cn_sub1" %in% names(mydata) &
       "major_cn_sub2" %in% names(mydata) &
       "minor_cn_sub2" %in% names(mydata) &
       "subclonal_cn" %in% names(mydata) &
       "frac_cn_sub1" %in% names(mydata) &
       "frac_cn_sub2" %in% names(mydata) ) {
    mydata$subclonal_cn <- mydata$frac_cn_sub1 < 1
    return(mydata)
  }

  if ( ! "frac_cn_sub1" %in% names(mydata) ) {
    message(sprintf("Missing column: frac_cn_sub1. Assuming input copy number profiles are clonal"))
    mydata$frac_cn_sub1 <- 1
  }

  if ( ! "frac_cn_sub2" %in% names(mydata) ) {
    message(sprintf("Missing column: frac_cn_sub2. Set frac_cn_sub2 as 1- frac_cn_sub1"))
    mydata$frac_cn_sub2 <- 1 - mydata$frac_cn_sub1
  }

  if ( ! "subclonal_cn" %in% names(mydata) ) {
    message(sprintf("Missing column: subclonal_cn. Set subclonal_cn as frac_cn_sub1 < 1"))
    mydata$subclonal_cn <- mydata$frac_cn_sub1 < 1
  }

  if ( ! "major_cn_sub1" %in% names(mydata) ) {
    message(sprintf("Missing column: major_cn_sub1. Assuming input copy number profiles are clonal, set major_cn_sub1 as major_cn"))
    stopifnot(all(c("major_cn") %in% names(mydata)))
    mydata <- dplyr::rename(mydata, major_cn_sub1 = major_cn)
  }

  if ( ! "minor_cn_sub1" %in% names(mydata) ) {
    message(sprintf("Missing column: minor_cn_sub1 Assuming input copy number profiles are clonal, set minor_cn_sub1 as minor_cn"))
    stopifnot(all(c("minor_cn") %in% names(mydata)))
    mydata <- dplyr::rename(mydata, minor_cn_sub1 = minor_cn)
  }

  if ( ! "major_cn_sub2" %in% names(mydata) ) {
    message(sprintf("Missing column: major_cn_sub2 Assuming input copy number profiles are clonal, set major_cn_sub2 as -100"))
    mydata$major_cn_sub2 <- -100
  }

  if ( ! "minor_cn_sub2" %in% names(mydata) ) {
    message(sprintf("Missing column: minor_cn_sub2 Assuming input copy number profiles are clonal, set minor_cn_sub2 as -100"))
    mydata$minor_cn_sub2 <- -100
  }

  mydata <- dplyr::mutate(mydata,
                          major_cn = frac_cn_sub1 * major_cn_sub1 + frac_cn_sub2 * major_cn_sub2,
                          minor_cn = frac_cn_sub1 * minor_cn_sub1 + frac_cn_sub2 * minor_cn_sub2,
                          total_cn = major_cn + minor_cn)

  if ( ! "normal_cn" %in% names(mydata) ) {
    message(sprintf("Missing column: normal_cn. Set normal_cn as 2"))
    mydata$normal_cn <- 2
  }

  if ( ! "purity" %in% names(mydata) & estimatePurity) {
    message(sprintf("Missing column: purity. Estimate purity with GetPurity"))
    mydata$purity <- GetPurity(mydata)
  }

  return(mydata)
}


#' Computing responsibilities (assignment probabilities) with modified Ccube VBEM-step
#' @param res A reference Ccube results list
#' @param ssm A data frame of SNVs to be assigned. Ideally, the data has been processed by Ccube model. So it should have the ccube_mult column.
#' @param tol convergence threshold
#' @param maxiter maximum iterations, default number is 100.
#' @param epi sequencing error, default is 1e-3
#' @param verbose show VBEM progress
#' @return a standard Ccube model with recomputed responsibilities and logResponsibilities
#' @export
AssignWithCcube <- function(res, ssm, tol = 1e-8, maxiter = 100, epi = 1e-3, verbose = F) {

  ssm <- CheckAndPrepareCcubeInupts(ssm)


  if ("ccube_mult" %in% names(ssm) ) {
    res$full.model$bv <- ssm$ccube_mult
  } else {
    ssm <- GetCcf(ssm, use = "use_base")
    res$full.model$bv <- ssm$rough_mult
  }

  # res$full.model$dirichletConcentration0 <- mean(res$full.model$dirichletConcentration/sum(res$full.model$dirichletConcentration))
  res$full.model$ccfMean <- t(as.matrix(as.numeric(res$full.model$ccfMean)))
  res$full.model$ccfCov <-  t(as.matrix(as.numeric(res$full.model$ccfCov)))

  ll = rep(-Inf, maxiter)
  vbiter = 1
  converged = F
  degenerated = F

  while (!converged & vbiter < maxiter & !degenerated) {

    vbiter = vbiter + 1

    if (vbiter == 2)  {
      no.weights = T
    } else {
      no.weights = F
    }

    res$full.model <- VarationalExpectationStep(bn = ssm$var_counts,
                                                dn = ssm$ref_counts + ssm$var_counts,
                                                cn = ssm$normal_cn,
                                                cr = ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
                                                  ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2),
                                                epi = epi,
                                                purity = unique(ssm$purity),
                                                model = res$full.model,
                                                no.weights = no.weights)

    res$full.model <- VariationalMaximimizationStep(bn = ssm$var_counts,
                                                    dn = ssm$ref_counts + ssm$var_counts,
                                                    cn = ssm$normal_cn,
                                                    cr = ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
                                                      ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2),
                                                    max_mult_cn_sub1 = ssm$major_cn_sub1,
                                                    max_mult_cn_sub2 = ssm$major_cn_sub2,
                                                    frac_cn_sub1 = ssm$frac_cn_sub1,
                                                    frac_cn_sub2 = ssm$frac_cn_sub2,
                                                    subclonal_cn = ssm$subclonal_cn,
                                                    epi = epi,
                                                    purity = unique(ssm$purity),
                                                    model = res$full.model,
                                                    fit_params_ccf = F)

    ll[vbiter] = VariationalLowerBound(bn = ssm$var_counts,
                                       dn = ssm$ref_counts + ssm$var_counts,
                                       cn = ssm$normal_cn,
                                       cr = ssm$frac_cn_sub1 * (ssm$major_cn_sub1 + ssm$minor_cn_sub1) +
                                         ssm$frac_cn_sub2 *(ssm$major_cn_sub2 + ssm$minor_cn_sub2),
                                       epi = epi,
                                       purity = unique(ssm$purity),
                                       model = res$full.model)/length(res$full.model$bv)

    converged <- abs(ll[vbiter] - ll[vbiter-1]) < (tol * abs(ll[vbiter]))
    degenerated <- (ll[vbiter] - ll[vbiter-1]) < 0
    if(verbose) cat(sprintf("\rVB-EM-%d: L = %.8f \r", vbiter, ll[vbiter]))
  }

  res$full.model <- SortClusters(res$full.model)

  res$label <- apply(res$full.model$responsibility, 1, which.max)

  Epi <- (res$full.model$dirichletConcentration + colSums(res$full.model$responsibility)) /
    (length(res$full.model$ccfMean) * res$full.model$dirichletConcentration0 + length(res$label))
  res$full.model$Epi <- Epi/sum(Epi)

  return(res)
}

#' Annotate mydata with Ccube results
#' @param ssm Ccube input data
#' @param res Ccube results list
#' @return ssm Annotated input data
#' @export
AnnotateCcubeResults <- function(ssm, res) {

  ssm <- CheckAndPrepareCcubeInupts(ssm)
  ssm <- GetCcf(ssm, use = "use_base")


  ssm$ccube_ccf_mean <- res$full.model$ccfMean[res$label]
  ssm$ccube_mult <- res$full.model$bv

  ssm <- dplyr::mutate( dplyr::rowwise(ssm),
                        ccube_ccf = MapVaf2CcfPyClone(vaf,
                                                      purity,
                                                      normal_cn,
                                                      total_cn,
                                                      total_cn,
                                                      ccube_mult,
                                                      constraint=F) )

  return(ssm)
}


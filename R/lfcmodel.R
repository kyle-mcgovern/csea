verify <- function(m, ...) {
	UseMethod("verify", m)
}

model_sample <- function(m, ...) {
	UseMethod("model_sample", m)
}

model_sample_parallel <- function(m, ...) {
	UseMethod("model_sample_parallel", m)
}

new_seainput <- function(obs_lfc, shuffled_lfc) {
	seainput <- structure(
		list(
			obs_lfc=obs_lfc,
			shuffled_lfc=shuffled_lfc
		),
		class=c("seainput")
	)
	verify(seainput)
	seainput
}

verify.seainput <- function(object, sample_type=NULL) {
	if(!is.null(sample_type)) {
		if(!(sample_type%in%c("column", "entity"))) {
			stop("sample_type must be 'column' or 'entity'")
		}
		if(sample_type=="column") {
			if(!is.matrix(object$shuffled_lfc)) {
				stop("If sample_type is 'column', object$shuffled_lfc must be a matrix")
			}
			if(!identical(dim(object$obs_lfc), dim(object$shuffled_lfc))) {
				stop("dim(object$obs_lfc) must equal dim(object$shuffled_lfc)")
			}
		}
	}
	if(!is.matrix(object$obs_lfc)) {
		stop("object$obs_lfc must be a matrix")
	}
}

#' Construct a lfcmodel object
#' 
#' @param D matrix (M x N) of counts
#' @param X matrix (N x Q) design matrix
#' @param lfc_column integer representing the column of X to shuffle for permutations
#' @param prior Either a vector (length M) or a matrix (M x N) representing dirichlet priors
#' @return object of class lfcmodel
new_lfcmodel <- function(D, X, lfc_column, prior=NULL) {
	if(is.null(prior)) {
		prior <- rep(1, nrow(D))
	}
	lfcmodel <- structure(
		  list(
		       # D is a matrix of counts of size mxn
		       D=D,
		       # Design matrix
		       X=X,
		       # Column labels is a factor of length n
		       lfc_column=lfc_column,
		       # Prior is a vector of length m
		       prior=prior
		  ),
		  class=c("lfcmodel")
	)
	verify(lfcmodel)
	lfcmodel
}

verify.lfcmodel <- function(m) {
	if(!is.matrix(m$D)||!is.matrix(X)) {
		stop("X and D must both be matrices")
	}
	if(ncol(m$D)!=nrow(m$X)) {
		stop("ncol(D) must equal nrow(X)")
	}
	if(m$lfc_column>ncol(m$X)) {
		stop("lfc_column larger than ncol(X)")
	}
	if(!all(m$X[,m$lfc_column]%in%c(0,1))) {
		stop("Values in X[,lfc_column] must be 0 or 1 only")
	}
	if(is.matrix(m$prior)) {
		if(!identical(dim(m$prior), dim(m$D))) {
			stop("if prior is a matrix, dim(prior) must equal dim(D)")
		}
	} else {
		if(length(m$prior)!=nrow(m$D)) {
			stop("if prior is a vector, length(prior) must equal nrow(D)")
		}
	}
}

model_sample.lfcmodel <- function(m, sample_type="column") {
	# Multinomial-Dirichlet Posterior Sample
	p_sample <- rdirichlet_cpp(m$D+m$prior)

	# ALR coefficitients
	alr_sample <- log(t(t(p_sample)/p_sample[nrow(p_sample),]))[1:(nrow(p_sample)-1),]

	# Linear Model LFC
	B <- ginv(t(m$X)%*%m$X)%*%t(m$X)%*%t(alr_sample)
	obs_lfc <- t(B)[,m$lfc_column]
	obs_lfc <- c(obs_lfc, 0)
	obs_lfc <- exp(obs_lfc)/sum(exp(obs_lfc))
	obs_lfc <- log(obs_lfc/gm_mean(obs_lfc))

	# Linear Model Shuffled Columns
	if (sample_type=="column") {
		X <- m$X
		X[,m$lfc_column] <- sample(X[,m$lfc_column], replace=F)
		B <- ginv(t(X)%*%X)%*%t(X)%*%t(alr_sample)
		shuffled_lfc <- t(B)[,m$lfc_column]
		shuffled_lfc <- c(shuffled_lfc, 0)
		shuffled_lfc <- exp(shuffled_lfc)/sum(exp(shuffled_lfc))
		shuffled_lfc <- log(shuffled_lfc/gm_mean(shuffled_lfc))

	} else if (sample_type=="entity") {
		shuffled_lfc <- NULL
	} else {
		stop("sample_type must be 'column' or 'entity'")
	}
	return(list(obs_lfc=obs_lfc, shuffled_lfc=shuffled_lfc))
}

model_sample_parallel.lfcmodel <- function(m, permutations,
					   sample_type="column",
					   cores=detectCores()) {
	# Set up parallel processing
	cl <- makeCluster(cores)
	registerDoSNOW(cl)

	# Set up progress bar
	pb <- txtProgressBar(max=permutations, style=3)
	progress <- function(n) setTxtProgressBar(pb, n)
	opts <- list(progress = progress)

	comb <- function(results, r) {
		results$obs_lfc <- cbind(results$obs_lfc, r$obs_lfc)
		if (is.null(r$shuffled_lfc)) {
			results$shuffled_lfc <- NULL
		} else {
			results$shuffled_lfc <- cbind(results$shuffled_lfc,
						      r$shuffled_lfc)
		}
		results
	}

	# Take samples in parallel
	res <- foreach(i=1:permutations,
		       .combine=combine_lfcmodel_sample,
		       .options.snow=opts,
		       .export=c("model_sample",
				 "model_sample.lfcmodel",
				 "gm_mean")) %dopar% {
		s <- model_sample(m, sample_type)
		s
	}
	close(pb)
	stopCluster(cl)
	gc()
	res <- new_seainput(res$obs_lfc, res$shuffled_lfc)
	verify(res, sample_type)
	res
}


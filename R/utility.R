gm_mean = function(x, na.rm=TRUE){
	exp(sum(log(x)) / length(x))
}

build_txt_pb_opts <- function(iterations, style=3) {
	# Set up progress bar
	pb <- txtProgressBar(max=iterations, style=style)
	progress <- function(n) setTxtProgressBar(pb, n)
	opts <- list(progress=progress)
	list(opts=opts, pb=pb)
}

combine_lfcmodel_sample <- function(results, r) {
	results$obs_lfc <- cbind(results$obs_lfc, r$obs_lfc)
	if (is.null(r$shuffled_lfc)) {
		results$shuffled_lfc <- NULL
	} else {
		results$shuffled_lfc <- cbind(results$shuffled_lfc,
					      r$shuffled_lfc)
	}
	results
}


combine_csea_lfcmodel <- function(results, r) {
	results$shuffled_scores <- cbind(results$shuffled_scores,
					 r$shuffled_scores)
	results$obs_scores <- cbind(results$obs_scores, r$obs_scores)
	results
}


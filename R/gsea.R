gsea_base <- function(scores, inds) {
	sorted_obj <- sort(scores, decreasing=T, index.return=T)
	inds <- which(sorted_obj$ix%in%inds)
	inds <- sort(inds) 
	scores <- sorted_obj$x
	num_gt <- 0
	hit_scores <- abs(scores[inds])/sum(abs(scores[inds])) 
	return(gsea_score_C(hit_scores, inds, length(scores)) )
}

gsea <- function(object, entity_sets, iterations=2000, cores=detectCores()) {
	cl <- makeCluster(cores)
	registerDoSNOW(cl)
	pbd <- build_txt_pb_opts(iterations)
	output <- foreach(i=1:iterations, .options.snow=pbd$opts, .export=c("gsea_base"),
		       .combine=combine_csea_lfcmodel) %dopar% {
		# Take sample from model
		ms <- model_sample(object, "column") # FIXME
		# Run through all gene sets
		obs_scores <- c()
		shuffled_scores <- c()
		for(j in 1:length(entity_sets)) {
			obs_scores <- c(obs_scores, gsea_base(ms$obs_lfc, entity_sets[[j]]))
			shuffled_scores <- c(shuffled_scores, gsea_base(ms$shuffled_lfc, entity_sets[[j]]))
		}
		list(shuffled_scores=shuffled_scores, obs_scores=obs_scores)
	}
	res <- csea_lfcmodel_agg(output, names(entity_sets))
	close(pbd$pb)
	stopCluster(cl)
	gc()
	res

}

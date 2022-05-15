calculate_W <- function(Y, X, lfc_column, epsilon) {
	W_par <- t(t(Y)/colSums(Y))
	W_perp <- log(apply(W_par, 2, gm_mean))
	error <- epsilon * X[,lfc_column]
	W <- log(W_par)
	W <- sweep(W, 2, W_perp, "+")
	W <- sweep(W, 2, error, "+")
	W
}

ols_solution <- function(Y, X, lfc_column) {
	t(ginv(t(X)%*%X)%*%t(X)%*%t(Y))[,lfc_column]
}

scale_sensitivity_base <- function(Y, X, ind, lfc_column, pseudo_count, epsilons, iterations=2000) {
	# Loop over values of epsilon
	Y <- Y + pseudo_count
	
	p_values <- c()
	obs_scores <- c()
	for(epsilon in epsilons) {
		# Compute observed LFC gsea score
		W <- calculate_W(Y, X, lfc_column, epsilon)
		obs_lfc <- ols_solution(W, X, lfc_column)
		obs_score <- gsea_base(obs_lfc, ind)
		obs_scores <- c(obs_scores, obs_score)
		# Compute permuted scores
		shuffled_scores <- c()
		s_X <- X
		for(i in 1:iterations) {
			s_X[,lfc_column] <- sample(s_X[,lfc_column], replace=F)
			shuffled_lfc <- ols_solution(W, s_X, lfc_column)
			shuffled_scores <- c(shuffled_scores, gsea_base(shuffled_lfc, ind))
		}
		# Compute p-value
		shuffled_scores <- shuffled_scores[sign(shuffled_scores)==sign(obs_score)]
		p_value <- sum(abs(obs_score) <= abs(shuffled_scores)) / length(shuffled_scores)
		p_values <- c(p_values, p_value)
	}
	return(list(obs_scores=obs_scores, p_values=p_values))
}

#scale_sensitivity_analysis
scale_sensitivity_analysis <- function(Y, X, entity_sets, lfc_column, pseudo_count, epsilons, iterations=2000, cores=detectCores()) {
	# Set up parallel processing
	cl <- makeCluster(cores)
	registerDoSNOW(cl)

	combn <- function(results, r) {
		results$obs_scores <- cbind(results$obs_scores, r$obs_scores)
		results$p_values <- cbind(results$p_values, r$p_values)
		results
	}

	pbd <- build_txt_pb_opts(length(entity_sets))
	res <- foreach(i=1:length(entity_sets), .options.snow=pbd$opts,
		       .export=c("scale_sensitivity_base"), .combine=combn) %dopar% {
		ss_res <- scale_sensitivity_base(Y, X, entity_sets[[i]], lfc_column, pseudo_count, epsilons, iterations)
		ss_res
	}
	close(pbd$pb)
	stopCluster(cl)
	gc()
	res
}

sst <- function(object) {
	;
}

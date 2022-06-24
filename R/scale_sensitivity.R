calculate_W <- function(Y, X, lfc_column, epsilon) {
	W_par <- t(t(Y)/colSums(Y))
	W_perp <- log(apply(W_par, 2, function(col) 1/gm_mean(col)))
	error <- epsilon * X[,lfc_column]
	W <- log(W_par)
	W <- sweep(W, 2, W_perp, "+")
	W <- sweep(W, 2, error, "+")
	W
}

ols_solution <- function(Y, X, lfc_column) {
	t(ginv(t(X)%*%X)%*%t(X)%*%t(Y))[,lfc_column]
}

clr_inv <- function(x) {
  exp(x)/sum(exp(x))
}

gsea_compositional_weighting <- function(scores, inds) {
  weights <- c()
  for(ind in inds) {
    weights <- c(weights, (length(inds)-1) * sum(abs(scores - scores[ind])))
  }
  weights/sum(weights)
}

gsea_base_inds <- function(scores, inds, cw=F) {
	sorted_obj <- sort(scores, decreasing=T, index.return=T)
	scores <- sorted_obj$x
	gsea_scores <- rep(0, length(inds))
	for (i in 1:length(inds)) {
		idx <- sort(which(sorted_obj$ix%in%(inds[[i]])))
		if(cw) {
		  hit_scores <- gsea_compositional_weighting(scores, idx)
		} else {
		  hit_scores <- abs(scores[idx])/sum(abs(scores[idx]))
		}
		gsea_scores[i] <- gsea_score_C(hit_scores, idx, length(scores))
	}
	gsea_scores
}

add_to_totals <- function(total, gt, obs_scores, shuffled_scores) {
  obs_sign <- sign(obs_scores)
  shuffled_sign <- sign(shuffled_scores)
  total <- total + as.numeric(shuffled_sign==obs_sign)
  gt <- gt + as.numeric( (shuffled_sign==obs_sign) & (abs(obs_scores) <= abs(shuffled_scores)) )
  return(list(total=total, gt=gt))
}

draw_shuffled_col <- function(W, X, lfc_column) {
  X[,lfc_column] <- sample(X[,lfc_column], replace=F)
  ols_solution(W, X, lfc_column)
}

draw_perm_col <- function(W, X, X_shuf, lfc_column) {
  X[,lfc_column] <- X_shuf
  ols_solution(W, X, lfc_column)
}

scale_sensitivity_analysis <- function(Y, X, inds, lfc_column, pseudo_count, epsilons, iterations=2000, cores=detectCores(), cw=F) {
	Y <- Y + pseudo_count
	
	cl <- makeCluster(cores)
	registerDoSNOW(cl)

	combn <- function(results, r) {
		results$obs_scores <- cbind(results$obs_scores, r$obs_scores)
		results$p_values <- cbind(results$p_values, r$p_values)
		results
	}

	pbd <- build_txt_pb_opts(length(epsilons))
  res <- foreach(j=1:length(epsilons), .options.snow=pbd$opts,
		.export=c("gsea_base_inds", "calculate_W",
		          "ols_solution", "draw_shuffled_col"), .combine=combn) %dopar% {
		epsilon <- epsilons[j]
		# Calculate observed
		W <- calculate_W(Y, X, lfc_column, epsilon)
		obs_lfc <- ols_solution(W, X, lfc_column)
		obs_scores <- gsea_base_inds(obs_lfc, inds, cw)

		gt <- rep(0, length(obs_scores))
		total <- rep(0, length(obs_scores))
		for(i in 1:iterations) {
			shuffled_lfc <- draw_shuffled_col(W, X, lfc_column)
			shuffled_scores <- gsea_base_inds(shuffled_lfc, inds, cw)
			r <- add_to_totals(total, gt, obs_scores, shuffled_scores)
			total <- r$total
			gt <- r$gt
		}
		p_values <- gt / total
		list(obs_scores=obs_scores, p_values=p_values)
	}
	res
}

gsea_parallel_matrix <- function(S, X, inds, lfc_column, iterations=2000, permutation_matrix=NULL, cores=detectCores(), cw=F) {
	cl <- makeCluster(cores)
	registerDoSNOW(cl)

	# Calculate observed
	obs_lfc <- ols_solution(S, X, lfc_column)
	obs_scores <- gsea_base_inds(obs_lfc, inds, cw)
	
	if(!is.null(permutation_matrix)) {
	  iterations <- ncol(permutation_matrix) 
	}
	
	pbd <- build_txt_pb_opts(iterations)
	res <- foreach(j=1:iterations, .options.snow=pbd$opts,
	               .export=c("gsea_base_inds", "calculate_W",
	                         "ols_solution", "draw_shuffled_col"), .combine=rbind) %dopar% {
	  if(is.null(permutation_matrix)) {
	    shuffled_lfc <- draw_shuffled_col(S, X, lfc_column)
	  } else {
	    shuffled_lfc <- draw_perm_col(S, X, permutation_matrix[,j], lfc_column) 
	  }
		shuffled_scores <- gsea_base_inds(shuffled_lfc, inds, cw)
		shuffled_scores
	}
	
	p_values <- c()
	for(i in 1:ncol(res)) {
	  shuff_filt <- res[,i][sign(res[,i])==sign(obs_scores[i])]
	  p_values <- c(p_values, sum(abs(shuff_filt) >= abs(obs_scores[i]))/length(shuff_filt))
	}

	close(pbd$pb)
	stopCluster(cl)
	gc()

	return(list(obs_scores=obs_scores, p_values=p_values))
}

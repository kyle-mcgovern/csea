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

gsea_base_inds <- function(scores, inds) {
	sorted_obj <- sort(scores, decreasing=T, index.return=T)
	scores <- sorted_obj$x
	gsea_scores <- rep(0, length(inds))
	for (i in 1:length(inds)) {
		idx <- sort(which(sorted_obj$ix%in%(inds[[i]])))
		hit_scores <- abs(scores[idx])/sum(abs(scores[idx]))
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

scale_sensitivity_analysis <- function(Y, X, inds, lfc_column, pseudo_count, epsilons, iterations=2000, cores=detectCores()) {
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
		obs_scores <- gsea_base_inds(obs_lfc, inds)

		gt <- rep(0, length(obs_scores))
		total <- rep(0, length(obs_scores))
		for(i in 1:iterations) {
			shuffled_lfc <- draw_shuffled_col(W, X, lfc_column)
			shuffled_scores <- gsea_base_inds(shuffled_lfc, inds)
			r <- add_to_totals(total, gt, obs_scores, shuffled_scores)
			total <- r$total
			gt <- r$gt
		}
		p_values <- gt / total
		list(obs_scores=obs_scores, p_values=p_values)
	}
	res
}


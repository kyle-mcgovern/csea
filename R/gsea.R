gsea <- function(scores, inds, iters=2000) {
	sorted_obj <- sort(scores, decreasing=T, index.return=T)
	inds <- which(sorted_obj$ix%in%inds)
	inds <- sort(inds) 
	scores <- sorted_obj$x
	num_gt <- 0
	hit_scores <- abs(scores[inds])/sum(abs(scores[inds])) 
	obs_score <- gsea_score_C(hit_scores, inds, length(scores)) 
	found <- 0
	for(i in 1:iters) {
		random_inds <- sample(1:length(scores), length(inds), replace=F)
		random_inds <- sort(random_inds) 
		hit_scores <- abs(scores[inds])/sum(abs(scores[inds])) 
		rand_score <- gsea_score_C(hit_scores, random_inds, length(scores)) 
		if (obs_score < 0 & rand_score > 0) {
			next
		} else if (obs_score > 0 & rand_score < 0) {
			next
		} else {
			found <- found + 1
		}
		if (abs(rand_score) >= abs(obs_score)) {
			num_gt <- num_gt + 1
		}
	}
	pval <- num_gt/found
	if(is.na(pval)) { # Only occurs if found and num_gt are 0
		pval <- 0
	}
	return(list(obs_stat=obs_score, pval=pval))
}

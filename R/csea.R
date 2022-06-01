csea_base <- function(s, i, statistic, power=1) {
	in_set <- s[i]
	not_in_set <- s[-i]
	if(statistic=="cvm") {
		return(cvm_stat(in_set, not_in_set, power))
	} else if (statistic=="ks") {
		return(ks_stat(in_set, not_in_set, power))
	} else if (statistic=="ad") {
		return(ad_stat(in_set, not_in_set, power))
	} else if (statistic=="dts") {
		return(dts_stat(in_set, not_in_set, power))
	} else {
		stop(paste0(statistic, " not a valid statistic. Must be ks, cvm, ad, dts."))
	}
}

csea <- function(object, entity_set, perm_type, statistic, power=1, iterations=2000) {
	if (!(perm_type%in%c("entity", "column"))) {
		stop("perm_type must be entity or column")
	}
	res <- list(obs_scores=c(), shuffled_scores=c())
	if(is(object, "seainput")) {
		for(i in 1:ncol(object$obs_lfc)) {
			obs_score <- csea_base(object$obs_lfc[,i], entity_set, statistic, power)
			if (perm_type=="column") {
				shuffled_score <- csea_base(object$shuffled_lfc[,i], entity_set, statistic, power)
			} else {
				shuffled_set <- sample(1:length(object$obs_lfc[,i]), length(entity_set), replace=F)
				shuffled_score <- csea_base(object$obs_lfc[,i], shuffled_set, statistic, power)
			}
			res$obs_scores <- c(res$obs_scores, obs_score)
			res$shuffled_scores <- c(res$shuffled_scores, shuffled_score)
		}
		res$p_value <- sum(res$obs_scores <= res$shuffled_scores) / length(res$shuffled_scores)
	} else if(is.vector(object) & !is.list(object)) {
		res$obs_scores <- csea_base(object, entity_set, statistic, power)
		for(i in 1:iterations) {
			shuffled_set <- sample(1:length(object), length(entity_set), replace=F)
			res$shuffled_scores <- c(res$shuffled_scores, csea_base(object, shuffled_set, statistic, power))
		}
		res$p_value <- sum(res$obs_scores <= res$shuffled_scores) / length(res$shuffled_scores)

	} else {
		stop("object must be seainput or vector")
	}
	return(res)
}

csea_lfcmodel_agg <- function(output, set_names) {
	res <- c()
	p_values <- c()
	for(i in 1:nrow(output$obs_scores)) {
		p_values <- c(p_values, sum(output$obs_scores[i,] <= output$shuffled_scores[i,])/ncol(output$obs_scores))
	}
	res <- cbind(res, p_values)
	res <- cbind(res, apply(output$obs_scores, 1, mean))
	res <- cbind(res, apply(output$obs_scores, 1, sd))
	res <- cbind(res, apply(output$shuffled_scores, 1, mean))
	res <- cbind(res, apply(output$shuffled_scores, 1, sd))
	colnames(res) <- c("p_value", "mean_obs_score", "sd_obs_scores",
			   "mean_shuffled_scores", "sd_shuffled_scores")
	row.names(res) <- set_names
	res
}

csea_parallel <- function(object, entity_sets, perm_type, statistic, cores=detectCores(), power=1, iterations=2000) {
	# Set up parallel processing
	cl <- makeCluster(cores)
	registerDoSNOW(cl)

	# Handle different input types
	if(is(object, "seainput")) {
		pbd <- build_txt_pb_opts(length(entity_sets))
		res <- foreach(i=1:length(entity_sets), .options.snow=pbd$opts, .export=c("csea"), .combine=rbind) %dopar% {
			res <- csea(object, entity_sets[[i]], perm_type, statistic, power, iterations)
			res$pathway_name <- names(entity_sets)[i]
			c(res$p_value, mean(res$obs_scores), sd(res$obs_scores), mean(res$shuffled_scores), sd(res$shuffled_scores))
		}
		colnames(res) <- c("p_value", "mean_obs_score", "sd_obs_scores",
				   "mean_shuffled_scores", "sd_shuffled_scores")
		row.names(res) <- names(entity_sets)
	} else if(is(object, "lfcmodel")) {
		pbd <- build_txt_pb_opts(iterations)
		output <- foreach(i=1:iterations, .options.snow=pbd$opts, .export=c("csea"),
			       .combine=combine_csea_lfcmodel) %dopar% {
			# Take sample from model
			ms <- model_sample(object, perm_type)
			# Run through all gene sets
			obs_scores <- c()
			shuffled_scores <- c()
			for(j in 1:length(entity_sets)) {
				if(perm_type=="column") {
					obs_scores <- c(obs_scores, csea_base(ms$obs_lfc, entity_sets[[j]], statistic, power))
					shuffled_scores <- c(shuffled_scores, csea_base(ms$shuffled_lfc, entity_sets[[j]], statistic, power))
				} else {
					obs_scores <- c(obs_scores, csea_base(ms$obs_lfc, entity_sets[[j]], statistic, power))
					shuffled_set <- sample(1:length(ms$obs_lfc), length(entity_sets[[j]]), replace=F)
					shuffled_scores <- c(shuffled_scores, csea_base(ms$obs_lfc, shuffled_set, statistic, power))
				}

			}
			list(shuffled_scores=shuffled_scores, obs_scores=obs_scores)
		}
		res <- csea_lfcmodel_agg(output, names(entity_sets))
	} else if(is.vector(object) & !is.list(object)) {
		pbd <- build_txt_pb_opts(length(object))
		if(perm_type!="entity") {
			stop("If object is a vector perm_type must be entity")
		}
		res <- foreach(i=1:length(entity_sets), .options.snow=pbd$opts, .export=c("csea"),
			       .combine=rbind) %dopar% {
			csea(object, entity_sets[[i]], perm_type, statistic, power, iterations)
		}
	} else {
		stop("object must be of class seainput, vector, or lfcmodel")
	}

	close(pbd$pb)
	stopCluster(cl)
	gc()
	res
}

csea_parallel_matrix <- function(S, X, inds, lfc_column, statistic="ks", iterations=2000, permutation_matrix=NULL, cores=detectCores()) {
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  
  combn <- function(results, r) {
    results$obs_score <- c(results$obs_score, r$obs_score)
    results$p_value <- c(results$p_value, r$p_value)
    results
  }

  if(!is.null(permutation_matrix)) {
    iterations <- ncol(permutation_matrix) 
  }
  
  pbd <- build_txt_pb_opts(length(inds))
  res <- foreach(j=1:length(inds), .options.snow=pbd$opts, .export=c("csea_base", "calculate_W", "ols_solution", "draw_shuffled_col"), .combine=combn) %dopar% {
    # Calculate observed
    obs_lfc <- ols_solution(S, X, lfc_column)
    obs_score <- csea_base(obs_lfc, inds[[j]], statistic)

    gt <- 0
    for(i in 1:iterations) {
      if(is.null(permutation_matrix)) {
        shuffled_lfc <- draw_shuffled_col(S, X, lfc_column)
      } else {
        shuffled_lfc <- draw_perm_col(S, X, permutation_matrix[,i], lfc_column) 
      }
      shuffled_score <- csea_base(shuffled_lfc, inds[[j]], statistic)
      if(abs(shuffled_score)>=abs(obs_score)) {
        gt <- gt + 1
      }
    }
    list(obs_score=obs_score, p_value=gt/iterations)
  }
  close(pbd$pb)
  stopCluster(cl)
  gc()
  
  res
}

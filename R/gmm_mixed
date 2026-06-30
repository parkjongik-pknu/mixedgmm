library(mvtnorm)
library(mclust)


gmm_mixed <- function(X_Q = NULL, X_W = NULL, X_D = NULL, K = NULL, K_max = NULL, max_iter = 100, tol = 1e-6) {

  if (!is.null(K_max)) {
    if (K_max < 2) stop("K_max는 2 이상이어야 합니다.")
    k_seq <- 2:K_max
    is_search <- TRUE
  } else if (!is.null(K)) {
    if (K < 2) stop("K는 2 이상이어야 합니다.")
    k_seq <- K
    is_search <- FALSE
  } else {
    stop("K 또는 K_max 값을 입력해야 합니다. (예: K = 3 또는 K_max = 5)")
  }
  
  N <- max(nrow(X_Q), nrow(X_W), nrow(X_D))
  if (is.null(X_Q)) X_Q <- matrix(0, nrow = N, ncol = 0)
  if (is.null(X_D)) X_D <- matrix(0, nrow = N, ncol = 0)
  if (is.null(X_W)) X_W <- data.frame(matrix(NA, nrow = N, ncol = 0))
  
  P_Q <- ncol(X_Q)
  P_D <- ncol(X_D)
  P_W <- ncol(X_W)
  
  all_models <- list()
  models_summary <- data.frame(K = integer(), LogLik = numeric(), BIC = numeric(), ICL = numeric(), n_params = integer())
  best_bic <- Inf
  best_model <- NULL
  
  if (is_search) cat(sprintf("최적 군집 탐색 시작 (K = 2 ~ %d)...\n", K_max))
  

  for (current_K in k_seq) {
    if (is_search) cat(sprintf("-> K = %d 적합 중... ", current_K))
    
    fit_K <- tryCatch({
      X_num <- cbind(X_Q, X_D)
      Z <- matrix(0, nrow = N, ncol = current_K)
      
      if (ncol(X_num) > 0) {
        gmm_init <- suppressWarnings(mclust::Mclust(X_num, G = current_K, verbose = FALSE))
        if (!is.null(gmm_init) && !is.null(gmm_init$z)) {
          Z <- gmm_init$z
        } else {
          km <- kmeans(scale(X_num), centers = current_K, nstart = 10)
          for (k in 1:current_K) Z[, k] <- ifelse(km$cluster == k, 1, 0)
        }
      } else {
        rand_mat <- matrix(runif(N * current_K), N, current_K)
        Z <- rand_mat / rowSums(rand_mat)
      }
      
      mu_k <- list(); sigma_k <- list(); nu_k <- list(); alpha_k <- list()
      ll_history <- numeric(); curr_ll <- -Inf
      
      # EM 알고리즘
      for (iter in 1:max_iter) {
        # M-Step
        N_k <- colSums(Z)
        pi_k <- N_k / N
        
        for (k in 1:current_K) {
          if (P_Q > 0) { 
            mu_k[[k]] <- colSums(Z[, k] * X_Q) / max(N_k[k], 1e-10)
            diff_k <- sweep(X_Q, 2, mu_k[[k]])
            sigma_k[[k]] <- t(diff_k) %*% (Z[, k] * diff_k) / max(N_k[k], 1e-10) + diag(1e-5, P_Q) 
          }
          if (P_D > 0) {
            nu_k[[k]] <- colSums(Z[, k] * X_D) / max(N_k[k], 1e-10)
          }
          if (P_W > 0) { 
            alpha_list <- list()
            for (j in 1:P_W) { 
              lvls <- levels(X_W[, j])
              alpha_vec <- sapply(lvls, function(l) sum(Z[, k] * (X_W[, j] == l)) / max(N_k[k], 1e-10))
              alpha_list[[j]] <- alpha_vec / sum(alpha_vec) 
            }
            alpha_k[[k]] <- alpha_list 
          }
        }
        
        # E-Step
        f_k <- matrix(0, nrow = N, ncol = current_K)
        for (k in 1:current_K) {
          prob_X <- rep(1, N)
          if (P_Q > 0) prob_X <- prob_X * dmvnorm(X_Q, mean = mu_k[[k]], sigma = sigma_k[[k]])
          if (P_D > 0) { 
            for (j in 1:P_D) prob_X <- prob_X * dpois(X_D[, j], lambda = nu_k[[k]][j]) 
          }
          if (P_W > 0) { 
            for (j in 1:P_W) prob_X <- prob_X * alpha_k[[k]][[j]][as.character(X_W[, j])] 
          }
          f_k[, k] <- pi_k[k] * prob_X
        }
        
        f_total <- rowSums(f_k)
        f_total_safe <- pmax(f_total, 1e-300)
        Z <- sweep(f_k, 1, f_total_safe, "/")
        
        prev_ll <- curr_ll
        curr_ll <- sum(log(f_total_safe))
        ll_history <- c(ll_history, curr_ll)
        
        if (iter > 1 && abs(curr_ll - prev_ll) < tol) break
      }
      
      n_pi <- current_K - 1
      n_mu <- ifelse(P_Q > 0, current_K * P_Q, 0)
      n_sigma <- ifelse(P_Q > 0, current_K * (P_Q * (P_Q + 1)) / 2, 0)
      n_nu <- ifelse(P_D > 0, current_K * P_D, 0)
      n_alpha <- ifelse(P_W > 0, current_K * sum(sapply(1:P_W, function(j) length(levels(X_W[, j])) - 1)), 0)
      
      total_params <- n_pi + n_mu + n_sigma + n_nu + n_alpha
      BIC_val <- -2 * curr_ll + total_params * log(N)
      Z_safe <- pmax(Z, 1e-300)
      ICL_val <- BIC_val - 2 * sum(Z * log(Z_safe))
      
      list(
        K = current_K, pi = pi_k, mu = mu_k, sigma = sigma_k, nu = nu_k, alpha = alpha_k,
        Z = Z, loglik = ll_history, BIC = BIC_val, ICL = ICL_val, n_params = total_params, iterations = iter
      )
      
    }, error = function(e) {
      return(e$message)
    })

    if (is.list(fit_K)) {
      if (is_search) cat(sprintf("[성공] BIC: %.2f | ICL: %.2f\n", fit_K$BIC, fit_K$ICL))
      
      all_models[[as.character(current_K)]] <- fit_K
      models_summary <- rbind(models_summary, data.frame(
        K = fit_K$K, LogLik = tail(fit_K$loglik, 1), 
        BIC = fit_K$BIC, ICL = fit_K$ICL, n_params = fit_K$n_params
      ))
      
      if (fit_K$BIC < best_bic) {
        best_bic <- fit_K$BIC
        best_model <- fit_K
      }
    } else {
      if (is_search) cat(sprintf("[실패] %s\n", fit_K))
      models_summary <- rbind(models_summary, data.frame(
        K = current_K, LogLik = NA, BIC = NA, ICL = NA, n_params = NA
      ))
    }
  }
  

  if (is_search) {
    cat("\n탐색 완료!\n")
    return(list(
      best_model = best_model,
      models_summary = models_summary,
      all_models = all_models
    ))
  } else {
    if (is.null(best_model)) stop("모델 적합에 실패했습니다. (데이터 문제 혹은 수렴 실패)")
    # 단일 K 모드일 때는 리스트를 벗겨내고 결과 자체만 깔끔하게 반환
    return(best_model)
  }
}

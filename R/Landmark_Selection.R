#-------------- 5. Landmark Point Selection--------------#

choose_landmarks <- function(P,
                             center = NULL,
                             K = NULL,
                             dtheta = NULL,
                             angle_start = 0,
                             tol = 1.0,
                             max_tol = 5.0) {
  stopifnot(is.matrix(P) || is.data.frame(P))
  P <- as.matrix(P)
  stopifnot(ncol(P) == 2)
  colnames(P) <- c("x", "y")
  
  if (is.null(center)) center <- colMeans(P)
  stopifnot(length(center) == 2, is.finite(center[1]), is.finite(center[2]))
  
  if (is.null(K) && is.null(dtheta))
    stop("Provide either K or dtheta (degrees).")
  if (!is.null(dtheta)) {
    stopifnot(dtheta > 0)
    dtheta_rad <- dtheta * pi / 180
    K <- max(1L, round(2 * pi / dtheta_rad))
  }
  stopifnot(K >= 1)
  angles <- angle_start + 2 * pi * (0:(K - 1)) / K
  V <- sweep(P, 2, center, FUN = "-")     
  n <- nrow(P)
  
  pick_for_angle <- function(theta) {
    d <- c(cos(theta), sin(theta))
    nperp <- c(-sin(theta), cos(theta))
    t <- V %*% d          
    s <- V %*% nperp      
    tol_k <- tol
    chosen <- NA_integer_
    while (is.na(chosen) && tol_k <= max_tol) {
      cand <- which(abs(s) <= tol_k & t > 0)
      if (length(cand) > 0) {
        chosen <- cand[which.max(t[cand])]   
      } else {
        tol_k <- tol_k * 1.5
      }
    }
    if (is.na(chosen)) {
      phi <- atan2(V[, 2], V[, 1])
      adiff <- abs(atan2(sin(phi - theta), cos(phi - theta)))  
      chosen <- which.min(adiff)
    }
    chosen
  }
  
  idx <- vapply(angles, pick_for_angle, integer(1))
  L <- P[idx, , drop = FALSE]
  rownames(L) <- NULL
  attr(L, "center") <- center
  attr(L, "angles") <- angles
  attr(L, "indices") <- idx
  L
}

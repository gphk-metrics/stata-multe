cap mata mata drop multe_helper_ols()
cap mata mata drop multe_helper_olsw()

mata
void function MulTE(real colvector Y, real colvector X, real matrix W)
{
    n       = rows(Y)
    xlevels = uniqrows(X)
    X0      = (X :== xlevels[1])
    Xm      = designmatrix(X)
    alpha0  = multe_helper_ols(select(Y, X0), select(Wm, X0))

    psi_alpha0  = (X0 :* Wm :* (Y - Wm * alpha0)) * qrinv(cross(X0 :* Wm, Wm)/n)
    est = se_or = se_po = J(length(xlevels) - 1, 3, .)

    // common weights
    ps        = Wm * multe_helper_ols(Xm, Wm) // propensity scores
    Xt        = Xm - ps                       // residuals
    lam       = 1 :/ rowsum(1:/ps)
    alphac    = multe_helper_olsw(Y, Xm, lam :/ rowsum(ps :* Xm))
    resc      = Y - Xm * alphac
    est[., 3] = alphac[2::length(alphac)] :- alphac[1]

    // Fitted t values
    ts = Wm * multe_helper_ols(resc :* Xm :/ (ps:^2)  :* lam, Wm)
    s0 = Wm * multe_helper_ols(resc :* X0 :/ ps[., 1] :* (lam:^2) :/ (ps:^2), Wm)

    // loop (you are here)
}

real matrix function multe_helper_ols(real matrix Y, real matrix X)
{
    return(invsym(cross(X, X)) * cross(X, Y))
}

real matrix function multe_helper_olsw(real matrix Y, real matrix X, real matrix W)
{
    return(qrinv(cross(X :* W, X)) * cross(X :* W, Y))
}
end

// te_estimates <- function(Y, X, Wm) {
//     n <- length(Y)
//     X0 <- X==levels(X)[1]
//     Xm <- outer(X, levels(X)[-1], `==`)+0    ## X matrix
//     alpha0 <- lm(Y~0+Wm, subset=X0)$coefficients
//     psi_alpha0 <- (X0*Wm*(Y-drop(Wm %*% alpha0))) %*%
//         solve(crossprod(X0*Wm, Wm)/n)
//     est <- se_or <- se_po <- matrix(ncol=3, nrow=(length(levels(X))-1))
//
//     ## common weights
//     ps <- predict(lm(cbind(X0, Xm)~Wm)) # propensity scores
//     Xt <- cbind(X0, Xm)-ps              # residuals
//     lam <- 1/rowSums(1/ps)
//     rc <- lm(Y~0+X, weights=lam/rowSums(ps*cbind(X0, Xm)))
//     alphac <- rc$coefficients
//     est[, 3] <- alphac[-1]-alphac[1]
//     ## Fitted t values
//     ts <- predict(lm(rc$residuals * cbind(X0, Xm) / ps^2 * lam ~ 0+Wm))
//     s0 <- predict(lm(rc$residuals*X0/ps[, 1]*lam^2/ps^2 ~ 0+Wm))
//
//     for (k in seq_along(levels(X)[-1])) {
//         alphak <- lm(Y~0+Wm, subset=(X==levels(X)[k+1]))$coefficients
//         psi_alphak <- (Xm[, k]*(Y-drop(Wm %*% alphak))*Wm) %*%
//             solve(crossprod(Xm[, k]*Wm, Wm)/n)
//         ## ATE
//         est[k, 1] <- drop(colMeans(Wm) %*% (alphak-alpha0))
//         psi_or <- drop((psi_alphak-psi_alpha0) %*% colMeans(Wm))
//         se_or[k, 1] <- sqrt(var(psi_or)*(n-1)/n^2)
//         se_po[k, 1] <- sqrt(var(psi_or+drop(scale(Wm, scale = FALSE) %*%
//                                             (alphak-alpha0)))*(n-1)/n^2)
//         ## 1 at a time
//         s <- (X0 | X==levels(X)[k+1])
//         Xdot <- lm(Xm[, k]~Wm, subset=s)$residuals
//         rk <- lm(Y[s]~0+Xdot+Wm[s, ])
//         est[k, 2] <- rk$coefficients[1]
//         se_po[k, 2] <- sqrt(sum(rk$residuals^2*Xdot^2)/sum(Xdot^2)^2)
//         eps <- Xm[s, k]*(Y[s]-Wm[s, ] %*% alphak) +
//             X0[s]*(Y[s]-Wm[s, ] %*% alpha0)
//         se_or[k, 2] <- sqrt(sum(eps^2*Xdot^2)/sum(Xdot^2)^2)
//         ## common weights
//         psi_or <- lam*(Xm[, k]*(Y-Wm %*% alphak)/ps[, k+1]-
//                        X0*(Y-Wm %*% alpha0)/ps[, 1])/mean(lam)
//         sk <- predict(lm(rc$residuals*Xm[, k]/ps[, k+1]*lam^2/ps^2 ~ 0+Wm))
//         psi_po <- (lam*rc$residuals*(Xm[, k]/ps[, k+1]-X0/ps[, 1]) +
//             Xt[, 1]*ts[, 1]-Xt[, k+1]*ts[, k+1]+rowSums(Xt * (sk-s0)))/mean(lam)
//         se_po[k, 3] <- sqrt(var(psi_po)*(n-1)/n^2)
//         se_or[k, 3] <- sqrt(var(drop(psi_or))*(n-1)/n^2)
//
//     }
//     rownames(est) <- levels(X)[-1]
//     rownames(se_po) <- rep("se", nrow(est))
//     rownames(se_or) <- rep("oracle_se", nrow(est))
//     ret <- rbind(est, se_po, se_or)
//     colnames(ret) <- c("ATE", "1-at-a-time", "common weights")
//     rbind(ret[seq_len(nrow(ret)) %% 2 == 1, ],
//           ret[seq_len(nrow(ret)) %% 2 == 0, ])
// }

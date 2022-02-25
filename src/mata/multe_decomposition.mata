cap mata mata drop multe_decomposition()
cap mata mata drop multe_decomposition2()
cap mata mata drop multe_helper_antiselect()

mata
real vector function multe_helper_antiselect(real vector x, real vector ix)
{
    real vector sel
    sel = J(length(x), 1, 1)
    sel[ix] = J(length(ix), 1, 0)
    return(x[selectindex(sel)])
}

// From R> matrix of controls Wm including intercept; X must be a
// From R> factor, first level to be dropped
void function multe_decomposition(string scalar Yvar, string scalar Tvar, real matrix Wm, string scalar touse)
{
Xm = Xm[., 2::k]

    tXm      = Xm - Wm * multe_helper_ols(Xm, Wm) // Xtilde
    rlcoef   = multe_helper_ols(Y, (Xm, Wm))
    rlres    = Y - (Xm, Wm) * rlcoef
    psi_beta = (rlres :* tXm) * invsym(cross(tXm, tXm))
    beta     = rlcoef[1::(k-1)]

    WmXm = J(n, kw + (k - 1) * kw, 0)
    WmXm[|1, 1 \ n, kw|] = Wm
    for(j = 1; j < k; j++) {
        j1 = kw + (j - 1) * kw + 1
        j2 = kw + (j - 0) * kw
        WmXm[|1, j1 \ n, j2|] = Xm[., j] :* Wm
    }
    ricoef   = multe_helper_ols(Y, WmXm) // interactions
    rires    = Y - WmXm * ricoef
    avar_own = avar_cb = own = beta

    for(j = 1; j < k; j++) {
        selk   = (j * kw + 1)::(j * kw)
        gammak = ricoef[selk]
        selk   = multe_helper_antiselect(1::cols(WmXm), selk)
        Xg     = multe_helper_olsr(Xm[., j] :* Wm, WmXm[., selk]).residuals
        rk     = multe_helper_olsr(Xm[., j] :* Wm, (Xm, Wm))
        deltak = rk.coefficients[j, .]
        selk   = multe_helper_antiselect(1::cols(tXm), j)
        dtX    = multe_helper_olsr(tXm[., j], tXm[., selk]).residuals // doubletilde(X)

        psi_ownk = (rires :* Xg) * qrsolve(cross(Xg, Xg), deltak) +
            (dtX :* rk.residuals :/ sum(dtX:^2)) * gammak

        own[j]      = sum(gammak :* deltak)
        avar_own[j] = variance(psi_ownk)
        avar_cb[j]  = variance(psi_beta[., j] - psi_ownk)
    }

    // est <- cbind(beta, own, "cont. bias"=beta-own)
    // rownames(est) <- levels(X)[-1]
    // se <- sqrt((n-1)*cbind(diag(var(psi_beta)), avar_own, avar_cb))
    // rownames(se) <- rep("se", nrow(se))
    // rbind(est, se)[rep(seq_len(nrow(est)), each=2) + c(0, nrow(est)), ]

}
end

// Also computes worst-case contamination bias
// decomposition2 <- function(Y, X, Wm) {
//     Xm <- outer(X, levels(X)[-1], `==`) + 0  # X matrix
//     ri <- lm(Y~0+Wm+Xm:Wm)                   # interactions
//     K <- ncol(Xm)
//     L <- ncol(Wm)
// 
//     gamma <- ri$coefficients[-seq_len(ncol(Wm))]
//     psi_gamma <- ((ri$residuals*model.matrix(ri)) %*%
//                   solve(crossprod(model.matrix(ri))))[, -seq_len(ncol(Wm))]
//     rd <- lm(model.matrix(lm(Y~0+Wm:Xm))~0+Xm+Wm)                   # delta
//     ## Sort columns by size
//     sortc <- function(v, decreasing=FALSE) {
//         as.vector(apply(matrix(v, nrow=L), 2,
//                         function(x) sort(x, index.return=TRUE,
//                                          decreasing=decreasing)$ix)) +
//             rep((0:(K-1))*L, each=L)
//     }
//     gi <- sortc(gamma)
//     gd <- sortc(gamma, TRUE)
//     est <- se <- matrix(nrow=K, ncol=5)
//     ## Standard errors
//     for (k in seq_len(K)) {
//         ddX <- lm(Xm[, k]~0+Xm[, -k]+Wm)$residuals # ddot(X)
//         deltak <- rd$coefficients[k, ]
//         psi_deltak <- ddX*rd$residuals / sum(ddX^2)
// 
//         M <- kronecker(diag(K), rep(1, L))
//         psi <- psi_deltak %*% (gamma*M) + psi_gamma %*% (deltak*M)
//         estk <- drop(gamma %*% (deltak*M))
//         ## Sort delta columns
//         di <- sortc(deltak)
// 
//         psimax <- psi_deltak[, di] %*% (gamma[gi]*M) +
//             psi_gamma[, gi] %*% (deltak[di]*M)
//         psimin <- psi_deltak[, di] %*% (gamma[gd]*M) +
//             psi_gamma[, gd] %*% (deltak[di]*M)
//         sqrt(colSums(psimax^2))
//         sqrt(colSums(psimin^2))
//         est[k, ] <- c(sum(estk), estk[k], sum(estk[-k]),
//                       sum(drop(gamma[gi] %*% (deltak[di]*M))[-k]),
//                       sum(drop(gamma[gd] %*% (deltak[di]*M))[-k]))
//         se[k, ] <- sqrt(c(sum(rowSums(psi)^2), sum(psi[, k]^2),
//                           sum(rowSums(psi[, -k, drop=FALSE])^2),
//                           sum(rowSums(psimax[, -k, drop=FALSE])^2),
//                           sum(rowSums(psimin[, -k, drop=FALSE])^2)))
//     }
//     rownames(se) <- rep("se", K)
//     rownames(est) <- levels(X)[-1]
//     colnames(est) <- c("beta", "own", "cont. bias", "maxbias", "minbias")
//     rbind(est, se)[rep(seq_len(nrow(est)), each=2) + c(0, nrow(est)), ]
// }

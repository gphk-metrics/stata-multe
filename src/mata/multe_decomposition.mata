cap mata mata drop multe_decomposition()

mata
// From R> matrix of controls Wm including intercept; X must be a
// From R> factor, first level to be dropped
// void function multe_decomposition(string scalar Yvar, string scalar Tvar, real matrix Wm, string scalar touse)
// {
//     Xm = Xm[., 2::k]
// 
//     tXm      = Xm - Wm * multe_helper_ols(Xm, Wm) // Xtilde
//     rlcoef   = multe_helper_ols(Y, (Xm, Wm))
//     rlres    = Y - (Xm, Wm) * rlcoef
//     psi_beta = (rlres :* tXm) * invsym(cross(tXm, tXm))
//     beta     = rlcoef[1::(k-1)]
// 
//     WmXm = J(n, kw + (k - 1) * kw, 0)
//     WmXm[|1, 1 \ n, kw|] = Wm
//     for(j = 1; j < k; j++) {
//         j1 = kw + (j - 1) * kw + 1
//         j2 = kw + (j - 0) * kw
//         WmXm[|1, j1 \ n, j2|] = Xm[., j] :* Wm
//     }
//     ri = multe_helper_olsr(Y, WmXm) // interactions
//     avar_own = avar_cb = own = beta
// 
//     for(j = 1; j < k; j++) {
//         selk   = (j * kw + 1)::((j + 1) * kw)
//         gammak = ri.coefficients[selk]
//         selk   = multe_helper_antiselect(1::cols(WmXm), selk)
//         rXg    = multe_helper_olsr(Xm[., j] :* Wm, WmXm[., selk])
//         Xg     = rXg.residuals
//         rk     = multe_helper_olsr(Xm[., j] :* Wm, (Xm, Wm))
//         deltak = rk.coefficients[j, .]'
//         selk   = multe_helper_antiselect(1::cols(tXm), j)
//         rdtX   = multe_helper_olsr(tXm[., j], tXm[., selk])
//         dtX    = rdtX.residuals // doubletilde(X)
// 
//         psi_ownk = (ri.residuals :* Xg) * qrsolve(cross(Xg, Xg), deltak) +
//             (dtX :* rk.residuals :/ sum(dtX:^2)) * gammak
// 
//         own[j]      = sum(gammak :* deltak)
//         avar_own[j] = variance(psi_ownk)
//         avar_cb[j]  = variance(psi_beta[., j] - psi_ownk)
//     }
// 
// // TODO: xx "beta", "own", "cont. bias"
// // TODO: xx rownames are labels or "se"
//     est = beta, own, beta - own
//     se  = sqrt((n - 1) :* (diagonal(variance(psi_beta)), avar_own, avar_cb))
//     rowshape((est, se), 2 * rows(est))
// }

// TODO: xx some objects here are defined in the function above
// TODO: xx incorporate this into MulTE()

// void function multe_decomposition(string scalar Yvar, string scalar Tvar, real matrix Wm, string scalar touse)
// {
//     gamma = ri.coefficients[|(kw+1) \ cols(WmXm)|]
//     psi_gamma = ((ri.residuals :* WmXm) *
//                   invsym(cross(WmXm, WmXm)))[|(1, kw+1) \ (n, cols(WmXm))|]
//     rd = multe_helper_olsr(WmXm[|(1, kw+1) \ (n, cols(WmXm))|], (Xm, Wm)) // delta
// 
//     // Sort columns by size
//     ghelper = rowshape(colshape(J(kw, 1, 1::(k-1)), k-1)', kw * (k - 1))
//     gi  = order((ghelper, gamma), (1, 2))
//     gd  = order((ghelper, gamma), (1, -2))
//     est = se = J(k-1, 5, .)
// 
//     // Standard errors
//     for (j = 1; j < k; j++) {
//         rddX = multe_helper_olsr(Xm[., j], (multe_helper_antiselect(Xm, j), Wm))
//         ddX  = rddX.residuals // ddot(X)
// 
//         deltak     = rd.coefficients[j, .]'
//         psi_deltak = ddX :* rd.residuals / sum(ddX:^2)
// 
//         M    = I(k - 1)#J(kw, 1, 1)
//         psi  = psi_deltak * (gamma :* M) + psi_gamma * (deltak :* M)
//         estk = gamma' * (deltak :* M)
// 
//         // Sort delta columns
//         di = order((ghelper, deltak), (1, 2))
// 
//         psimax = psi_deltak[., di] * (gamma[gi] :* M) + psi_gamma[., gi] * (deltak[di] :* M)
//         psimin = psi_deltak[., di] * (gamma[gd] :* M) + psi_gamma[., gd] * (deltak[di] :* M)
// 
//         est[j,.] = (
//             sum(estk),
//             estk[j],
//             sum(multe_helper_antiselect(estk, j)),
//             sum(gamma[gi]' * multe_helper_antiselect(deltak[di] :* M, j)),
//             sum(gamma[gd]' * multe_helper_antiselect(deltak[di] :* M, j))
//         )
// 
//         se[j,.] = sqrt((
//             sum(rowsum(psi):^2),
//             sum(psi[., j]:^2),
//             sum(rowsum(multe_helper_antiselect(psi, j)):^2),
//             sum(rowsum(multe_helper_antiselect(psimax, j)):^2),
//             sum(rowsum(multe_helper_antiselect(psimin, j)):^2)
//         ))
//     }
// 
// // TODO: xx "beta", "own", "cont. bias", "maxbias", minbias"
// // TODO: xx rownames are labels or "se"
//     rowshape((est, se), 2 * rows(est))
// }
end

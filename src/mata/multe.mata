cap mata mata drop MulTE()
cap mata mata drop MulTE_Estimates()
cap mata mata drop MulTE_Decomposition()
cap mata mata drop MulTE_Results()

mata
// Yvar  = "`depvar'"
// Tvar  = "`treatment'"
// touse = "`touse'"

struct MulTE_Results
{
    class MulTE_Estimates scalar estimates
    class MulTE_Decomposition scalar decomposition
}

class MulTE_Estimates
{
    real matrix est
    real matrix se_po
    real matrix se_or
    real   vector Tvalues
    string vector Tlabels
    string vector colnames

    void new()
    void print()
    void save()
}

class MulTE_Decomposition
{
    real matrix tmp
    real matrix est
    real matrix se
    real matrix tauhat
    real matrix lambda
    real   vector Tvalues
    string vector Tlabels
    string vector colnames

    void new()
    void print()
    void save()
}

struct MulTE_Results scalar MulTE(string scalar Yvar, string scalar Tvar, real matrix Wm, string scalar touse)
{
    struct MulTE_Results scalar results
    struct multe_helper_results scalar rk, ri, rd, rddX

    string scalar Tlab

    real scalar j, n, k, kw, j1, j2
    real vector Y, X
    real vector xlevels, X0, alpha0, lam, rcalpha, rcres, ts, s0, Wmean
    real matrix Xm, psi_alpha0, ps, Xt, WmXm
    real matrix gamma, psi_gamma, deltak, psi_deltak, ddX, M, psi
    real matrix est, estk, se_or, se_po, se, psimin, psimax
    real matrix alphak, psi_alphak, psi_or, psi_po
    real vector s, Xdot, eps, sk, ghelper, gi, gd, di

    // -----------------------------------------------------------------
    // Setup
    // -----------------------------------------------------------------

    Y       = st_data(., Yvar, touse)
    X       = st_data(., Tvar, touse)
    n       = rows(Y)
    xlevels = uniqrows(X)
    k       = length(xlevels)
    kw      = cols(Wm)
    X0      = (X :== xlevels[1])
    Xm      = designmatrix(X)
    alpha0  = multe_helper_ols(select(Y, X0), select(Wm, X0))

    psi_alpha0  = (X0 :* Wm :* (Y - Wm * alpha0)) * qrinv(cross(X0 :* Wm, Wm)/n)
    est = se_or = se_po = J(k - 1, 3, .)
// TODO: xx why isn't this 2 by k instead of k-2 by 3? think abt it

    // common weights
    ps        = Wm * multe_helper_ols(Xm, Wm) // propensity scores
    Xt        = Xm - ps                       // residuals
    lam       = 1 :/ rowsum(1:/ps)
    rcalpha   = multe_helper_olsw(Y, Xm, lam :/ rowsum(ps :* Xm))
    rcres     = Y - Xm * rcalpha
    est[., 3] = rcalpha[2::length(rcalpha)] :- rcalpha[1]

    // Fitted t values
    ts = Wm * multe_helper_ols(rcres :* Xm :/ (ps:^2)  :* lam, Wm)
    s0 = Wm * multe_helper_ols(rcres :* X0 :/ ps[., 1] :* (lam:^2) :/ (ps:^2), Wm)

// TODO: xx figure out if selectindex(s) is faster
// TODO: xx or maybe make subset vars at start of loop

    // -----------------------------------------------------------------
    // TE estimates
    // -----------------------------------------------------------------

    Wmean = colsum(Wm) / n
    for(j = 2; j <= k; j++) {
        alphak  = multe_helper_ols(select(Y, Xm[., j]), select(Wm, Xm[., j]))
        psi_alphak = (Xm[.,j] :* (Y - Wm * alphak) :* Wm) * qrinv(cross(Xm[., j] :* Wm, Wm) / n)

        // ATE
        est[j - 1, 1] = (Wmean * (alphak - alpha0))
        psi_or = ((psi_alphak - psi_alpha0) * Wmean')
        se_or[j - 1, 1] = sqrt(variance(psi_or) * (n - 1) / n^2)
        se_po[j - 1, 1] = sqrt(variance(psi_or + ((Wm :-  Wmean) * (alphak - alpha0))) * (n - 1) / n^2)

        // One treatment at a time

        s               = (X0 :| Xm[., j])
        Xdot            = select(Xm[., j], s) - select(Wm, s) * multe_helper_ols(select(Xm[., j], s), select(Wm, s))
        rk              = multe_helper_olsr(select(Y, s), (Xdot, select(Wm, s)))
        est[j - 1, 2]   = rk.coefficients[1]
        se_po[j - 1, 2] = sqrt(sum((rk.residuals:^2) :* (Xdot:^2)) / sum(Xdot:^2):^2)

        eps = select(Xm[., j], s) :* (select(Y, s) - select(Wm, s) * alphak) +
              select(X0, s) :* (select(Y, s) - select(Wm, s) * alpha0)
        se_or[j - 1, 2] = sqrt(sum((eps:^2) :* (Xdot:^2))/sum(Xdot:^2)^2)

        // common weights
        psi_or = lam :* (Xm[., j] :* (Y - Wm * alphak) :/ ps[., j] -
                 X0 :* (Y - Wm * alpha0) :/ ps[., 1]) :/ mean(lam)
        sk = Wm * multe_helper_ols(rcres :* Xm[., j] :/ ps[., j] :* (lam:^2) :/ (ps:^2), Wm)
        psi_po = (lam :* rcres :* (Xm[., j] :/ ps[., j] - X0 :/ ps[., 1]) +
                  Xt[., 1] :* ts[., 1] - Xt[., j] :* ts[., j] + rowsum(Xt :* (sk - s0))) / mean(lam)
        se_po[j-1, 3] = sqrt(variance(psi_po)*(n-1)/n^2)
        se_or[j-1, 3] = sqrt(variance(psi_or)*(n-1)/n^2)
    }

    Tlab = st_varvaluelabel(Tvar)
    results.estimates.est   = est
    results.estimates.se_po = se_po
    results.estimates.se_or = se_or
    results.estimates.Tvalues = xlevels
    results.estimates.Tlabels = Tlab == ""? strofreal(xlevels): st_vlmap(Tlab, xlevels)

// TODO: xx decide whether to make decomposition a separate function

    // -----------------------------------------------------------------
    // Decomposition (not run)
    // -----------------------------------------------------------------

    // From R> matrix of controls Wm including intercept; X must be a
    // From R> factor, first level to be dropped
    Xm   = Xm[., 2::k]
    WmXm = J(n, kw + (k - 1) * kw, 0)
    WmXm[|1, 1 \ n, kw|] = Wm
    for(j = 1; j < k; j++) {
        j1 = kw + (j - 1) * kw + 1
        j2 = kw + (j - 0) * kw
        WmXm[|1, j1 \ n, j2|] = Xm[., j] :* Wm
    }
    ri = multe_helper_olsr(Y, WmXm) // interactions

    gamma = ri.coefficients[|(kw+1) \ cols(WmXm)|]
    psi_gamma = ((ri.residuals :* WmXm) *
                  invsym(cross(WmXm, WmXm)))[|(1, kw+1) \ (n, cols(WmXm))|]
    rd = multe_helper_olsr(WmXm[|(1, kw+1) \ (n, cols(WmXm))|], (Xm, Wm)) // delta

    // Sort columns by size
    ghelper = rowshape(colshape(J(kw, 1, 1::(k-1)), k-1)', kw * (k - 1))
    gi  = order((ghelper, gamma), (1, 2))
    gd  = order((ghelper, gamma), (1, -2))
    est = se = J(k-1, 5, .)

    // Standard errors
    for (j = 1; j < k; j++) {
        rddX = multe_helper_olsr(Xm[., j], (multe_helper_antiselect(Xm, j), Wm))
        ddX  = rddX.residuals // ddot(X)

        deltak     = rd.coefficients[j, .]'
        psi_deltak = ddX :* rd.residuals / sum(ddX:^2)

        M    = I(k - 1)#J(kw, 1, 1)
        psi  = psi_deltak * (gamma :* M) + psi_gamma * (deltak :* M)
        estk = gamma' * (deltak :* M)

        // Sort delta columns
        di = order((ghelper, deltak), (1, 2))

        psimax = psi_deltak[., di] * (gamma[gi] :* M) + psi_gamma[., gi] * (deltak[di] :* M)
        psimin = psi_deltak[., di] * (gamma[gd] :* M) + psi_gamma[., gd] * (deltak[di] :* M)

        est[j,.] = (
            sum(estk),
            estk[j],
            sum(multe_helper_antiselect(estk, j)),
            sum(gamma[gi]' * multe_helper_antiselect(deltak[di] :* M, j)),
            sum(gamma[gd]' * multe_helper_antiselect(deltak[di] :* M, j))
        )

        se[j,.] = sqrt((
            sum(rowsum(psi):^2),
            sum(psi[., j]:^2),
            sum(rowsum(multe_helper_antiselect(psi, j)):^2),
            sum(rowsum(multe_helper_antiselect(psimax, j)):^2),
            sum(rowsum(multe_helper_antiselect(psimin, j)):^2)
        ))
 
    }

    // Control-specific TEs and weights
    1
    gammam = rowshape(gamma, k-1) // w x k matrix
    2
    tauhat = Wm * gammam' // N x k matrix
    // confirmed this is good. N x 2 matrix, with gammas populated for appropriate groups for the whole sample (N)
    3
    // TODO: add lambdas. I'm not sure these lambdas are the same as the ones in my Stata code.
    lambda = Xt[.,2::k] :* Xm * invsym(Xt[.,2::k]' * Xt[.,2::k])
    4
    beta_hat = colsum(tauhat :* lambda)
    5
    beta_hat

// TODO: xx "beta", "own", "cont. bias", "maxbias", minbias"
// TODO: xx rownames are labels or "se"
    results.decomposition.est = est
    results.decomposition.se  = se
    results.decomposition.tmp = rowshape((est, se), 2 * rows(est))
    results.decomposition.Tvalues = xlevels
    results.decomposition.Tlabels = results.estimates.Tlabels
    6
    results.decomposition.tauhat = tauhat
    7
    results.decomposition.lambda = lambda

    results.decomposition.tmp
    return(results)
}

void function MulTE_Estimates::new()
{
    colnames = ("ATE", "1-at-a-time", "common weights")'
}

void function MulTE_Estimates::save(string scalar outmatrix)
{
    string vector rownames
    rownames = Tlabels[2::length(Tvalues)], J(rows(est), 1, "se"), J(rows(est), 1, "oracle_se")
    rownames = rowshape(rownames, rows(est) * 3)

    st_matrix(outmatrix, rowshape((est, se_po, se_or), rows(est) * 3))
    st_matrixcolstripe(outmatrix, (J(cols(est), 1, ""), colnames))
    st_matrixrowstripe(outmatrix, (J(3 * rows(est), 1, ""), rownames))
}

void function MulTE_Estimates::print(| real scalar digits)
{
    real scalar i, j
    real vector lengths
    string matrix fmt_res, fmt_est, fmt_se_po, fmt_se_or
    string vector rownames, formats

    rownames = Tlabels[2::length(Tvalues)], J(rows(est), 1, "se"), J(rows(est), 1, "oracle_se")
    rownames = rowshape(rownames, rows(est) * 3)

    if ( args() < 1 ) digits = 6

    fmt_est   = J(rows(est), cols(est), "")
    fmt_se_po = J(rows(est), cols(est), "")
    fmt_se_or = J(rows(est), cols(est), "")

    for (i = 1; i <= rows(est); i++) {
        for (j = 1; j <= cols(est); j++) {
            fmt_est[i, j]   = strtrim(sprintf("%21." + strofreal(digits) + "f", est[i, j]))
            fmt_se_po[i, j] = "(" + strtrim(sprintf("%21." + strofreal(digits) + "f", se_po[i, j])) + ")"
            fmt_se_or[i, j] = "(" + strtrim(sprintf("%21." + strofreal(digits) + "f", se_or[i, j])) + ")"
        }
    }

    fmt_res = rowshape((fmt_est, fmt_se_po, fmt_se_or), rows(est) * 3)
    lengths = max(strlen(rownames)), colmax(strlen(colnames' \ fmt_res))
    formats = " %" :+ strofreal(lengths) :+ "s "

    printf("|")
    printf(formats[1], "")
    for (j = 1; j <= length(colnames); j++) {
        printf("|")
        printf(formats[j + 1], colnames[j])
    }
    printf("|\n")

    printf("|")
    for (j = 1; j <= length(lengths); j++) {
        printf(formats[j], "-" * lengths[j])
        printf("|")
    }
    printf("\n")

    for (i = 1; i <= rows(fmt_res); i++) {
        printf("|")
        printf(formats[1], rownames[i])
        for (j = 1; j <= cols(fmt_res); j++) {
            printf("|")
            printf(formats[j + 1], fmt_res[i, j])
        }
        printf("|\n")
    }
}


void function MulTE_Decomposition::new()
{
    colnames = ("beta", "own", "cont bias", "maxbias", "minbias")'
}

void function MulTE_Decomposition::save(string scalar outmatrix)
{
    string vector rownames
    rownames = Tlabels[2::length(Tvalues)], J(rows(est), 1, "se")
    rownames = rowshape(rownames, rows(est) * 2)

    st_matrix(outmatrix, rowshape((est, se), rows(est) * 2))
    st_matrixcolstripe(outmatrix, (J(cols(est), 1, ""), colnames))
    st_matrixrowstripe(outmatrix, (J(2 * rows(est), 1, ""), rownames))
}

void function MulTE_Decomposition::print(| real scalar digits)
{
}
end

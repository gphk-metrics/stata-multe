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
    real scalar n
    real scalar k
    real matrix est
    real matrix se_po
    real matrix se_or
    real   vector Tvalues
    string vector Tlabels
    string vector colnames
    string scalar Tvar
    string scalar Yvar

    real matrix po_vcov
    real matrix or_vcov

    void new()
    void print()
    void save()
    void post()
}

class MulTE_Decomposition
{
    real scalar n
    real scalar k
    real matrix est
    real matrix se
    real matrix tmp
    real matrix gammam
    real matrix delta_pr

    real   vector Tvalues
    string vector Tlabels
    string vector Tvar
    string vector Yvar
    string vector tauhat_names
    string vector lambda_names
    string vector colnames

    real matrix tauhat()
    real matrix lambda()
    void new()
    void print()
    void save()
}

struct MulTE_Results scalar MulTE(string scalar Yvar, string scalar Tvar, real matrix Wm, string scalar touse)
{
    struct MulTE_Results scalar results
    struct multe_helper_results scalar rk, ri, rd, rddX

    string vector lambda_names, tauhat_names
    string scalar Tlab

    real scalar i, l, j, n, k, kw, j1, j2
    real vector Y, X
    real vector xlevels, X0, alpha0, lam, rcalpha, rcres, ts, s0, Wmean
    real matrix Xm, psi_alpha0, ps, Xt, WmXm
    real matrix gamma, psi_gamma, deltak, psi_deltak, ddX, M, psi
    real matrix est, estk, se_or, se_po, se, psimin, psimax
    real matrix alphak, psi_alphak, psi_or, psi_po
    real matrix delta_kl, delta_pr, gammam
    real matrix psi_pom, psi_orm, var_po_onem, var_or_onem, psi_pom_tl, po_vcov, psi_orm_tl, or_vcov
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
    Xm      = J(n, k, 0)
    for (j = 1; j <= k; j++) {
        Xm[., j] = (X :== xlevels[j])
    }
    X0      = Xm[., 1]
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

    var_po_onem = var_or_onem = J(k, 1, 0)    
    psi_pom = psi_orm = J(n, k*3, 0)
    Wmean = mean(Wm)
    for(j = 2; j <= k; j++) {
        alphak  = multe_helper_ols(select(Y, Xm[., j]), select(Wm, Xm[., j]))
        psi_alphak = (Xm[.,j] :* (Y - Wm * alphak) :* Wm) * qrinv(cross(Xm[., j] :* Wm, Wm) / n)

        // ATE
        est[j - 1, 1] = (Wmean * (alphak - alpha0))
        psi_or = ((psi_alphak - psi_alpha0) * Wmean')
        se_or[j - 1, 1] = sqrt(variance(psi_or) * (n - 1) / n^2)
        se_po[j - 1, 1] = sqrt(variance(psi_or + ((Wm :-  Wmean) * (alphak - alpha0))) * (n - 1) / n^2)
        psi_orm[.,j] = psi_or * sqrt((n - 1) / n^2)
        psi_pom[.,j] = (psi_or + ((Wm :-  Wmean) * (alphak - alpha0))) * sqrt((n - 1) / n^2)

        // One treatment at a time
        s               = (X0 :| Xm[., j])
        Xdot            = select(Xm[., j], s) - select(Wm, s) * multe_helper_ols(select(Xm[., j], s), select(Wm, s))
        rk              = multe_helper_olsr(select(Y, s), (Xdot, select(Wm, s)))
        est[j - 1, 2]   = rk.coefficients[1]
        se_po[j - 1, 2] = sqrt(sum((rk.residuals:^2) :* (Xdot:^2)) / sum(Xdot:^2):^2)

        eps = select(Xm[., j], s) :* (select(Y, s) - select(Wm, s) * alphak) +
              select(X0, s) :* (select(Y, s) - select(Wm, s) * alpha0)
        se_or[j - 1, 2] = sqrt(sum((eps:^2) :* (Xdot:^2))/sum(Xdot:^2)^2)
        var_po_onem[j, 1] = sum((rk.residuals:^2) :* (Xdot:^2)) / sum(Xdot:^2):^2
        var_or_onem[j, 1] = sum((eps:^2) :* (Xdot:^2)) / sum(Xdot:^2)^2

        // common weights
        psi_or = lam :* (Xm[., j] :* (Y - Wm * alphak) :/ ps[., j] -
                 X0 :* (Y - Wm * alpha0) :/ ps[., 1]) :/ mean(lam)
        sk = Wm * multe_helper_ols(rcres :* Xm[., j] :/ ps[., j] :* (lam:^2) :/ (ps:^2), Wm)
        psi_po = (lam :* rcres :* (Xm[., j] :/ ps[., j] - X0 :/ ps[., 1]) +
                  Xt[., 1] :* ts[., 1] - Xt[., j] :* ts[., j] + rowsum(Xt :* (sk - s0))) / mean(lam)
        se_po[j-1, 3] = sqrt(variance(psi_po)*(n-1)/n^2)
        se_or[j-1, 3] = sqrt(variance(psi_or)*(n-1)/n^2)
        psi_pom[., j+(k*2)] = psi_po*sqrt((n-1)/n^2)
        psi_orm[., j+(k*2)] = psi_or*sqrt((n-1)/n^2)

    }

    // Compute vcov matrices
    // TODO: confirm that the variance formula is correct
    psi_pom_tl = psi_pom :- (colsum(psi_pom) :/ n)
    po_vcov = (psi_pom_tl' * psi_pom_tl) / n
    psi_orm_tl = psi_orm :- (colsum(psi_orm) :/ n)
    or_vcov = (psi_orm_tl' * psi_orm_tl) / n
    
    for(j=2; j<=k; j++) {
        po_vcov[j+k, j+k] = var_po_onem[j,1]
        or_vcov[j+k, j+k] = var_or_onem[j,1]
    }

    po_vcov = blockdiag(po_vcov[|1,1 \ k, k|], blockdiag(po_vcov[|k+1, k+1 \ 2*k, 2*k|], po_vcov[|2*k+1, 2*k+1 \ 3*k, 3*k|]))
    or_vcov = blockdiag(or_vcov[|1,1 \ k, k|], blockdiag(or_vcov[|k+1, k+1 \ 2*k, 2*k|], or_vcov[|2*k+1, 2*k+1 \ 3*k, 3*k|]))

    Tlab = st_varvaluelabel(Tvar)
    results.estimates.n     = n
    results.estimates.k     = k
    results.estimates.est   = est
    results.estimates.se_po = se_po
    results.estimates.se_or = se_or
    results.estimates.Tvalues = xlevels
    results.estimates.Tlabels = Tlab == ""? strofreal(xlevels): st_vlmap(Tlab, xlevels)
    results.estimates.Tvar    = Tvar
    results.estimates.Yvar    = Yvar

    // TODO: not sure why when I remove these two lines the mata function has a conformability error. I'm pretty sure these objects aren't used.
    results.estimates.po_vcov = po_vcov
    results.estimates.or_vcov = or_vcov

// TODO: xx decide whether to make decomposition a separate function

    // -----------------------------------------------------------------
    // Decomposition (not run)
    // -----------------------------------------------------------------

    // From R> matrix of controls Wm including intercept; X must be a
    // From R> factor, first level to be dropped

    // NB> Wm does not have a constant because it contains the full
    // NB> matrix of dummies.

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
            sum(gamma[gd]' * multe_helper_antiselect(deltak[di] :* M, j)),
            sum(gamma[gi]' * multe_helper_antiselect(deltak[di] :* M, j))
        )

        se[j,.] = sqrt((
            sum(rowsum(psi):^2),
            sum(psi[., j]:^2),
            sum(rowsum(multe_helper_antiselect(psi, j)):^2),
            sum(rowsum(multe_helper_antiselect(psimin, j)):^2),
            sum(rowsum(multe_helper_antiselect(psimax, j)):^2)
        ))

    }

    // Control-specific TEs and weights
    //     - rd has coefficients of reg of X * W on X, W
    //     - 1..(k-1) rows are coefs of X[., 2..k]
    //     - reshape so lth set of 1..(k-1) columns has coefs
    //       for X[., l+1] corresponding to X[., 2..k] * W
    //       (i.e. X[., l+1] for X[., 2] * W, ..., X[., k] * W)

    i = 0
    tauhat_names = J(1, k-1, "")
    lambda_names = J(1, (k-1)*(k-1), "")
    for (l = 1; l < k; l++) {
        tauhat_names[l] = sprintf("%g", l)
        for (j = 1; j < k; j++) {
            lambda_names[++i] = sprintf("%g%g", l, j)
        }
    }

    delta_kl = rowshape(rowshape(rd.coefficients[1..(k-1),.], 1), (k-1)*(k-1))'
    delta_pr = delta_kl :/ (mean(Wm)')
    gammam   = rowshape(gamma, k-1)

    results.decomposition.n   = n
    results.decomposition.k   = k
    results.decomposition.est = est
    results.decomposition.se  = se
    results.decomposition.tmp = rowshape((est, se), 2 * rows(est))
    results.decomposition.Tvalues  = xlevels
    results.decomposition.Tlabels  = results.estimates.Tlabels
    results.decomposition.Tvar     = Tvar
    results.decomposition.Yvar     = Yvar
    results.decomposition.delta_pr = delta_pr
    results.decomposition.gammam   = gammam  
    results.decomposition.tauhat_names = tauhat_names
    results.decomposition.lambda_names = lambda_names

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

void function MulTE_Estimates::post(string scalar b, string scalar V,| string scalar vce)
{
    real scalar uselevels
    real matrix best, sest
    string vector rownames, eqnames

    uselevels = any(st_vartype(Tvar) :== ("byte", "int", "long"))
    uselevels = uselevels | all(floor(Tvalues) :== Tvalues)
    if ( uselevels ) {
        rownames    = strofreal(Tvalues)
        rownames[1] = rownames[1] + "b"
        rownames    = rownames :+ "." :+ Tvar
    }
    else {
        printf("{bf:note:} table coefficients correspond to treatment values\n")
        rownames = Tlabels
    }
    rownames = J(cols(est), 1, rownames)
    best     = J(1, cols(est), 0) \ est
    eqnames  = rowshape(J(1, rows(best), colnames), 1)'

    if ( vce == "oracle" ) {
        // changed sest = J(1, cols(se_or), 0) \ se_or to sest = or_vcov
        sest = or_vcov
    }
    else {
        // changed sest = J(1, cols(se_po), 0) \ se_po to sest = po_vcov
        sest = po_vcov
    }

    st_matrix(b, rowshape(best', 1))
    st_matrixcolstripe(b, (eqnames, rownames))
    st_matrixrowstripe(b, ("", Yvar))

    // changed st_matrix(V, diag(rowshape(sest', 1):^2)) to: st_matrix(V, sest) 
    st_matrix(V, sest)
    st_matrixcolstripe(V, (eqnames, rownames))
    st_matrixrowstripe(V, (eqnames, rownames))
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
    colnames = ("Coef", "Own Effect", "Bias", "Min Bias", "Max Bias")'
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

void function MulTE_Decomposition::print(| real scalar minmax, real scalar digits)
{
    real scalar i, j, kcolprint
    real vector lengths
    string scalar colsep
    string matrix fmt_res, fmt_est, fmt_se
    string vector rownames, formats, formats2

    rownames = Tlabels[2::length(Tvalues)], J(rows(est), 1, "")
    rownames = rowshape(rownames, rows(est) * 2)

    if ( args() < 1 ) minmax = 0
    if ( args() < 2 ) digits = 6

    fmt_est = J(rows(est), cols(est), "")
    fmt_se  = J(rows(est), cols(est), "")

    for (i = 1; i <= rows(est); i++) {
        for (j = 1; j <= cols(est); j++) {
            fmt_est[i, j] = strtrim(sprintf("%21." + strofreal(digits) + "f", est[i, j]))
            fmt_se[i, j]  = "(" + strtrim(sprintf("%21." + strofreal(digits) + "f", se[i, j])) + ")"
        }
    }

    fmt_res  = rowshape((fmt_est, fmt_se), rows(est) * 2)
    lengths  = max((strlen(rownames) \ strlen(Tvar))), colmax(strlen(colnames' \ fmt_res))
    lengths  = lengths :+ 3
    formats  = " %" :+ strofreal(lengths) :+ "s"
    formats2 = "-%" :+ strofreal(lengths) :+ "s"

    colsep    = "|"
    colsep    = " "
    colsep    = ""
    kcolprint = minmax? length(colnames): (length(colnames) - 2)
    printf("\nContamination Bias Decomposition\n")
    printf("-%s\n", "-" * (sum(lengths[1::(kcolprint+1)]) + (strlen(colsep) + 1) * (kcolprint + 1) - 1))
    printf(colsep)

    printf(formats[1], Tvar)
    for (j = 1; j <= kcolprint; j++) {
        printf(colsep)
        printf(formats[j + 1], colnames[j])
    }
    // printf("|\n")
    printf("\n")

    for (j = 1; j <= (kcolprint + 1); j++) {
        printf(colsep)
        printf(formats2[j], "-" * lengths[j])
    }
    // printf("|\n")
    printf("\n")

    for (i = 1; i <= rows(fmt_res); i++) {
        printf(colsep)
        printf(formats[1], rownames[i])
        for (j = 1; j <= kcolprint; j++) {
            printf(colsep)
            printf(formats[j + 1], fmt_res[i, j])
        }
        // printf("|\n")
        printf("\n")
    }
    printf("-%s\n", "-" * (sum(lengths[1::(kcolprint+1)]) + (strlen(colsep) + 1) * (kcolprint + 1) - 1))
    printf("SE in parentheses; bias estimates stored in e(decomposition).\n")
}

real matrix function MulTE_Decomposition::lambda(real matrix Wm)
{
    return(Wm * delta_pr)
}

real matrix function MulTE_Decomposition::tauhat(real matrix Wm)
{
    return(Wm * gammam')
}
end

cap mata mata drop MulTE()
cap mata mata drop MulTE_Estimates()
cap mata mata drop MulTE_Decomposition()

mata
// Yvar  = "`depvar'"
// Tvar  = "`treatment'"
// Wvar  = "`W'"
// touse = "`touse'"
// wgt   = "`wgt'"
// wtype = "`weight'"

class MulTE
{
    class MulTE_Estimates scalar estimates
    class MulTE_Decomposition scalar decomposition

    real scalar cache
    real vector Y
    real vector w
    real vector xlevels
    real matrix Xm
    real matrix Wm

    void new()
    void cache_load()
    void cache_drop()
    void estimates()
    void decomposition()
}

class MulTE_Estimates
{
    real scalar run
    real scalar n
    real scalar nw
    real scalar k
    real matrix est
    real matrix se_po
    real matrix se_or
    real   vector Tvalues
    string vector Tlabels
    string vector colnames
    string scalar Tvar
    string scalar Yvar
    string scalar wgtvar
    string scalar wtype

    real matrix po_vcov
    real matrix or_vcov

    void new()
    void print()
    void save()
    void post()
}

class MulTE_Decomposition
{
    real scalar run
    real scalar n
    real scalar nw
    real scalar k
    real matrix est
    real matrix se
    real matrix tmp
    real matrix gammam
    real matrix delta_pr
    string scalar wgtvar
    string scalar wtype

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

// ----------------------------------------------------------------- //
//                           Main function                           //
// ----------------------------------------------------------------- //

void function MulTE::new()
{
    cache = 0
}

void MulTE::estimates(
    string scalar Yvar,
    string scalar Tvar,
    string scalar Wvar,
    string scalar touse,
    string scalar wgt,
    string scalar wtype)
{
    struct multe_helper_results scalar rk
    string scalar Tlab

    real scalar j, nobs, nobsw, k
    real vector X0, alpha0, lam, rcalpha, rcres, ts, s0, Wmean
    real vector Y_s, X0_s, w_s
    real matrix Xm_s, Wm_s
    real matrix psi_alpha0, ps, Xt
    real matrix est, se_or, se_po
    real matrix alphak, psi_alphak, psi_or, psi_po
    real matrix psi_pom, psi_orm, var_po_onem, var_or_onem, psi_pom_tl, po_vcov, psi_orm_tl, or_vcov
    real vector s, Xdot, eps, sk

    if ( this.estimates.run ) return

    // -----------------------------------------------------------------
    // Setup
    // -----------------------------------------------------------------

    cache_load(Yvar, Tvar, Wvar, touse, wgt)
    nobs   = length(Y)
    if ( wtype == "aweight" ) w = nobs * w / sum(w)
    nobsw  = (wtype == "fweight")? sum(w): nobs
    k      = cols(Xm)
    X0     = Xm[., 1]
    alpha0 = multe_helper_olsw(select(Y, X0), select(Wm, X0), select(w, X0))

    psi_alpha0  = (X0 :* Wm :* (Y - Wm * alpha0)) * qrinv(cross(X0 :* Wm, w, Wm)/nobsw)
    est = se_or = se_po = J(k - 1, 3, .)

    // common weights
    ps        = Wm * multe_helper_olsw(Xm, Wm, w) // propensity scores
    Xt        = Xm - ps                           // residuals
    lam       = 1 :/ rowsum(1:/ps)
    rcalpha   = multe_helper_olsw(Y, Xm, w :* lam :/ rowsum(ps :* Xm))
    rcres     = Y - Xm * rcalpha
    est[., 3] = rcalpha[2::length(rcalpha)] :- rcalpha[1]

    // Fitted t values
    ts = Wm * multe_helper_olsw(rcres :* Xm :/ (ps:^2)  :* lam, Wm, w)
    s0 = Wm * multe_helper_olsw(rcres :* X0 :/ ps[., 1] :* (lam:^2) :/ (ps:^2), Wm, w)

    // -----------------------------------------------------------------
    // TE estimates
    // -----------------------------------------------------------------

    var_po_onem = var_or_onem = J(k, 1, 0)
    psi_pom = psi_orm = J(nobs, k*3, 0)
    Wmean = mean(Wm, w)

    for(j = 2; j <= k; j++) {
        alphak = multe_helper_olsw(select(Y, Xm[., j]), select(Wm, Xm[., j]), select(w, Xm[., j]))
        psi_alphak = (Xm[.,j] :* (Y - Wm * alphak) :* Wm) * qrinv(cross(Xm[., j] :* Wm, w, Wm)/nobsw)

        // ATE
        est[j - 1, 1] = (Wmean * (alphak - alpha0))
        psi_or        = ((psi_alphak - psi_alpha0) * Wmean')
        psi_orm[.,j]  = psi_or
        psi_pom[.,j]  = psi_or + ((Wm :-  Wmean) * (alphak - alpha0))
        if ( wtype == "fweight" ) {
            se_or[j - 1, 1] = sqrt(variance(psi_or, w) * (nobsw - 1) / nobsw^2)
            se_po[j - 1, 1] = sqrt(variance(psi_pom[.,j], w) * (nobsw - 1) / nobsw^2)
        }
        else {
            se_or[j - 1, 1] = sqrt(variance(w :* psi_or) * (nobsw - 1) / nobsw^2)
            se_po[j - 1, 1] = sqrt(variance(w :* psi_pom[.,j]) * (nobsw - 1) / nobsw^2)
        }

        // One treatment at a time
        s    = selectindex(X0 :| Xm[., j])
        Xm_s = Xm[s, j]
        Wm_s = Wm[s, .]
        X0_s = X0[s]
        Y_s  = Y[s]
        w_s  = w[s]
        if ( wtype == "aweight" ) w_s = length(Y_s) * w_s / sum(w_s)

        Xdot            = Xm_s - Wm_s * multe_helper_olsw(Xm_s, Wm_s, w_s)
        rk              = multe_helper_olswr(Y_s, (Xdot, Wm_s), w_s)
        est[j - 1, 2]   = rk.coefficients[1]
        eps             = Xm_s :* (Y_s - Wm_s * alphak) +
                          X0_s :* (Y_s - Wm_s * alpha0)
        if ( wtype == "fweight" ) {
            se_po[j - 1, 2] = sqrt(sum(w_s :* (rk.residuals:^2) :* (Xdot:^2)) / sum(w_s :* (Xdot:^2)):^2)
            se_or[j - 1, 2] = sqrt(sum(w_s :* (eps:^2)          :* (Xdot:^2)) / sum(w_s :* (Xdot:^2))^2)
            var_po_onem[j, 1] = sum(w_s :* (rk.residuals:^2) :* (Xdot:^2)) / sum(w_s :* (Xdot:^2)):^2
            var_or_onem[j, 1] = sum(w_s :* (eps:^2)          :* (Xdot:^2)) / sum(w_s :* (Xdot:^2))^2
        }
        else {
            se_po[j - 1, 2] = sqrt(sum((rk.residuals:^2) :* ((w_s :* Xdot):^2)) / sum(w_s :* (Xdot:^2)):^2)
            se_or[j - 1, 2] = sqrt(sum((eps:^2)          :* ((w_s :* Xdot):^2)) / sum(w_s :* (Xdot:^2))^2)
            var_po_onem[j, 1] = sum((rk.residuals:^2) :* ((w_s :* Xdot):^2)) / sum(w_s :* (Xdot:^2)):^2
            var_or_onem[j, 1] = sum((eps:^2)          :* ((w_s :* Xdot):^2)) / sum(w_s :* (Xdot:^2))^2
        }

        // common weights
        psi_or = lam :* (Xm[., j] :* (Y - Wm * alphak) :/ ps[., j] -
                 X0 :* (Y - Wm * alpha0) :/ ps[., 1]) :/ mean(lam, w)
        sk = Wm * multe_helper_olsw(rcres :* Xm[., j] :/ ps[., j] :* (lam:^2) :/ (ps:^2), Wm, w)
        psi_po = (lam :* rcres :* (Xm[., j] :/ ps[., j] - X0 :/ ps[., 1]) +
                  Xt[., 1] :* ts[., 1] - Xt[., j] :* ts[., j] + rowsum(Xt :* (sk - s0))) / mean(lam, w)
        if ( wtype == "fweight" ) {
            se_po[j-1, 3] = sqrt(variance(psi_po, w)*(nobsw-1)/nobsw^2)
            se_or[j-1, 3] = sqrt(variance(psi_or, w)*(nobsw-1)/nobsw^2)
        }
        else {
            se_po[j-1, 3] = sqrt(variance(psi_po :* w)*(nobsw-1)/nobsw^2)
            se_or[j-1, 3] = sqrt(variance(psi_or :* w)*(nobsw-1)/nobsw^2)
        }
        psi_pom[., j+(k*2)] = psi_po
        psi_orm[., j+(k*2)] = psi_or
    }

    // Compute vcov matrices
    // NOTE: Var(beta) = Var(psi)/n. That's why po_vcov has n^2 in the denominator and not n.
    if ( wtype == "fweight" ) {
        psi_pom_tl = psi_pom :- mean(psi_pom, w)
        psi_orm_tl = psi_orm :- mean(psi_orm, w)
        po_vcov    = cross(psi_pom_tl, w, psi_pom_tl) / (nobsw^2)
        or_vcov    = cross(psi_orm_tl, w, psi_orm_tl) / (nobsw^2)
    }
    else {
        psi_pom    = psi_pom :* w
        psi_orm    = psi_orm :* w
        psi_pom_tl = psi_pom :- mean(psi_pom)
        psi_orm_tl = psi_orm :- mean(psi_orm)
        po_vcov    = cross(psi_pom_tl, psi_pom_tl) / (nobsw^2)
        or_vcov    = cross(psi_orm_tl, psi_orm_tl) / (nobsw^2)
    }

    for(j=2; j<=k; j++) {
        po_vcov[j+k, j+k] = var_po_onem[j,1]
        or_vcov[j+k, j+k] = var_or_onem[j,1]
    }

    po_vcov = blockdiag(po_vcov[|1,1 \ k, k|], blockdiag(po_vcov[|k+1, k+1 \ 2*k, 2*k|], po_vcov[|2*k+1, 2*k+1 \ 3*k, 3*k|]))
    or_vcov = blockdiag(or_vcov[|1,1 \ k, k|], blockdiag(or_vcov[|k+1, k+1 \ 2*k, 2*k|], or_vcov[|2*k+1, 2*k+1 \ 3*k, 3*k|]))

    Tlab = st_varvaluelabel(Tvar)
    this.estimates.n        = nobs
    this.estimates.nw       = nobsw
    this.estimates.k        = k
    this.estimates.est      = est
    this.estimates.se_po    = se_po
    this.estimates.se_or    = se_or
    this.estimates.Tvalues  = xlevels
    this.estimates.Tlabels  = Tlab == ""? strofreal(xlevels): st_vlmap(Tlab, xlevels)
    this.estimates.Tvar     = Tvar
    this.estimates.Yvar     = Yvar
    this.estimates.wgtvar   = wgt
    this.estimates.wtype    = wtype
    this.estimates.po_vcov  = po_vcov
    this.estimates.or_vcov  = or_vcov
    this.estimates.run      = 1
}

void function MulTE::cache_load(
    string scalar Yvar,
    string scalar Tvar,
    string scalar Wvar,
    string scalar touse,
    string scalar wgt)
{
    real scalar j
    real vector X

    if ( cache ) return

    Wm = designmatrix(st_data(., Wvar, touse))
    Y  = st_data(., Yvar, touse)
    X  = st_data(., Tvar, touse)
    w  = wgt == ""? J(length(Y), 1, 1): st_data(., wgt, touse)
    xlevels = uniqrows(X)
    Xm = J(rows(Y), length(xlevels), 0)
    for (j = 1; j <= length(xlevels); j++) {
        Xm[., j] = (X :== xlevels[j])
    }

    cache = 1
}

void function MulTE::cache_drop()
{
    Y       = .
    Xm      = .
    Wm      = .
    xlevels = .
    cache   = 0
}

void function MulTE::decomposition(
    string scalar Yvar,
    string scalar Tvar,
    string scalar Wvar,
    string scalar touse,
    string scalar wgt,
    string scalar wtype)
{
    struct multe_helper_results scalar ri, rd, rddX
    string vector lambda_names, tauhat_names

    real scalar i, l, j, nobs, nobsw, k, kw, kwx, j1, j2
    real matrix WmXm, gamma, psi_gamma, deltak, psi_deltak, ddX, M, psi
    real matrix est, estk, se, psimin, psimax
    real matrix delta_kl, delta_pr, gammam
    real matrix ggi, ggd, psigi, psigd
    real vector ghelper, gi, gd, di

    if ( this.decomposition.run ) return

    cache_load(Yvar, Tvar, Wvar, touse, wgt)
    nobs  = length(Y)
    if ( wtype == "aweight" ) w = nobs * w / sum(w)
    nobsw = (wtype == "fweight")? sum(w): nobs
    k     = cols(Xm)
    kw    = cols(Wm)

    // From R> matrix of controls Wm including intercept; X must be a
    // From R> factor, first level to be dropped

    // NB> Wm does not have a constant because it contains the full
    // NB> matrix of dummies.

    kwx  = kw + (k - 1) * kw
    M    = I(k - 1)#J(kw, 1, 1)
    Xm   = Xm[., 2::k]
    WmXm = J(nobs, kwx, 0)
    WmXm[|1, 1 \ nobs, kw|] = Wm
    for(j = 1; j < k; j++) {
        j1 = kw + (j - 1) * kw + 1
        j2 = kw + (j - 0) * kw
        WmXm[|1, j1 \ nobs, j2|] = Xm[., j] :* Wm
    }

    ri = multe_helper_olswr(Y, WmXm, w)
    gamma = ri.coefficients[|(kw+1) \ cols(WmXm)|]
    rd = multe_helper_olswr(WmXm[|(1, kw+1) \ (nobs, cols(WmXm))|], (Xm, Wm), w)
    psi_gamma = ((ri.residuals :* WmXm) *
                  invsym(cross(WmXm, w, WmXm)))[|(1, kw+1) \ (nobs, cols(WmXm))|]

    // NB: Given Wm and Xm are collections of non-overlapping
    // indicators, cross(WmXm, WmXm) and subsequent calculations
    // simplify. In teting, however, this was not necessarily faster.

    // Sort columns by size
    ghelper = rowshape(colshape(J(kw, 1, 1::(k-1)), k-1)', kw * (k - 1))
    gi      = order((ghelper, gamma), (1, 2))
    gd      = order((ghelper, gamma), (1, -2))
    est     = se = J(k-1, 5, .)
    ggi     = gamma[gi] :* M
    ggd     = gamma[gd] :* M
    psigi   = psi_gamma[., gi]
    psigd   = psi_gamma[., gd]

    // Standard errors
    for (j = 1; j < k; j++) {
        rddX   = multe_helper_olswr(Xm[., j], (multe_helper_antiselect(Xm, j), Wm), w)
        ddX    = rddX.residuals // ddot(X)
        deltak = rd.coefficients[j, .]'
        if ( wtype == "fweight" ) {
            psi_deltak = ddX :* rd.residuals / sum(w :* (ddX:^2))
        }
        else {
            psi_deltak = (w :* ddX) :* rd.residuals / sum(w :* (ddX:^2))
        }

        psi  = psi_deltak * (gamma :* M) + psi_gamma * (deltak :* M)
        estk = gamma' * (deltak :* M)

        // Sort delta columns
        di = order((ghelper, deltak), (1, 2))

        psi_deltak = psi_deltak[., di]
        deltak = deltak[di] :* M
        psimax = psi_deltak * ggi + psigi * deltak
        psimin = psi_deltak * ggd + psigd * deltak

        // NB: These are coefficients and needn't be weighted summed
        est[j,.] = (
            sum(estk),
            estk[j],
            sum(multe_helper_antiselect(estk, j)),
            sum(ggd' * multe_helper_antiselect(deltak, j)),
            sum(ggi' * multe_helper_antiselect(deltak, j))
        )

        if ( wtype == "fweight" ) {
            se[j,.] = sqrt((
                sum(w :* (rowsum(psi):^2)),
                sum(w :* (psi[., j]:^2)),
                sum(w :* (rowsum(multe_helper_antiselect(psi, j)):^2)),
                sum(w :* (rowsum(multe_helper_antiselect(psimin, j)):^2)),
                sum(w :* (rowsum(multe_helper_antiselect(psimax, j)):^2))
            ))
        }
        else {
            se[j,.] = sqrt((
                sum(rowsum(psi):^2),
                sum(psi[., j]:^2),
                sum(rowsum(multe_helper_antiselect(psi, j)):^2),
                sum(rowsum(multe_helper_antiselect(psimin, j)):^2),
                sum(rowsum(multe_helper_antiselect(psimax, j)):^2)
            ))
        }
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
    delta_pr = delta_kl :/ (mean(Wm, w)')
    gammam   = rowshape(gamma, k-1)

    this.decomposition.n            = nobs
    this.decomposition.nw           = nobsw
    this.decomposition.k            = k
    this.decomposition.est          = est
    this.decomposition.se           = se
    this.decomposition.tmp          = rowshape((est, se), 2 * rows(est))
    this.decomposition.Tvalues      = this.estimates.Tvalues
    this.decomposition.Tlabels      = this.estimates.Tlabels
    this.decomposition.Tvar         = this.estimates.Tvar
    this.decomposition.Yvar         = this.estimates.Yvar
    this.decomposition.wgtvar       = wgt
    this.decomposition.wtype        = wtype
    this.decomposition.delta_pr     = delta_pr
    this.decomposition.gammam       = gammam
    this.decomposition.tauhat_names = tauhat_names
    this.decomposition.lambda_names = lambda_names
    this.decomposition.run          = 1
}

// ----------------------------------------------------------------- //
//                         Estimates Helpers                         //
// ----------------------------------------------------------------- //

void function MulTE_Estimates::new()
{
    run = 0
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
        sest = or_vcov
    }
    else {
        sest = po_vcov
    }

    st_matrix(b, rowshape(best', 1))
    st_matrixcolstripe(b, (eqnames, rownames))
    st_matrixrowstripe(b, ("", Yvar))

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

// ----------------------------------------------------------------- //
//                       Decomposition Helpers                       //
// ----------------------------------------------------------------- //

void function MulTE_Decomposition::new()
{
    run= 0
    colnames = ("Coef", "Own Effect", "Bias", "Min Bias", "Max Bias")'
}

void function MulTE_Decomposition::save(string scalar outmatrix)
{
    string vector rownames

    if ( run == 0 ) return
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

    if ( run == 0 ) return
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

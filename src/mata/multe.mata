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

    real scalar linear
    real scalar cache
    real vector Y
    real vector w
    real vector xlevels
    real matrix Xm
    real matrix Wm
    real vector omitW
    real vector omitTW
    string vector varTW
    string matrix controls

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
    cache  = 0
    linear = st_local("linear") != ""
    omitW  = linear? st_matrix(st_local("Womit")): .
    omitTW = st_matrix(st_local("TWomit"))
    varTW  = tokens(st_local("TW"))'
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

    real scalar j, nobs, nobsw, k, kw
    real vector X0, alpha0, lam, rcalpha, rcres, ts, s0, Wmean
    real vector Y_s, X0_s, w_s
    real matrix Xm_s, Wm_s
    real matrix psi_alpha0, ps, Xt
    real matrix est
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
    kw     = cols(Wm)
    X0     = Xm[., 1]
    s      = selectindex(X0)
    Wm_s   = Wm[s,.]
    w_s    = w[s]
    rk     = multe_helper_olswr(Y[s], Wm_s, w_s)
    alpha0 = rk.coefficients

    psi_alpha0      = J(nobs, kw, 0)
    psi_alpha0[s,.] = (rk.residuals :* Wm_s) * invsym(quadcross(Wm_s, w_s, Wm_s))
    est = J(k - 1, 3, .)

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
        s      = selectindex(Xm[., j])
        Wm_s   = Wm[s,.]
        w_s    = w[s]
        rk     = multe_helper_olswr(Y[s], Wm_s, w_s)
        alphak = rk.coefficients

        psi_alphak      = J(nobs, kw, 0)
        psi_alphak[s,.] = (rk.residuals :* Wm_s) * invsym(quadcross(Wm_s, w_s, Wm_s))

        // ATE
        est[j - 1, 1] = (Wmean * (alphak - alpha0))
        psi_or        = (psi_alphak - psi_alpha0) * Wmean'
        psi_pom[.,j]  = nobsw * (psi_or + ((Wm :- Wmean) * (alphak - alpha0))/sum(w))
        psi_orm[.,j]  = nobsw * psi_or
        // NB: Both psi_or and psi_pom should have weighted mean 0:
        //     - For psi_alphak, psi_alpha0, the residuals have to be orthogonal
        //       to Wm_s mechanically, so the weighted sum must be 0
        //     - Since (Wm :- Wmean) is centered, it's weighted mean is 0 as well

        // One treatment at a time
        s    = selectindex(X0 :| Xm[., j])
        Xm_s = Xm[s, j]
        Wm_s = Wm[s, .]
        X0_s = X0[s]
        Y_s  = Y[s]
        w_s  = w[s]

        Xdot            = Xm_s - Wm_s * multe_helper_olsw(Xm_s, Wm_s, w_s)
        rk              = multe_helper_olswr(Y_s, (Xdot, Wm_s), w_s)
        est[j - 1, 2]   = rk.coefficients[1]
        eps             = Xm_s :* (Y_s - Wm_s * alphak) +
                          X0_s :* (Y_s - Wm_s * alpha0)
        if ( wtype == "fweight" ) {
            var_or_onem[j, 1] = sum(w_s :* (eps:^2)          :* (Xdot:^2)) / sum(w_s :* (Xdot:^2))^2
            var_po_onem[j, 1] = sum(w_s :* (rk.residuals:^2) :* (Xdot:^2)) / sum(w_s :* (Xdot:^2)):^2
        }
        else {
            var_or_onem[j, 1] = sum((eps:^2)          :* ((w_s :* Xdot):^2)) / sum(w_s :* (Xdot:^2))^2
            var_po_onem[j, 1] = sum((rk.residuals:^2) :* ((w_s :* Xdot):^2)) / sum(w_s :* (Xdot:^2)):^2
        }

        // common weights
        psi_or = lam :* (Xm[., j] :* (Y - Wm * alphak) :/ ps[., j] -
                 X0 :* (Y - Wm * alpha0) :/ ps[., 1]) :/ ((wtype == "fweight")? mean(lam, w): mean(lam :* w))
        sk = Wm * multe_helper_olsw(rcres :* Xm[., j] :/ ps[., j] :* (lam:^2) :/ (ps:^2), Wm, w)
        psi_po = (lam :* rcres :* (Xm[., j] :/ ps[., j] - X0 :/ ps[., 1]) +
                  Xt[., 1] :* ts[., 1] - Xt[., j] :* ts[., j] + rowsum(Xt :* (sk - s0))) /
                  ((wtype == "fweight")? mean(lam, w): mean(lam :* w))
        psi_pom[., j+(k*2)] = psi_po
        psi_orm[., j+(k*2)] = psi_or
    }

// xx if ( st_local("debug") != "" ) {
// xx     debug_common = colsum(Y :* w :* lam :* Xm :/ ps) :/ colsum(w :* lam :* Xm :/ ps)
// xx     debug_common = (debug_common :- debug_common[1])'
// xx
// xx     xx_ate = colsum(Y :* w :* Xm :/ ps) :/ colsum(w :* Xm :/ ps)
// xx     xx_ate = (xx_ate :- xx_ate[1])'
// xx
// xx     xx_ate = multe_helper_olsw(Y, Xm, w :/ rowsum(ps :* Xm))
// xx     xx_ate = (xx_ate :- xx_ate[1])
// xx
// xx     xx_oaat = J(k, 1, .)
// xx     for (j = 2; j <= k; j++)  {
// xx         s    = selectindex(X0 :| Xm[., j])
// xx         Y_s  = Y[s]
// xx         Xm_s = Xm[s, j]
// xx         Wm_s = Wm[s, .]
// xx         w_s  = w[s]
// xx         xx_oaat[j] = multe_helper_olsw(Y_s, (Xm_s, Wm_s), w_s)[1]
// xx     }
// xx
// xx     // lam = (ps :* ps[., 1]) :/ (ps :+ ps[., 1])
// xx     // debug_oaat = colsum(Y :* w :* lam :* Xm :/ ps) :/ colsum(w :* lam :* Xm :/ ps)
// xx     // debug_oaat = (debug_oaat :- debug_oaat[1])'
// xx     // debug_oaat
// xx     //
// xx     // debug_oaat = multe_helper_olsw(Y, Xm, w :* rowsum(lam :* Xm) :/ rowsum(ps :* Xm))
// xx     // debug_oaat = (debug_oaat :- debug_oaat[1])
// xx     // debug_oaat
// xx }
// xx
// xx debug_oaat = J(k, 1, .)
// xx debug_ate  = J(k, 1, .)
// xx for (j = 2; j <= k; j++)  {
// xx     s     = selectindex(X0 :| Xm[., j])
// xx     Xm_s  = Xm[s, j]
// xx     Wm_s  = Wm[s, .]
// xx     Y_s   = Y[s]
// xx     w_s   = w[s]
// xx     Xdot  = Xm_s - Wm_s * multe_helper_olsw(Xm_s, Wm_s, w_s)
// xx     alpha0 = multe_helper_olsw(select(Y, X0), select(Wm, X0), select(w, X0))
// xx     alphak = multe_helper_olsw(select(Y, Xm[., j]), select(Wm, Xm[., j]), select(w, Xm[., j]))
// xx     Wmean  = mean(Wm, w)
// xx     debug_oaat[j] = multe_helper_olsw(Y_s, (Xdot, Wm_s), w_s)[1]
// xx     debug_ate[j]  = (Wmean * (alphak - alpha0))
// xx }
// xx
// xx    debug_common[2::k]  :- est[.,3]
// xx    debug_ate  :- xx_ate
// xx    debug_oaat :- xx_oaat

    // Compute vcov matrices
    // NOTE: Var(beta) = Var(psi)/n. That's why po_vcov has n^2 in the denominator and not n.
    if ( wtype == "fweight" ) {
        psi_pom_tl = psi_pom :- mean(psi_pom, w)
        psi_orm_tl = psi_orm :- mean(psi_orm, w)
        po_vcov    = quadcross(psi_pom_tl, w, psi_pom_tl) / (nobsw^2)
        or_vcov    = quadcross(psi_orm_tl, w, psi_orm_tl) / (nobsw^2)
    }
    else {
        psi_pom_tl = psi_pom :* w :- mean(psi_pom, w)
        psi_orm_tl = psi_orm :* w :- mean(psi_orm, w)
        po_vcov    = quadcross(psi_pom_tl, psi_pom_tl) / (nobsw^2)
        or_vcov    = quadcross(psi_orm_tl, psi_orm_tl) / (nobsw^2)
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
    this.estimates.se_po    = sqrt((colshape(diagonal(po_vcov), k)')[2..k,.])
    this.estimates.se_or    = sqrt((colshape(diagonal(or_vcov), k)')[2..k,.])
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

    Y  = st_data(., Yvar, touse)
    X  = st_data(., Tvar, touse)
    w  = wgt == ""? J(length(Y), 1, 1): st_data(., wgt, touse)
    xlevels = uniqrows(X)
    Xm = J(rows(Y), length(xlevels), 0)
    for (j = 1; j <= length(xlevels); j++) {
        Xm[., j] = (X :== xlevels[j])
    }

    if ( linear ) {
        if ( length(omitW) ) {
            Wm = J(length(Y), 1, 1), select(st_data(., Wvar, touse), !omitW)
        }
        else {
            Wm = J(length(Y), 1, 1), st_data(., Wvar, touse)
        }
    }
    else {
        Wm = designmatrix(st_data(., Wvar, touse))
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
    struct multe_helper_results scalar rv, ru, rd, Xt, Xtt
    string vector lambda_names, tauhat_names

    real vector varord
    real scalar i, l, j, nobs, nobsw, k, kw, kwx, j1, j2, tol, qc
    real matrix Xj, WmXm, deltak, Wm_s
    real matrix Xdot, XdX, Xdotdot, est, se
    real matrix psi_beta, psi_alpha0, psi_alphak, psi_gammak, psi_deltakk, psi_min, psi_max, psi_ownk, psi_biask
    real matrix delta_kl, delta_pr
    real matrix ggi, ggd, psigi, psigd
    real vector alpha, gamma, gammak, licols, ghelper, gi, gd, di, selk, inv, dsel, gsel, flags
    real vector beta, beta_bias, beta_own, beta_min, beta_max
    real vector var_beta, var_bias, var_own, var_min, var_max

    if ( this.decomposition.run ) return

    this.cache_load(Yvar, Tvar, Wvar, touse, wgt)
    nobs  = length(Y)
    nobsw = (wtype == "fweight")? sum(w): nobs
    k     = cols(Xm)
    kw    = cols(Wm)
    Xj    = Xm[., 2..k]
    kwx   = k * kw
    tol   = epsilon(1)

    varord   = J(length(varTW), 1, .)
    for (j = 1; j <= length(varTW); j++) {
        if ( ustrregexm(varTW[j], "(\d).*?\.") ) {
            varord[j] = strtoreal(ustrregexs(1))
        }
    }

    if ( wtype == "aweight" ) w = w * (nobs / sum(w))

    // NB: W is the paper notation for a single categorical variable;
    // Z is for an arbitrary collection. Z requires a linearity assumption
    // Whereas for W it's not requried. Here I use W throughout for both.

    // WmXm = st_data(., st_local("TW"), touse)
    WmXm = J(nobs, kwx, 0)
    for(j = 1; j <= k; j++) {
        j1 = (j - 1) * kw + 1
        j2 = (j - 0) * kw
        WmXm[|1, j1 \ nobs, j2|] = Xm[., j] :* Wm
    }

    // You have to reshape because of the way Stata inevitably sorts
    // it in this weird way. This runs
    //     - Equation (10) in the note, Y_i = sum_k X_ik * W_i alpha_k + V_i
    //     - Equation (15) in the note, Y_i = X_i' * beta + W_i + U_i
    licols = omitTW[order((varord, (1::length(varord))), (1, 2))]
    rv     = multe_helper_olswr(Y, select(WmXm, !licols), w)
    ru     = multe_helper_olswr(Y, (Xm[., 2..k], Wm), w)
    rd     = multe_helper_olswr(WmXm[|(1, kw+1) \ (nobs, kwx)|], (Xj, Wm), w)
    beta   = ru.coefficients[1..(k-1)]
    WmXm   = .

    // It is useful to regress X on W to find Xdot
    Xt   = multe_helper_olswr(Xj, Wm, w)
    Xdot = Xt.residuals
    XdX  = invsym(quadcross(Xdot, w, Xdot))

    // We can also compute \hat{alpha}
    alpha = J(kwx, 1, .)
    alpha[selectindex(!licols)] = rv.coefficients
    gamma = (alpha[|(kw+1) \ kwx|] :- J(k-1, 1, alpha[|1 \ kw|]))

    // psi_beta is just for regular estimate
    psi_beta = (ru.residuals :* Xdot) * XdX

    // psi_gammak = psi_alphak - psi_alpha0
    selk = selectindex(Xm[.,1])
    Wm_s = Wm[selk, .]
    psi_alpha0 = J(nobs, kw, 0)
    psi_alpha0[selk,.] = (rv.residuals[selk] :* Wm[selk, .]) * invsym(quadcross(Wm_s, w[selk,.], Wm_s))
    psi_gammak = J(nobs, kw * (k - 1), 0)
    for(j = 2; j <= k; j++) {
        j1   = (j - 1) * kw + 1
        j2   = (j - 0) * kw
        selk = selectindex(Xm[.,j])
        Wm_s = Wm[selk, .]
        psi_alphak = J(nobs, kw, 0)
        psi_alphak[selk,.] = (rv.residuals[selk] :* Wm[selk, .]) * invsym(quadcross(Wm_s, w[selk], Wm_s))
        psi_gammak[|1,j1-kw \ nobs,j2-kw|] = psi_alphak :- psi_alpha0
    }

    psi_alphak = .
    psi_alpha0 = .

    // min and max bias helpers
    // xx [take out worst-case for linear?] if ( linear == 0 ) {
        ghelper = rowshape(colshape(J(kw, 1, 1::(k-1)), k-1)', kw * (k - 1))
        gd      = order((ghelper, gamma), (1, -2))
        gi      = order((ghelper, gamma), (1, 2))
        ggd     = gamma[gd]
        ggi     = gamma[gi]
        psigd   = psi_gammak[., gd]
        psigi   = psi_gammak[., gi]
    // xx [take out worst-case for linear?] }

    // deltakk and psi_alphak
    beta_min = J(k-1,1,.)
    beta_max = J(k-1,1,.)
    beta_own = J(k-1,1,.)
    var_own  = J(k-1,1,.)
    var_bias = J(k-1,1,.)
    var_min  = J(k-1,1,.)
    var_max  = J(k-1,1,.)
    for(j = 2; j <= k; j++) {
        j1  = (j - 2) * kw + 1
        j2  = (j - 1) * kw
        inv = selectindex(!e(j-1, k-1))

        // deltak
        deltak = rd.coefficients[j-1,.]'
        gammak = gamma[|j1\j2|]
        flags  = (abs(deltak[|j1\j2|]) :< tol) :| (gammak :< .)
        if ( !all(flags) ) {
            var_own[j-1]  = .
            var_bias[j-1] = .
            beta_own[j-1] = .
            beta_min[j-1] = .
            beta_max[j-1] = .
            var_min[j-1]  = .
            var_max[j-1]  = .
        }
        else {
            // Own bias
            _editmissing(gammak, 0)
            beta_own[j-1] = quadcross(deltak[|j1\j2|], gammak)

            // alphak
            selk = selectindex(Xm[., j])
            Wm_s = Wm[selk, .]

            // X_k on X_-k and Z (Xdot has projected X on Z )
            Xtt     = multe_helper_olswr(Xdot[.,j-1], Xdot[., inv], w)
            Xdotdot = Xtt.residuals

            // SE via psi
            psi_deltakk = (Xdotdot :* rd.residuals) / quadsum(w :* (Xdotdot:^2))
            psi_ownk    = psi_gammak[|1,j1 \ nobs,j2|] * deltak[|j1\j2|] + psi_deltakk[|1,j1 \ nobs,j2|] * gammak
            psi_biask   = psi_beta[., j-1] :- psi_ownk

            var_own[j-1]  = (wtype == "fweight")? quadsum(w :* (psi_ownk:^2)):  quadsum((w :* psi_ownk):^2)
            var_bias[j-1] = (wtype == "fweight")? quadsum(w :* (psi_biask:^2)): quadsum((w :* psi_biask):^2)

            // Min and max bias, with SE
            // xx [take out worst-case for linear?] (linear == 0) |
            di   = order((ghelper, deltak), (1, 2))
            dsel = multe_helper_antiselect(di, j1::j2)
            gsel = multe_helper_antiselect(1::length(gamma), j1::j2)
            if ( missing(gamma[gsel]) ) {
                beta_min[j-1] = .
                beta_max[j-1] = .
                var_min[j-1]  = .
                var_max[j-1]  = .
            }
            else {
                deltak        = deltak[dsel]
                psi_deltakk   = psi_deltakk[., dsel]
                beta_min[j-1] = quadsum(ggd[gsel] :* deltak)
                beta_max[j-1] = quadsum(ggi[gsel] :* deltak)
                psi_min       = psigd[., gsel] * deltak :+ psi_deltakk * ggd[gsel]
                psi_max       = psigi[., gsel] * deltak :+ psi_deltakk * ggi[gsel]
                var_min[j-1]  = (wtype == "fweight")? quadsum(w :* (psi_min:^2)): quadsum((w :* psi_min):^2)
                var_max[j-1]  = (wtype == "fweight")? quadsum(w :* (psi_max:^2)): quadsum((w :* psi_max):^2)
            }
        }
    }

    qc        = nobs / (nobs - (k + kw - 1))
    beta_bias = beta - beta_own
    var_beta  = (wtype == "fweight")? quadcolsum(w :* (psi_beta:^2)): quadcolsum((w :* psi_beta):^2)

    est = beta, beta_own, beta_bias, beta_min, beta_max
    se  = sqrt((var_beta', var_own, var_bias, var_min, var_max))

    if ( missing(est[., 2..3]) ) {
        if ( linear ) {
            printf("(note: unable to estimate bias for all groups; collinear covariates within treatments)\n")
        }
        else {
            printf("(note: unable to estimate bias for all groups; consider option {opt overlap})\n")
        }
    }

    // Direct computation of bias
    // lam = Xdot * XdX
    // tau = Wm * gammam'
    // debug_own  = colsum(lam :* tau :* Xj :* w)'
    // debug_bias = colsum((lam :* (1 :- Xj)) :* rowsum(tau :* Xj) :* w)'

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
    delta_pr = delta_kl :/ mean(Wm, w)'

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
    this.decomposition.gammam       = rowshape(gamma, k-1)'
    this.decomposition.tauhat_names = tauhat_names
    this.decomposition.lambda_names = lambda_names
    this.decomposition.run          = 1
}

// xx void function MulTE::decomposition(
// xx     string scalar Yvar,
// xx     string scalar Tvar,
// xx     string scalar Wvar,
// xx     string scalar touse,
// xx     string scalar wgt,
// xx     string scalar wtype)
// xx {
// xx     struct multe_helper_results scalar ri, rd, rddX, check
// xx     string vector lambda_names, tauhat_names
// xx
// xx     real scalar i, l, j, nobs, nobsw, k, kw, kwx, j1, j2
// xx     real matrix WmXm, gamma, psi_gamma, deltak, psi_deltak, ddX, M, psi
// xx     real matrix est, estk, se, psimin, psimax
// xx     real matrix delta_kl, delta_pr, gammam
// xx     real matrix ggi, ggd, psigi, psigd
// xx     real matrix lam, tau
// xx     real vector ghelper, gi, gd, di
// xx
// xx     if ( this.decomposition.run ) return
// xx
// xx     cache_load(Yvar, Tvar, Wvar, touse, wgt)
// xx     nobs  = length(Y)
// xx     if ( wtype == "aweight" ) w = nobs * w / sum(w)
// xx     nobsw = (wtype == "fweight")? sum(w): nobs
// xx     k     = cols(Xm)
// xx     kw    = cols(Wm)
// xx
// xx     // From R> matrix of controls Wm including intercept; X must be a
// xx     // From R> factor, first level to be dropped
// xx
// xx     // NB> Wm does not have a constant because it contains the full
// xx     // NB> matrix of dummies.
// xx
// xx     kwx  = kw + (k - 1) * kw
// xx     M    = I(k - 1)#J(kw, 1, 1)
// xx     Xm   = Xm[., 2::k]
// xx     WmXm = J(nobs, kwx, 0)
// xx     WmXm[|1, 1 \ nobs, kw|] = Wm
// xx     for(j = 1; j < k; j++) {
// xx         j1 = kw + (j - 1) * kw + 1
// xx         j2 = kw + (j - 0) * kw
// xx         WmXm[|1, j1 \ nobs, j2|] = Xm[., j] :* Wm
// xx     }
// xx
// xx     ri = multe_helper_olswr(Y, WmXm, w)
// xx     gamma = ri.coefficients[|(kw+1) \ cols(WmXm)|]
// xx     rd = multe_helper_olswr(WmXm[|(1, kw+1) \ (nobs, cols(WmXm))|], (Xm, Wm), w)
// xx     psi_gamma = ((ri.residuals :* WmXm) *
// xx                   invsym(cross(WmXm, w, WmXm)))[|(1, kw+1) \ (nobs, cols(WmXm))|]
// xx
// xx     // Check using equations 25-27 directly
// xx     if ( st_local("debug") != "" ) {
// xx         check = multe_helper_olswr(Xm, Wm, w)
// xx         lam   = check.residuals * invsym(cross(check.residuals, w, check.residuals))
// xx         tau   = Wm * colshape(gamma, kw)'
// xx         this.decomposition.debug_own  = colsum(lam :* tau :* Xm :* w)'
// xx         this.decomposition.debug_bias = colsum((lam :* (1 :- Xm)) :* rowsum(tau :* Xm) :* w)'
// xx     }
// xx
// xx     // NB: Given Wm and Xm are collections of non-overlapping
// xx     // indicators, cross(WmXm, WmXm) and subsequent calculations
// xx     // simplify. In teting, however, this was not necessarily faster.
// xx
// xx     // Sort columns by size
// xx     ghelper = rowshape(colshape(J(kw, 1, 1::(k-1)), k-1)', kw * (k - 1))
// xx     gi      = order((ghelper, gamma), (1, 2))
// xx     gd      = order((ghelper, gamma), (1, -2))
// xx     est     = se = J(k-1, 5, .)
// xx     ggi     = gamma[gi] :* M
// xx     ggd     = gamma[gd] :* M
// xx     psigi   = psi_gamma[., gi]
// xx     psigd   = psi_gamma[., gd]
// xx
// xx     // Standard errors
// xx     for (j = 1; j < k; j++) {
// xx         rddX   = multe_helper_olswr(Xm[., j], (multe_helper_antiselect(Xm, j), Wm), w)
// xx         ddX    = rddX.residuals // ddot(X)
// xx         deltak = rd.coefficients[j, .]'
// xx         if ( wtype == "fweight" ) {
// xx             psi_deltak = ddX :* rd.residuals / sum(w :* (ddX:^2))
// xx         }
// xx         else {
// xx             psi_deltak = (w :* ddX) :* rd.residuals / sum(w :* (ddX:^2))
// xx         }
// xx
// xx         psi  = psi_deltak * (gamma :* M) + psi_gamma * (deltak :* M)
// xx         estk = gamma' * (deltak :* M)
// xx
// xx         // Sort delta columns
// xx         di = order((ghelper, deltak), (1, 2))
// xx
// xx         psi_deltak = psi_deltak[., di]
// xx         deltak = deltak[di] :* M
// xx         psimax = psi_deltak * ggi + psigi * deltak
// xx         psimin = psi_deltak * ggd + psigd * deltak
// xx
// xx         // NB: These are coefficients and needn't be weighted summed
// xx         est[j,.] = (
// xx             sum(estk),
// xx             estk[j],
// xx             sum(multe_helper_antiselect(estk, j)),
// xx             sum(ggd' * multe_helper_antiselect(deltak, j)),
// xx             sum(ggi' * multe_helper_antiselect(deltak, j))
// xx         )
// xx
// xx         if ( wtype == "fweight" ) {
// xx             se[j,.] = sqrt((
// xx                 sum(w :* (rowsum(psi):^2)),
// xx                 sum(w :* (psi[., j]:^2)),
// xx                 sum(w :* (rowsum(multe_helper_antiselect(psi, j)):^2)),
// xx                 sum(w :* (rowsum(multe_helper_antiselect(psimin, j)):^2)),
// xx                 sum(w :* (rowsum(multe_helper_antiselect(psimax, j)):^2))
// xx             ))
// xx         }
// xx         else {
// xx             se[j,.] = sqrt((
// xx                 sum(rowsum(psi):^2),
// xx                 sum(psi[., j]:^2),
// xx                 sum(rowsum(multe_helper_antiselect(psi, j)):^2),
// xx                 sum(rowsum(multe_helper_antiselect(psimin, j)):^2),
// xx                 sum(rowsum(multe_helper_antiselect(psimax, j)):^2)
// xx             ))
// xx         }
// xx     }
// xx
// xx     // Control-specific TEs and weights
// xx     //     - rd has coefficients of reg of X * W on X, W
// xx     //     - 1..(k-1) rows are coefs of X[., 2..k]
// xx     //     - reshape so lth set of 1..(k-1) columns has coefs
// xx     //       for X[., l+1] corresponding to X[., 2..k] * W
// xx     //       (i.e. X[., l+1] for X[., 2] * W, ..., X[., k] * W)
// xx
// xx     i = 0
// xx     tauhat_names = J(1, k-1, "")
// xx     lambda_names = J(1, (k-1)*(k-1), "")
// xx     for (l = 1; l < k; l++) {
// xx         tauhat_names[l] = sprintf("%g", l)
// xx         for (j = 1; j < k; j++) {
// xx             lambda_names[++i] = sprintf("%g%g", l, j)
// xx         }
// xx     }
// xx
// xx     delta_kl = rowshape(rowshape(rd.coefficients[1..(k-1),.], 1), (k-1)*(k-1))'
// xx     delta_pr = delta_kl :/ ((wtype == "fweight")? mean(Wm, w)': mean(Wm :* w)')
// xx     gammam   = rowshape(gamma, k-1)
// xx
// xx     this.decomposition.n            = nobs
// xx     this.decomposition.nw           = nobsw
// xx     this.decomposition.k            = k
// xx     this.decomposition.est          = est
// xx     this.decomposition.se           = se
// xx     this.decomposition.tmp          = rowshape((est, se), 2 * rows(est))
// xx     this.decomposition.Tvalues      = this.estimates.Tvalues
// xx     this.decomposition.Tlabels      = this.estimates.Tlabels
// xx     this.decomposition.Tvar         = this.estimates.Tvar
// xx     this.decomposition.Yvar         = this.estimates.Yvar
// xx     this.decomposition.wgtvar       = wgt
// xx     this.decomposition.wtype        = wtype
// xx     this.decomposition.delta_pr     = delta_pr
// xx     this.decomposition.gammam       = gammam
// xx     this.decomposition.tauhat_names = tauhat_names
// xx     this.decomposition.lambda_names = lambda_names
// xx     this.decomposition.run          = 1
// xx }

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

real matrix function MulTE_Decomposition::lambda(real matrix Xm, real matrix Wm, real vector w)
{
    struct multe_helper_results scalar Xt
    real vector Xdot
    real matrix XdX, Xj
    Xj   = Xm[., 2..cols(Xm)]
    Xt   = multe_helper_olswr(Xj, Wm, w)
    Xdot = Xt.residuals
    XdX  = invsym(quadcross(Xdot, w, Xdot))
    return(w :* Xdot[., (1..cols(Xj))#J(1, cols(Xj), 1)] :* J(1, cols(Xj), Xdot * XdX))
    // return(Wm * delta_pr)
}

real matrix function MulTE_Decomposition::tauhat(real matrix Wm)
{
    return(Wm * gammam)
}
end

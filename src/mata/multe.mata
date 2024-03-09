cap mata mata drop MulTE()
cap mata mata drop MulTE_Info()
cap mata mata drop MulTE_Return()
cap mata mata drop MulTE_Results()
cap mata mata drop MulTE_Tests()
cap mata mata drop multe_helper_post()

mata
class MulTE
{
    real scalar cache
    string rowvector labels

    struct MulTE_Info scalar info
    void new()
    void setup()
    void debug()
    struct MulTE_Results scalar decomposition()
}

struct MulTE_Tests
{
    real scalar stat
    real scalar df
    real scalar pval
}

struct MulTE_Return
{
    struct MulTE_Results scalar full, overlap
    real scalar has_overlap, Tk
    string scalar Tvar, Yvar
    string rowvector Tlevels
}

struct MulTE_Results
{
    string rowvector labels
    real matrix estA, estB, seO, seP, seB
    struct MulTE_Tests scalar LM, Wa
    real matrix Vpop_PL, Vpop_OWN, Vpop_ATE, Vpop_EW, Vpop_CW
    real matrix Vo_ATE, Vo_EW, Vo_CW
    real matrix Vdiff_OWN, Vdiff_ATE, Vdiff_EW, Vdiff_CW
    real vector touse
    real scalar N
}

struct MulTE_Info
{
    string scalar cluster
    string scalar touse
    string scalar wcall
    string scalar wgtvar
    string scalar Y
    string scalar T
    string scalar X
    string scalar Wi
    string scalar cons2
    string scalar zreg
    string scalar zindep

    real scalar Tk
    real scalar N
    real scalar Zk
    real scalar wgt
    real scalar zerotol
    real scalar cw_uniform
    real scalar debug

    real vector zomit
    real vector zomitbase
}

// ----------------------------------------------------------------- //
//                           Main function                           //
// ----------------------------------------------------------------- //

void function MulTE::new()
{
    cache  = 0
    labels = "PL", "OWN", "ATE", "EW", "CW"
    this.info.zerotol = 1e-6
    this.setup()
}

void function MulTE::setup()
{
    this.info.cluster    = st_local("cluster")
    this.info.touse      = st_local("touse")
    this.info.wcall      = st_local("wcall")
    this.info.wgtvar     = st_local("wgt")
    this.info.Y          = st_local("Y")
    this.info.T          = st_local("T")
    this.info.X          = st_local("X")
    this.info.Wi         = st_local("Wi")
    this.info.zreg       = st_local("zreg")
    this.info.zindep     = st_local("zindep")
    this.info.cons2      = st_local("cons2")
    this.info.zomit      = st_matrix(st_local("zomit"))
    this.info.zomitbase  = st_matrix(st_local("zomitbase"))
    this.info.Tk         = strtoreal(st_local("Tk"))
    this.info.wgt        = (st_local("wgt") != "")
    this.info.debug      = (st_local("debug") != "")
    this.info.cw_uniform = (st_local("cw_uniform") != "")
}

struct MulTE_Results scalar MulTE::decomposition()
{
    struct MulTE_Results scalar res
    struct multe_helper_results scalar rk, ro, dtXr, Xhat
    string scalar regopts, regcall, regtol, mtol, ir_start
    string vector mlfit
    real scalar i, j, wsum, ir_diff
    real vector Y, T, X0, C, ws, LM, Wa
    real vector rl, ri, Zb, s, w_s, Y_s, lam, ipi, cw
    real vector seli, selj, idx1, idx1n, th1, th, si
    real matrix Zm, Xf, tX, deltak, gamk, dtX, pp, Z_s, X_s
    real matrix psi_beta, psi_al, psi_po, psi_ownk, psi_k, psi_1, psi_or, psi_rk, psi_ri
    real matrix ir_omit, ir_gam, ir_nam, ir_ord, ir_X, ate0, ate
    real matrix pis, vpi, xf1, M0, M, xf1l, xf1r, a
    real matrix Sc, He, He1112, Vu, pis0, Scr, Her, Her1112, Vr

    res.labels = this.labels
    res.estA = res.estB = res.seP = res.seB = res.seO = J(this.info.Tk-1, 5, .)

    regopts = sprintf("%s if %s, noconstant", this.info.wcall, this.info.touse)
    regcall = sprintf("qui reg %s b1.%s %s %s", this.info.Y, this.info.T, this.info.zreg, regopts)

    this.debug("run base reg")
    stata(regcall)
    stata("tempvar rl Xf")
    stata(sprintf("predict %s, resid", st_local("rl")))

    this.info.N    = st_numscalar("e(N)")
    this.info.Zk   = length(st_matrix("e(b)")) - this.info.Tk
    res.estA[., 1] = st_matrix("e(b)")[2..this.info.Tk]'

    this.debug("read in data")
    rl   = st_data(., st_local("rl"), this.info.touse)
    Y    = st_data(., this.info.Y,  this.info.touse)
    T    = st_data(., this.info.T,  this.info.touse)
    X0   = T :== 1
    Zm   = st_data(., sprintf("%s#c.(%s)", this.info.cons2, this.info.zreg), this.info.touse)
    Xf   = J(this.info.N, this.info.Tk-1, 0)
    tX   = J(this.info.N, this.info.Tk-1, 0)
    C    = this.info.cluster == ""? .: st_data(., this.info.cluster, this.info.touse)
    ws   = this.info.wgt? st_data(., this.info.wgtvar, this.info.touse): 1
    wsum = this.info.wgt? quadsum(ws): this.info.N

    for (j = 1; j < this.info.Tk; j++) {
        Xf[., j] = (T :== (j+1))
    }

    this.debug("run linear propensity")
    (void) st_addvar("byte", st_local("Xf"))
    for (j = 1; j < this.info.Tk; j++) {
        (void) st_store(., st_local("Xf"), Xf[., j])
        regcall = sprintf("qui reg %s %s %s", st_local("Xf"), this.info.zreg, regopts)
        stata(regcall)
        tX[., j] = Xf[., j] :- Zm * st_matrix("e(b)")'
    }

    psi_beta = (tX :* ws :* rl) * invsym(quadcross(tX, ws, tX))
    res.Vpop_PL = multe_helper_Vhat(psi_beta, C)
    res.seP[., 1] = sqrt(diagonal(res.Vpop_PL))

    this.debug("run fully interacted model")
    regcall = sprintf("qui reg %s ibn.%s#c.(%s) %s", this.info.Y, this.info.T, this.info.zreg, regopts)
    stata(regcall)
    stata("tempvar ri")
    stata(sprintf("predict %s, resid", st_local("ri")))
    ri = st_data(., st_local("ri"), this.info.touse)

    // TODO: If there are collinear variables within treatment (i.e. in
    // the interacted model) then their coefficient are not identified and
    // Stata loops forever. However, the choice probabilities ARE idnetified
    // (i.e. the likelihood has a maximum just not a unique maximizer). In
    // this case we check convergence in the choice probabilities. Current
    // implementation is ad hoc and neds work.

    this.debug("parse omit info")
    stata("_ms_omit_info e(b)")
    ir_omit = st_matrix("r(omit)")
    ir_gam  = st_matrix("e(b)")
    ir_nam  = st_matrixcolstripe("e(b)")[., 2]
    ir_ord  = order((strtoreal(ustrregexra(ir_nam, "^(\d+).+", "$1")), (1::length(ir_gam))), (1, 2))
    ir_nam  = colshape(ustrregexra(ir_nam[ir_ord], ".+#", ""), this.info.Zk)'
    ir_omit = colshape(ir_omit[ir_ord], this.info.Zk)'
    ir_diff = 0
    for(j = 2; j <= this.info.Tk; j++) {
        ir_diff = any(ir_nam[., 1]  :!= ir_nam[., j]) | any(ir_omit[., 1] :!= ir_omit[., j])
    }
    ir_gam  = ir_gam[ir_ord]
    ir_gam  = (colshape(ir_gam[(this.info.Zk+1)::cols(ir_gam)], this.info.Zk) :- ir_gam[1::this.info.Zk])'
    ir_omit = ir_omit[., 1]'
    ir_nam  = ir_nam[., 1]'
    Zb      = mean(Zm, ws)
    res.estA[., 3] = (Zb * ir_gam)'

    this.debug("interacted model results")
    ir_X   = st_data(., sprintf("ibn.%s#c.(%s)", this.info.T, this.info.zreg), this.info.touse)[., ir_ord]
    psi_al = (ir_X :* ws :* ri) * invsym(quadcross(ir_X, ws, ir_X))
    ir_X   = .
    ate0   = psi_al[., 1..this.info.Zk] * Zb'
    ate    = J(this.info.N, this.info.Tk-1, .)
    for(j = 1; j < this.info.Tk; j++) {
        ate[., j] = psi_al[., ((j * this.info.Zk)+1)..((j+1)*this.info.Zk)] * Zb'
    }
    ate           = ate :- ate0
    res.Vo_ATE    = multe_helper_Vhat(ate, C)
    res.seO[., 3] = sqrt(diagonal(res.Vo_ATE))
    psi_po        = ate :+ (ws :* ((Zm :- Zb) * ir_gam/wsum))
    res.Vpop_ATE  = multe_helper_Vhat(psi_po, C)
    res.Vdiff_ATE = multe_helper_Vhat(psi_beta :- psi_po, C)
    res.seP[., 3] = sqrt(diagonal(res.Vpop_ATE))
    res.seB[., 3] = sqrt(diagonal(res.Vdiff_ATE))

    psi_ownk = J(rows(psi_beta), cols(psi_beta), .)
    psi_rk   = J(rows(psi_beta), cols(psi_beta), 0)
    psi_ri   = J(rows(psi_beta), cols(psi_beta), 0)
    psi_1    = psi_beta
    this.debug("run one-at-a-time")
    for(j = 1; j < this.info.Tk; j++) {
        rk     = multe_helper_olswr((Xf[.,j] :* Zm), (Xf, Zm), ws)
        deltak = rk.coefficients[j,.]
        gamk   = ir_gam[.,j]
        res.estA[j, 2] = deltak * gamk
        if ( this.info.Tk > 2 ) {
            dtXr = multe_helper_olswr(tX[.,j], tX[., multe_helper_antiselect(1..(this.info.Tk-1), j)], ws)
            dtX  = dtXr.residuals
            dtXr = multe_helper_results()
        }
        else {
            dtX = tX
        }
        pp = psi_al[., ((j * this.info.Zk)+1)..((j + 1)*this.info.Zk)] :- psi_al[., 1..this.info.Zk]
        psi_ownk[., j] = pp * deltak' + (ws :* dtX :* rk.residuals :/ quadsum(ws:*(dtX:^2))) * gamk

        // 1 at a time
        s     = selectindex(X0 :| Xf[., j])
        Z_s   = Zm[s, .]
        X_s   = Xf[s, j]
        w_s   = this.info.wgt? ws[s]: 1
        Y_s   = Y[s]
        Xhat  = multe_helper_olswr(X_s, Z_s, w_s)
        rk    = multe_helper_olswr(Y_s, (Xhat.residuals, Z_s), w_s)

        res.estA[j, 4] = rk.coefficients[1, 1]
        psi_k          = w_s :* Xhat.residuals :/ quadsum(w_s:*(Xhat.residuals:^2))
        psi_rk[s, j]   = rk.residuals :* psi_k
        psi_ri[s, j]   = ri[s,.] :* psi_k

        // Test against 1 at a time
        psi_1[s, j] = psi_1[s, j] :- (rk.residuals:*psi_k)
    }
    res.Vpop_OWN  = multe_helper_Vhat(psi_ownk, C)
    res.Vdiff_OWN = multe_helper_Vhat(psi_beta :- psi_ownk, C)
    res.Vpop_EW   = multe_helper_Vhat(psi_rk, C)
    res.Vdiff_EW  = multe_helper_Vhat(psi_1,  C)
    res.Vo_EW     = multe_helper_Vhat(psi_ri, C)

    res.seP[., 2] = sqrt(diagonal(res.Vpop_OWN))
    res.seB[., 2] = sqrt(diagonal(res.Vdiff_OWN))
    res.seB[., 4] = sqrt(diagonal(res.Vpop_EW))
    res.seP[., 4] = sqrt(diagonal(res.Vdiff_EW))
    res.seO[., 4] = sqrt(diagonal(res.Vo_EW))

    this.debug("generalized overlap weights")
    mlfit = J(1, this.info.Tk, "")
    for(j = 1; j <= this.info.Tk; j++) {
        mlfit[j] = st_tempname()
    }
    regcall = sprintf("qui mlogit %s %s %s",
                      this.info.T,
                      this.info.zreg,
                      regopts + " base(1) " + (ir_diff? " iter(100) ": ""))
    stata(regcall)
    stata(sprintf("predict %s, pr", invtokens(mlfit)))
    pis = st_data(., mlfit, this.info.touse)
    regtol = max(pis) * this.info.zerotol
    _edittozero(pis, regtol)
    if ( ir_diff ) {
        ir_start = st_tempname()
        st_matrix(ir_start, st_matrix("e(b)"))
        stata(sprintf(regcall + " from(%s,  skip)", ir_start))
        stata(sprintf("drop %s", invtokens(mlfit)))
        stata(sprintf("predict %s, pr", invtokens(mlfit)))
        mtol = max(reldif(pis, st_data(., mlfit, this.info.touse)))
        if ( mtol > regtol ) {
            pis = J(rows(pis), cols(pis), .)
            errprintf("unable to compute generalized overlap weights\n")
        }
        else {
            pis = st_data(., mlfit, this.info.touse)
            ir_diff = 0
        }
        _edittozerotol(pis, regtol)
    }

    // xx TODO: check whether omitting is necessary; invsym is a generalized inverse
    // stata("_ms_omit_info e(b)")
    // ir_omit = st_matrix("r(omit)")
    // ir_gam  = st_matrix("e(b)")
    // ir_nam  = st_matrixcolstripe("e(b)")[., 2]
    // ir_omit = colshape(ir_omit, this.info.Zk)'
    // ir_nam  = colshape(ustrregexra(ir_nam, ".+#", ""), this.info.Zk)'
    // for(j = 2; j <= this.info.Tk; j++) {
    //     assert(all(ir_nam[., 2]  :== ir_nam[., j]))
    //     assert(all(ir_omit[., 2] :== ir_omit[., j]))
    // }
    // ir_gam  = (colshape(ir_gam, this.info.Zk) :- ir_gam[1::this.info.Zk])'
    // ir_omit = ir_omit[., 2]'
    // ir_nam  = ir_nam[., 2]'
    // Zm      = select(st_data(., invtokens(ir_nam), this.info.touse), !ir_omit)
    // Zm      = select(Zm, !ir_omit)
    // this.info.Zk = cols(Zm)

    Sc = J(this.info.N, (this.info.Tk-1)*this.info.Zk, .)
    He = J((this.info.Tk-1) * this.info.Zk, (this.info.Tk-1) * this.info.Zk, .)
    for (i = 1; i < this.info.Tk; i++) {
        seli = (this.info.Zk * (i-1) + 1)..(this.info.Zk*i)
        Sc[., seli] = ws :* (Xf[.,i] :- pis[., i+1]) :* Zm
        for (j = 1; j < this.info.Tk; j++) {
            selj = (this.info.Zk * (j-1) + 1)..(this.info.Zk*j)
            He[seli', selj] = quadcross(Zm, ws :* pis[., i+1] :* ((i == j) :- pis[, j+1]), Zm)
        }
    }

    this.debug("Wald and LM tests (assumes first column is intercept)")
    idx1   = (0..(this.info.Tk-2)):*this.info.Zk:+1
    idx1n  = multe_helper_antiselect(1..((this.info.Tk-1) * this.info.Zk), idx1)
    He1112 = invsym(He[idx1, idx1]) * He[idx1, idx1n]
    Vu     = multe_helper_Vhat(Sc[., idx1n] :- Sc[., idx1] * He1112, C)
    // th1    = vec(ir_gam[selectindex(!ir_omit), 2..this.info.Tk])
    th1    = st_matrix("e(b)")[(this.info.Zk+1)..(this.info.Zk*this.info.Tk)]
    th     = (He[idx1n, idx1n] :- He[idx1n, idx1] * He1112) * th1[idx1n]'

    // LM, calculate weighted colMeans(Xf)
    pis0  = mean(Xf, ws)
    Scr   = J(rows(Zm), (this.info.Tk-1)*this.info.Zk, .)
    Her   = J((this.info.Tk-1) * this.info.Zk, (this.info.Tk-1) * this.info.Zk, .)
    for (i = 1; i < this.info.Tk; i++) {
        seli = (this.info.Zk * (i-1) + 1)..(this.info.Zk*i)
        Scr[., seli] = ws :* (Xf[.,i] :- pis0[., i]) :* Zm
        for (j = 1; j < this.info.Tk; j++) {
            selj = (this.info.Zk * (j-1) + 1)..(this.info.Zk*j)
            Her[seli', selj] = quadcross(Zm, ws :* pis0[., i] :* ((i == j) :- pis0[, j]), Zm)
        }
    }

    Her1112 = invsym(Her[idx1, idx1]) * Her[idx1, idx1n]
    Vr = multe_helper_Vhat(Scr[., idx1n] :- Scr[., idx1] * Her1112, C)
    LM = multe_helper_quadform(Vr, quadcolsum(Scr[., idx1n])')
    Wa = ir_diff? (., .): multe_helper_quadform(Vu, th)
    res.LM.stat = LM[1]
    res.Wa.stat = Wa[1]
    res.LM.df   = LM[2]
    res.Wa.df   = Wa[2]
    res.LM.pval = chi2tail(LM[2], LM[1])
    res.Wa.pval = chi2tail(Wa[2], Wa[1])

    this.debug("generalized overlap weights + standard errors")
    pis0 = mean(X0, ws), pis0
    if ( !this.info.cw_uniform ) {
        vpi = pis0 :* (1 :- pis0)
    }
    else {
        vpi = 1
    }

    // Note ipi is always multiplied by lam; if ipi is set to 1/lam then it
    // will give a missing when lam is 0, but missing times 0 is missing and
    // we want it to be 0, so we set missing ipi values to 0.
    lam = editmissing(1 :/ quadrowsum(vpi :/ pis), 0)
    _edittozerotol(lam, max(lam) * this.info.zerotol)
    ipi = editmissing(1:/pis, 0)
    cw  = lam :/ quadrowsum((X0, Xf) :* pis)

    this.debug("efficient common weights")
    if ( quadsum(lam) > 0 ) {
        ro             = multe_helper_olswr(Y, (J(rows(Xf), 1, 1), Xf), cw :* ws)
        res.estA[., 5] = ro.coefficients[2::this.info.Tk]
        psi_or         = (ws :* lam :/ quadsum(ws :* lam)) :* ((Xf :* ipi[., 2..this.info.Tk]) :- (X0 :* ipi[, 1]))
        res.Vo_CW      = multe_helper_Vhat(ri :* psi_or, C)
        res.seO[., 5]  = sqrt(diagonal(res.Vo_CW))

        // Final step: calculate se_po for overlap weights. model.matrix
        // rather than ml$residuals to account for rounding
        xf1 = lam :* (vpi[., 2..this.info.Tk] :* ipi[., 2..this.info.Tk]) :* X0
        M0  = J(this.info.Zk, this.info.Tk-1, .)
        for (i = 1; i < this.info.Tk; i++) {
            M0[., i] = quadcolsum((xf1[., i] :* ws :* cw :* ro.residuals) :* Zm)'
        }
        M    = J(this.info.Zk * (this.info.Tk-1), this.info.Tk-1, .)
        si   = J(1, (this.info.Tk-2), 0), 1, J(1, (this.info.Tk-2), 0)
        xf1l = lam :* (vpi[., 2..this.info.Tk] :* ipi[., 2..this.info.Tk])
        for (j = 1; j < this.info.Tk; j++) {
            xf1r = J(rows(Xf), 1, si[(this.info.Tk-j)..(2*this.info.Tk-2-j)])
            xf1  = (xf1l :- xf1r) :* Xf[., j]
            for (i = 1; i < this.info.Tk; i++) {
                seli = (this.info.Zk * (i-1) + 1)::(this.info.Zk*i)
                M[seli, j] = quadcolsum((xf1[., i] :* ws :* cw :* ro.residuals) :* Zm)'
            }
        }

        a = (Sc * invsym(He)) * (M :- vec(M0)) :/ quadsum(ws:*lam)
        res.Vpop_CW   = multe_helper_Vhat(ro.residuals :* psi_or :+ a, C)
        res.Vdiff_CW  = multe_helper_Vhat(psi_beta :- ro.residuals :* psi_or :- a, C)
        res.seP[., 5] = sqrt(diagonal(res.Vpop_CW))
        res.seB[., 5] = sqrt(diagonal(res.Vdiff_CW))
    }
    else {
        errprintf("Sample for efficient common weights is empty.\n")
    }

    res.estB[., 2..5] = res.estA[., 1] :- res.estA[., 2..5]
    res.estB = res.estB[., 2..5]
    res.seB  = res.seB[., 2..5]

    return(res)
}

void function MulTE::debug(string scalar msg)
{
    if ( this.info.debug ) printf("debug %g: %s\n", this.info.debug, msg)
    else return
    (void) this.info.debug++
}

void function multe_helper_post(string scalar b, string scalar V, struct MulTE_Return scalar res)
{
    real scalar uselevels
    real matrix best, sest
    string vector rownames, eqnames

    uselevels = any(st_vartype(res.Tvar) :== ("byte", "int", "long"))
    // uselevels = uselevels | all(floor(Tvalues) :== Tvalues)
    rownames  = colshape(res.Tlevels, 1)
    if ( uselevels ) {
        rownames[1] = rownames[1] + "b"
        rownames    = colshape(rownames :+ "." :+ res.Tvar, 1)
    }
    else if ( !strpos(st_vartype(res.Tvar), "str") ) {
        printf("{bf:note:} table coefficients correspond to treatment values\n")
    }
    best    = 0, rowshape(st_matrix(b), 1)
    eqnames = J(1, res.Tk, "")'
    sest    = J(res.Tk, res.Tk, 0)
    sest[2::res.Tk, 2..res.Tk] = st_matrix(V)

    st_matrix(b, rowshape(best', 1))
    st_matrixcolstripe(b, (eqnames, rownames))
    st_matrixrowstripe(b, ("", res.Yvar))

    st_matrix(V, sest)
    st_matrixcolstripe(V, (eqnames, rownames))
    st_matrixrowstripe(V, (eqnames, rownames))
}
end

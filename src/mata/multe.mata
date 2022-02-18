cap mata mata drop MulTE()
cap mata mata drop MulTE_Results()
cap mata mata drop multe_helper_ols()
cap mata mata drop multe_helper_olsw()

mata
class MulTE_Results
{
    real matrix est
    real matrix se_po
    real matrix se_or
    real   vector Tvalues
    string vector Tlabels
    string vector colnames

    void print()
    void export()
}

class MulTE_Results scalar MulTE(string scalar Yvar, string scalar Tvar, real matrix Wm, string scalar touse)
{
    class MulTE_Results scalar results

    string scalar Tlab

    real scalar j, n, k
    real vector Y, X
    real vector xlevels, X0, alpha0, lam, rcalpha, rcres, rkalpha, rkres, ts, s0, Wmean
    real matrix Xm, psi_alpha0, ps, Xt
    real matrix est, se_or, se_po

    Y       = st_data(., Yvar, touse)
    X       = st_data(., Tvar, touse)
    n       = rows(Y)
    xlevels = uniqrows(X)
    k       = length(xlevels)
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

* TODO: xx figure out if selectindex(s) is faster
* TODO: xx or maybe make subset vars at start of loop
    // loop xx
    Wmean = colsum(Wm) / n
    for(j = 2; j <= k; j++) {
        alphak  = multe_helper_ols(select(Y, Xm[., j]), select(Wm, Xm[., j]))
        psi_alphak = (Xm[.,j] :* (Y - Wm * alphak) :* Wm) * qrinv(cross(Xm[., j] :* Wm, Wm) / n)

        // ATE
        est[j - 1, 1] = (Wmean * (alphak - alpha0)) / n
        psi_or = ((psi_alphak - psi_alpha0) * Wmean') / n
        se_or[j - 1, 1] = sqrt(variance(psi_or) * (n - 1) / n^2)
        se_po[j - 1, 1] = sqrt(variance(psi_or + ((Wm :-  Wmean) * (alphak - alpha0))) * (n - 1) / n^2)

        // One treatment at a time

        s               = (X0 :| Xm[., j])
        Xdot            = select(Xm[., j], s) - select(Wm, s) * multe_helper_ols(select(Xm[., j], s), select(Wm, s))
        rkalpha         = multe_helper_ols(select(Y, s), (Xdot, select(Wm, s)))
        est[j - 1, 2]   = rkalpha[1]
        rkres           = select(Y, s) - (Xdot, select(Wm, s)) * rkalpha
        se_po[j - 1, 2] = sqrt(sum((rkres:^2) :* (Xdot:^2)) / sum(Xdot:^2):^2)

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
    results.est   = est
    results.se_po = se_po
    results.se_or = se_or
    results.Tvalues = xlevels
    results.Tlabels = Tlab == ""? strofreal(xlevels): st_vlmap(Tlab, xlevels)
}

real matrix function multe_helper_ols(real matrix Y, real matrix X)
{
    return(invsym(cross(X, X)) * cross(X, Y))
}

real matrix function multe_helper_olsw(real matrix Y, real matrix X, real matrix W)
{
    return(qrinv(cross(X :* W, X)) * cross(X :* W, Y))
}

void function MulTE_Results::init()
{
    colnames = ("ATE", "1-at-a-time", "common weights")'
}

void function MulTE_Results::export(string scalar outmatrix)
{
    string vector rownames
    rownames = Tlabels[2::length(Tvalues)], J(rows(est), 1, "se"), J(rows(est), 1, "oracle_se")
    rownames = rowshape(rownames, rows(est) * 3)

    st_matrix(outmatrix, rowshape((est, se_po, se_or), rows(est) * 3))
    st_matrixcolstripe(outmatrix, (J(3, 1, ""), colnames))
    st_matrixrowstripe(outmatrix, (J(3 * rows(est), 1, ""), rownames))
}

void function MulTE_Results::print(| real scalar digits)
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

void function xx(formats, colnames)
{
    printf("|")
    printf(formats[1], "")
    for (j = 1; j <= length(colnames); j++) {
        printf("|")
        printf(formats[j + 1], colnames[j])
    }
    printf("|\n")
}

end

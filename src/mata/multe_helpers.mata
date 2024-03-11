cap mata mata drop multe_helper_ols()
cap mata mata drop multe_helper_olsw()
cap mata mata drop multe_helper_olsr()
cap mata mata drop multe_helper_olswr()
cap mata mata drop multe_helper_olse()
cap mata mata drop multe_helper_olswe()
cap mata mata drop multe_helper_results()
cap mata mata drop multe_helper_antiselect()
cap mata mata drop multe_helper_licols()
cap mata mata drop multe_helper_Vhat()
cap mata mata drop multe_helper_sehat()
cap mata mata drop multe_helper_quadform()

// without clustering multe in R doesn't do the small-sample adjustment
mata:
real matrix function multe_helper_Vhat(real matrix psi, real colvector cluster)
{
    real scalar G
    real matrix _psi, psisum
    _psi = editmissing(psi, 0)
    if ( missing(cluster) ) {
        // G = rows(_psi)
        // return((G / (G - 1)) :* quadcross(_psi, _psi))
        return(quadcross(_psi, _psi))
    }
    else {
        psisum = sort((cluster, _psi), 1)
        psisum = panelsum(psisum[., 2..cols(psisum)], panelsetup(psisum, 1))
        G = rows(psisum)
        return((G / (G - 1)) :* quadcross(psisum, psisum))
    }
}

real vector function multe_helper_sehat(real matrix psi, real colvector cluster)
{
    return(sqrt(diagonal(multe_helper_Vhat(psi, cluster))))
}

struct multe_helper_results {
    real matrix coefficients
    real matrix residuals
    real matrix vcov
}

// invsym _should_ work under collinearity; can't recall the example where it failed
real matrix function multe_helper_ols(real matrix Y, real matrix X)
{
    return(invsym(quadcross(X, X), 1..cols(X)) * quadcross(X, Y))
}

// return(qrsolve(quadcross(X, W, X), quadcross(X, W, Y)))
real matrix function multe_helper_olsw(real matrix Y, real matrix X, real matrix W)
{
    return(invsym(quadcross(X, W, X), 1..cols(X)) * quadcross(X, W, Y))
}

struct multe_helper_results scalar function multe_helper_olsr(real matrix Y, real matrix X)
{
    struct multe_helper_results scalar results
    results.coefficients = multe_helper_ols(Y, X)
    results.residuals    = Y - X * results.coefficients
    return(results)
}

struct multe_helper_results scalar function multe_helper_olswr(real matrix Y, real matrix X, real matrix W)
{
    struct multe_helper_results scalar results
    results.coefficients = multe_helper_olsw(Y, X, W)
    results.residuals    = Y - X * results.coefficients
    return(results)
}

struct multe_helper_results scalar function multe_helper_olse(real matrix Y, real matrix X)
{
    struct multe_helper_results scalar results
    real matrix XX, XDX
    results      = multe_helper_olsr(Y, X)
    XX           = invsym(quadcross(X, X), 1..cols(X))
    XDX          = quadcross(X, results.residuals:^2, X)
    results.vcov = XX * XDX * XX
    return(results)
}

struct multe_helper_results scalar function multe_helper_olswe(real matrix Y, real matrix X, real matrix W)
{
    struct multe_helper_results scalar results
    real matrix XX, XDX
    results      = multe_helper_olswr(Y, X, W)
    XX           = invsym(quadcross(X, W, X), 1..cols(X))
    XDX          = quadcross(X, (W :* results.residuals):^2, X)
    results.vcov = XX * XDX * XX
    return(results)
}

real matrix function multe_helper_antiselect(real matrix x, real vector ix)
{
    real vector sel
    if ( length(ix) == 1 ) {
        if ( min((rows(x), cols(x))) > 1 ) {
            errprintf("multe_helper_antiselect(): ambiguous scalar anti-selection of a matrix")
            _error(198)
        }
        sel = rows(x) > cols(x)? J(rows(x), 1, 1): J(1, cols(x), 1)
    }
    else {
        sel = rows(ix) > cols(ix)? J(rows(x), 1, 1): J(1, cols(x), 1)
    }
    sel[ix] = J(rows(ix) > cols(ix)? length(ix): 1, rows(ix) > cols(ix)? 1: length(ix), 0)
    return(select(x, sel))
}

real rowvector function multe_helper_licols(real matrix X)
{
    if ( (cols(X) == 0) | (rows(X) == 0) ) return(J(1, 0, 0))
    return(rowshape(diagonal(invsym(quadcross(X, X), 1..cols(X))) :!= 0, 1))
}

real vector function multe_helper_quadform(real matrix A, real vector b)
{
    real scalar rank, quadf
    real matrix Ainv
    if ( !issymmetric(A) ) {
        errprintf("multe_helper_quadform(): matrix must be symmetric")
        _error(198)
    }
    Ainv = invsym(A)
    rank = sum(diagonal(Ainv) :!= 0)
    if ( rows(b) > cols(b) ) {
        quadf = b' * Ainv * b
    }
    else {
        quadf = b * Ainv * b'
    }
    return((quadf, rank))
}
end

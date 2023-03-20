cap mata mata drop multe_helper_ols()
cap mata mata drop multe_helper_olsw()
cap mata mata drop multe_helper_olsr()
cap mata mata drop multe_helper_olswr()
cap mata mata drop multe_helper_olse()
cap mata mata drop multe_helper_olswe()
cap mata mata drop multe_helper_results()
cap mata mata drop multe_helper_antiselect()

mata
struct multe_helper_results {
    real matrix coefficients
    real matrix residuals
    real matrix vcov
}

// using qrsolve directly can at times be more precise and less precise...
// return(invsym(quadcross(X, X)) * quadcross(X, Y))
// return(qrsolve(X, Y))
real matrix function multe_helper_ols(real matrix Y, real matrix X)
{
    return(qrsolve(quadcross(X, X), quadcross(X, Y)))
}

// return(invsym(quadcross(X, W, X)) * quadcross(X, W, Y))
// return(min(W) < 0? invsym(quadcross(X, W, X)) * quadcross(X, W, Y): qrsolve(sqrt(W) :* X, sqrt(W) :* Y))
real matrix function multe_helper_olsw(real matrix Y, real matrix X, real matrix W)
{
    return(qrsolve(quadcross(X, W, X), quadcross(X, W, Y)))
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
    XX           = invsym(quadcross(X, X))
    XDX          = quadcross(X, results.residuals:^2, X)
    results.vcov = XX * XDX * XX
    return(results)
}

struct multe_helper_results scalar function multe_helper_olswe(real matrix Y, real matrix X, real matrix W)
{
    struct multe_helper_results scalar results
    real matrix XX, XDX
    results      = multe_helper_olswr(Y, X, W)
    XX           = invsym(quadcross(X, W, X))
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
    sel[ix] = J(length(ix), 1, 0)
    return(select(x, sel))
}
end

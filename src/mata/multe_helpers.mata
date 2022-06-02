cap mata mata drop multe_helper_ols()
cap mata mata drop multe_helper_olsw()
cap mata mata drop multe_helper_olsr()
cap mata mata drop multe_helper_olswr()
cap mata mata drop multe_helper_results()
cap mata mata drop multe_helper_antiselect()

mata
struct multe_helper_results {
    real matrix coefficients
    real matrix residuals
}

real matrix function multe_helper_ols(real matrix Y, real matrix X)
{
    return(invsym(cross(X, X)) * cross(X, Y))
}

real matrix function multe_helper_olsw(real matrix Y, real matrix X, real matrix W)
{
    return(qrinv(cross(X :* W, X)) * cross(X :* W, Y))
}

struct multe_helper_results scalar function multe_helper_olsr(real matrix Y, real matrix X)
{
    struct multe_helper_results scalar results
    results.coefficients = invsym(cross(X, X)) * cross(X, Y)
    results.residuals    = Y - X * results.coefficients
    return(results)
}

struct multe_helper_results scalar function multe_helper_olswr(real matrix Y, real matrix X, real matrix W)
{
    struct multe_helper_results scalar results
    results.coefficients = qrinv(cross(X :* W, X)) * cross(X :* W, Y)
    results.residuals    = Y - X * results.coefficients
    return(results)
}

real matrix function multe_helper_antiselect(real matrix x, real vector ix)
{
    real vector sel
    sel = rows(ix) > cols(ix)? J(rows(x), 1, 1): J(1, cols(x), 1)
    sel[ix] = J(length(ix), 1, 0)
    return(select(x, sel))
}
end

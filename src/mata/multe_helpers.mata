cap mata mata drop multe_helper_ols()
cap mata mata drop multe_helper_olsw()
cap mata mata drop multe_helper_olsr()
cap mata mata drop multe_helper_results()

mata
end

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
end

cap mata mata drop multe_helper_ols()
cap mata mata drop multe_helper_olsw()
cap mata mata drop multe_helper_olsr()
cap mata mata drop multe_helper_olswr()
cap mata mata drop multe_helper_results()
cap mata mata drop multe_helper_antiselect()
cap mata mata drop MulTE_RunningTimer()
cap mata mata drop MulTE_LoopTimer()

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

// ----------------------------------------------------------------- //
//                               TIMER                               //
// ----------------------------------------------------------------- //

mata
class MulTE_LoopTimer
{
    real scalar ntimers
    real scalar started
    real scalar timer_total
    real scalar timer_current
    real rowvector timer_step
    real rowvector timer_times
    string rowvector timer_msg

    void new()
    void start()
    void finish()
    void time()
    real scalar free_timer()
}

void function MulTE_LoopTimer::new()
{
    started = 0
}

void function MulTE_LoopTimer::start(real scalar _ntimers)
{
    real scalar i
    started = 1
    ntimers = _ntimers
    timer_step  = J(1, ntimers, .)
    timer_times = J(1, ntimers, 0)
    timer_msg   = J(1, ntimers, "")
    for(i = 1; i <= ntimers; i++) {
        timer_step[i] = free_timer()
        timer_off(timer_step[i])
        timer_clear(timer_step[i])
        timer_on(timer_step[i])
    }

    timer_total = free_timer()
    timer_off(timer_total)
    timer_clear(timer_total)
    timer_on(timer_total)

    for(i = 1; i <= ntimers; i++) {
        timer_off(timer_step[i])
    }

    timer_current = 1
    timer_on(timer_step[timer_current])
}

real scalar function MulTE_LoopTimer::free_timer()
{
    real scalar i
    i = 99
    while ( (i > 0) & timer_value(i)[2] ) {
        i = i - 1
    }
    return(i)
}

void function MulTE_LoopTimer::time(| string scalar msg)
{

    if ( args() == 0 ) msg = ""

    timer_off(timer_step[timer_current])
    timer_msg[timer_current] = msg
    timer_times[timer_current] = timer_times[timer_current] + timer_value(timer_step[timer_current])[1]
    timer_clear(timer_step[timer_current])
    timer_current = timer_current < ntimers? timer_current + 1: 1

    timer_off(timer_step[timer_current])
    timer_clear(timer_step[timer_current])
    timer_on(timer_step[timer_current])
}

void function MulTE_LoopTimer::finish(| string scalar msg)
{
    real scalar i
    string scalar s, s_total, s_step
    if ( started == 0 ) return

    for(i = 1; i <= ntimers; i++) {
        timer_off(timer_step[i])
        timer_clear(timer_step[i])
    }
    timer_off(timer_total)
    s_total = strtrim(sprintf("%21.1fc", timer_value(timer_total)[1]))
    timer_clear(timer_total)

    if ( args() == 0 ) msg = "Loop Timer"

    printf("\t%s (%ss)\n", msg, s_total)
    for (i = 1; i <= ntimers; i++) {
        s = (strlen(timer_msg[i]) & substr(timer_msg[i], strlen(timer_msg[i]), 1) != " ")? " ": ""
        s_step = strtrim(sprintf("%21.2fc", timer_times[i]))
        printf("\t\t%s%s(%ss)\n", timer_msg[i], s, s_step)
    }
}

// ----------------------------------------------------------------- //
//                               TIMER                               //
// ----------------------------------------------------------------- //


class MulTE_RunningTimer
{
    real scalar started
    real scalar timer_step
    real scalar timer_total

    void new()
    void start()
    void finish()
    void time()
    real scalar free_timer()
}

void function MulTE_RunningTimer::new()
{
    started = 0
}

void function MulTE_RunningTimer::start()
{
    started = 1
    timer_step = free_timer()
    timer_off(timer_step)
    timer_clear(timer_step)
    timer_on(timer_step)

    timer_total = free_timer()
    timer_off(timer_total)
    timer_clear(timer_total)
    timer_on(timer_total)
}

void function MulTE_RunningTimer::finish()
{
    if ( started ) {
        timer_off(timer_step)
        timer_clear(timer_step)

        timer_off(timer_total)
        timer_clear(timer_total)
    }
}

real scalar function MulTE_RunningTimer::free_timer()
{
    real scalar i
    i = 99
    while ( (i > 0) & timer_value(i)[2] ) {
        i = i - 1
    }
    return(i)
}

void function MulTE_RunningTimer::time(| string scalar msg)
{
    if ( args() == 0 ) msg = ""
    if ( started == 0 ) start()

    string scalar s, s_step, s_total
    timer_off(timer_step)
    timer_off(timer_total)

    s_step  = strtrim(sprintf("%21.2fc", timer_value(timer_step)[1]))
    s_total = strtrim(sprintf("%21.1fc", timer_value(timer_total)[1]))

    timer_clear(timer_step)
    timer_on(timer_step)
    timer_on(timer_total)

    s = (strlen(msg) & substr(msg, strlen(msg), 1) != " ")? " ": ""
    printf("\t%s%s(%ss / %ss)\n", msg, s, s_step, s_total)
}
end

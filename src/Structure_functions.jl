module Structure_functions

import StatsBase, CurveFit, LsqFit

export fct_sct
export fit_fctsct
export xhi_fct_p
export xhi_fct_p!

function fct_sct(cvicube,LAG,ORDERS)
    sct = Array{Float64}(undef,size(ORDERS)[1],size(LAG)[1])
    for ord=1:size(ORDERS)[1]    
        for lag=1:size(LAG)[1]
            sct[ord,lag]= StatsBase.mean(skipmissing(abs.(cvicube[:,:,lag])).^ORDERS[ord])
        end
    end
    return(sct)
end

"""
    fit_fctstruct(xdata,ydata,y0,confidinterv)

Fit the model ydata=A*(xdata)^B. Return in first B, then A.
"""
function fit_fctsct(xdata,ydata,y0 ; confinterv=true)
    model(x,xhi) = xhi[2].*x.^xhi[1]
    fit = LsqFit.curve_fit(model, ydata, xdata, y0)
    #confid_interval = LsqFit.standard_error(fit)
    #confinterv==false && return(fit.param[1],fit.param[2])
    return(fit.param[1],fit.param[2],fit.resid[1])
end

#= WORKING
function fit_fctsct(xdata,ydata)
    fitted = CurveFit.curve_fit(CurveFit.PowerFit, ydata, xdata)
    #confid_interval = confidence_interval(fit, confidinterv)
    return(fitted.coefs[2],fitted.coefs[1])
end
=#
"""
    xhi_fct_p(pvec,struct_lag,y0)   

Fit the model y=A*(x)^B for multiple order and lag of structure functions.
"""
function xhi_fct_p(pvec,struct_lag)
    zeta = Array{Float64}(undef,size(pvec)[1],3)
    p0 = [[0.37,1],[0.7,1],[1,1],[1.27,1],[1.53,1],[1.77,1]]
    for ix=1:size(pvec)[1]
        
        zeta[ix,:] .= fit_fctsct(struct_lag[ix,:],struct_lag[3,:],p0[ix])
        # WORKED with other version fit_fctsct
        #zeta[ix,:] .= fit_fctsct(struct_lag[ix,:],struct_lag[3,:])
    end
    return(zeta)
end

function xhi_fct_p!(pvec,struct_lag,lagvec,y0,confidinterv,zeta)
    for ix in eachindex(pvec)
        zeta[ix,:] .= fit_fctsct(struct_lag[ix,:],struct_lag[3,:])
    end
    return(zeta)
end


end #module
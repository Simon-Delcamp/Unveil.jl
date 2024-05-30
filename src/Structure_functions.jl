module Structure_functions

using StatsBase, CurveFit, Plots
export fct_sct
export fit_fctsct
export xhi_fct_p
export xhi_fct_p!


"""
   fct_sct(cvicube,LAG,ORDERS)

Compute the structure functions from Monin & Yaglom+75 method.
"""
function fct_sct(cvicube,LAG,ORDERS)
    sct = Array{Float64}(undef,size(ORDERS)[1],size(LAG)[1])
    for ord=1:size(ORDERS)[1]    
        for lag=1:size(LAG)[1]
            #sct[ord,lag]= mean(skipmissing(abs.(cvicube[:,:,lag])).^ORDERS[ord])
            sct[ord,lag]= mean(skipmissing(abs.(cvicube[:,lag])).^ORDERS[ord])
        end
    end
    return(sct)
end


"""
   fct_sct_int(cvicube,LAG,ORDERS)

Compute the structure functions from Hily-Blant+08 method.
"""
function fct_sct_int(cvicube,LAG,ORDERS)
    spl = Array{Float64}(undef,size(ORDERS)[1],size(LAG)[1])
    for lx=1:size(LAG)[1]
        cvire = reshape(cvicube[:,:,lx],size(cvicube)[1]*size(cvicube)[2])

        # Histogram
        hist = fit(Histogram,cvire,-30:0.01:30)
        deltax = abs(hist.edges[1][2]-hist.edges[1][1])
        xhist = hist.edges[1][1]+abs(hist.edges[1][2]-hist.edges[1][1])/2:abs(hist.edges[1][2]-hist.edges[1][1]):hist.edges[1][end]

        # Normalisation of the histogram (mean 0 and dispersion unity)
        sumi = sum(hist.weights)
        temp = hist.weights/(sumi*deltax)
        tx = xhist
        ty = temp
        model(x,xhi) = exp.(.-(x.-xhi[1]).^2 ./2 ./xhi[2].^2)./sqrt.(2pi)./xhi[2]
        fitted = curve_fit(model, tx, ty, [0.0,1.0])
        tx = (xhist.-fitted.param[1])./fitted.param[2]
        ty = temp.*fitted.param[2]
        deltatx = abs(tx[1]-tx[2])

        #Sp(l) computation
        s(p) = sum(abs.(xhist).^(p).*ty.*deltatx)
        for po = 1:size(ORDERS)[1]
            spl[po,lx] = s(po)
        end
    end
    return(spl)
end

"""
    fit_fctstruct(xdata,ydata,y0,confidinterv)

Fit the model ydata=A*(xdata)^B. Return in first B, then A.
"""
function fit_fctsct(xdata,ydata,y0 ; confinterv=true)
    model(x,xhi) = xhi[2].*x.^xhi[1]
    fit = curve_fit(model, ydata, xdata, y0)
    #confid_interval = standard_error(fit)
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
    pp = [(i-1)*0.5+0.37 for i=1:size(pvec)[1]] 
    pf = [1 for i=1:size(pvec)[1]]
    p0 = [[pp[ix],pf[ix]] for ix=1:size(pvec)[1]] # Guesses for exponents B and factor A. 
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
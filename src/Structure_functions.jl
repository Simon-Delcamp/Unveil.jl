###################################################################
# STILL WORK IN PROGRESS
# Functions for computing and using the Structures Functions. Though for use after CVI computing.
# Called these functions by calling (or writting on a script ):
#       >include("../src/Structure_fct.jl")
#       >using .Structure_fct
#       >output = Structure_fct.NameOfTheFunction(input)
###################################################################

module Structure_functions

include("Functionforcvi.jl") #Calculations of CVI
include("Data_preparation.jl") #Read and write fits

using .Functionforcvi
using .Data_preparation

using Plots, StatsPlots,MultivariateStats, Statistics, StatsBase, Distributions, LinearAlgebra
using KernelDensity,CurveFit#,LsqFit

export construct_fctstruct!
export construct_fctstruct
export fct_struct
export fit_fctstruct
export pdf_normed
export xhi_fct_p




function fct_sct(cvicube,LAG,ORDERS)
    sct = Array{Float64}(undef,size(ORDERS)[1],size(LAG)[1])
    for ord=1:size(ORDERS)[1]    
        for lag=1:size(LAG)[1]
            sct[ord,lag]=mean(skipmissing(abs.(cvicube[:,:,lag])).^ORDERS[ord])
        end
    end
    return(sct)
end








"""
    construct_fctstruct!

Construct a matrix with the values of the structure functions at every given lag and orders (lag in columns and orders in lines). Need a preallocated array.
"""
function construct_fctstruct!(cube,struct_fct_order,Lag,bins,cvilim::Tuple{Int64,Int64},blank,struct_lag)
   for ix in eachindex(Lag)
    tostruct          = Data_preparation.delete_allnotvalue(cube[:,:,ix],blank)
    pdf_cube          = Structure_fct.pdf_normed(tostruct,bins,cvilim,blank)
     for jx in eachindex(struct_fct_order)
            # IDEE : IF pdf_cube.density<0, put them =0
            struct_lag[jx,ix] = Structure_fct.fct_struct(abs.(tostruct),pdf_cube,struct_fct_order[jx])
            #struct_lag[jx,ix] = Structure_fct.fct_struct(abs.(tostruct),pdf_cube,struct_fct_order[jx],bins)
        end
    end
    return(struct_lag)
end

function OLDconstruct_fctstruct!(cube,struct_fct_order,Lag,bins,cvilim::Tuple{Int64,Int64},blank,struct_lag,rep)
    @views @inbounds for ix in eachindex(Lag)
        @views @inbounds for jx in eachindex(struct_fct_order)
            tostruct          = Data_preparation.delete_allnotvalue(cube[:,:,ix],blank)
            pdf_cube          = Structure_fct.pdf_normed(tostruct,bins,cvilim,blank)
            # IDEE : IF pdf_cube.density<0, put them =0
            rep = zeros(Float64,size(pdf_cube.x)[1])
            xar = Array{Float64}(pdf_cube.x[1]:dx:pdf_cube.x[end])
            struct_lag[jx,ix] = Structure_fct.fct_struct!(abs.(tostruct),pdf_cube,struct_fct_order[jx],rep,xar)
        end
    end
    return(struct_lag)
end




"""
    construct_fctstruct

Construct a matrix with the values of the structure functions at every given lag and orders (lag in lines and orders in columns).
"""
function construct_fctstruct(cube,struct_fct_order,Lag,bins,cvilim::Tuple{Int64,Int64},blank)
    struct_lag = Array{Float64}(undef,size(struct_fct_order)[1],size(Lag)[1])
    tostruct = Data_preparation.delete_allnotvalue(cube[:,:,2],blank)
    pdf_cube = kde(tostruct,npoints=bins)
    pdf_cube.density = Data_preparation.blank_inf(pdf_cube.density,0.0,blank)
    pdf_cube.density = Data_preparation.delete_allnotvalue(pdf_cube.density,blank)
    dx = abs(pdf_cube.x[1]-pdf_cube.x[2])
    xar = Array{Float64}(pdf_cube.x[1]:dx:pdf_cube.x[end])
    rep = zeros(Float64,size(pdf_cube.x)[1])
    for ix=1:size(Lag)[1]
        for jx=1:size(struct_fct_order)[1]
            tostruct = Data_preparation.delete_allnotvalue(cube[:,:,ix],blank)
            # Construct the cvi map at lag l on densities 
            pdf_cube = kde(tostruct,npoints=bins)
            #pdf_cube,sig,mu          = Structure_fct.pdf_normed(tostruct,bins,cvilim,blank)
            #mu       = moment(pdf_cube.x,1,aweights(pdf_cube.density),0)
            #sig      = sqrt(moment(pdf_cube.x,2,aweights(pdf_cube.density)))
            #pdf_cube.x = (pdf_cube.x.*sig).+mu

            pdf_cube.density = Data_preparation.blank_inf(pdf_cube.density,0.0,blank)
            pdf_cube.density = Data_preparation.delete_allnotvalue(pdf_cube.density,blank)
            struct_lag[jx,ix] = Structure_fct.fct_struct(tostruct,pdf_cube,struct_fct_order[jx],dx,xar,rep)
            #struct_lag[jx,ix] = Structure_fct.fct_struct(abs.(pdf_cube.density),pdf_cube,struct_fct_order[jx])
        end
    end
    return(struct_lag)
end





"""
    fct_struct(data,datapdf,order,bin)

Compute the structure function of a data at a specific order.
"""
function fct_struct(maplag,mappdf, order,dx,xar,rep)
    #init_step = findall(x->x>=0, xar)[1]
    for ix=1:size(xar)[1]-1
        for jx=1:size(maplag)[1]
            xar[ix+1]>=abs(maplag[jx]) && xar[ix]<=abs(maplag[jx]) && (rep[ix]+=mappdf.density[ix].*abs(maplag[jx]).^order)
        end
        #(rep[ix]+=mappdf.density[ix]*abs.(maplag).^order)
        #rep[ix] = rep[ix]*dx
    end

    #dcvi = abs(maximum(maplag)-minimum(maplag))/size(maplag)[1] 
    #return(sum(rep)*dcvi)
    return(sum(rep)*dx)
end



function OLDfct_struct(data,datapdf, order)
    dx = abs(datapdf.x[1]-datapdf.x[2])
    xar = Array{Float64}(datapdf.x[1]:dx:datapdf.x[end])
    rep = zeros(Float64,size(datapdf.x)[1])
    init_step = findall(x->x>=0, xar)[1]
    @inbounds @views for ix=init_step:size(xar)[1]-1
        @inbounds @views for jx=1:size(data)[1]
            xar[ix+1]>=abs(data[jx]) && xar[ix]<=abs(data[jx]) && (rep[ix]+=datapdf.density[ix].*abs(data[jx]).^order)
        end
    end
    ddcl = abs(maximum(data)-minimum(data))/size(data)[1] 
    return(sum(rep)*ddcl)
end

function fct_struct!(data,datapdf, order,rep,xar)
    init_step = findall(x->x>=0, xar)[1]
    @inbounds @views while init_step<size(xar)[1]-1
        xar[init_step+1].>=abs(data[:]) && xar[init_step].<=abs(data[:]) && (rep[init_step]+=datapdf.density[init_step].*abs(data[:]).^order)
        init_step+=1
    end
    ddcl = abs(maximum(data)-minimum(data))/size(data)[1] 
    return(sum(rep)*ddcl)
end




"""
    fit_fctstruct(xdata,ydata,y0,confidinterv ; confinterv=true)

Fit the model ydata=A*(xdata)^B. Can return the confidence interval too of the model. Return in first B, then A.
"""
function fit_fctstruct(xdata,ydata,confidinterv ; confinterv=true)
    #model(x,xhi) = xhi[2].*x.^xhi[1]
    fitted = curve_fit(PowerFit, ydata, xdata)
    #confid_interval = confidence_interval(fit, confidinterv)
    confinterv==false && return(fitted.coefs[2],fitted.coefs[1])
    return(fitted.coefs[2],fitted.coefs[1],confid_interval)
end

#function fit_fctstruct(xdata,ydata,y0,confidinterv ; confinterv=true)
#    model(x,xhi) = xhi[2].*x.^xhi[1]
#    fit = curve_fit(model, ydata, xdata, y0)
#    #confid_interval = confidence_interval(fit, confidinterv)
#    confinterv==false && return(fit.param[1],fit.param[2])
#    return(fit.param[1],fit.param[2],confid_interval)
#end



"""
    pdf_normed(data,bin,bound::Tuple)

Compute the PDF of the data given, and normalize it by centered at 0 mean, a standard deviation of unity and an area=1.
"""
function pdf_normed(data,bin,bound::Tuple,blank)
    #data     = Data_preparation.delete_allnotvalue(data,blank)
    #hist     = fit(Histogram,data[:])

    dens_dat = kde(data[:],npoints=bin)#,boundary=bound)
    mu       = moment(dens_dat.x,1,aweights(dens_dat.density),0)
    sig      = sqrt(moment(dens_dat.x,2,aweights(dens_dat.density)))

    dens_dat.density = Data_preparation.blank_inf(dens_dat.density,0.0,blank)
    #dens_dat.density = Data_preparation.delete_allnotvalue(dens_dat.density,blank)

    dens_dat.x       = (dens_dat.x.-mu)./sig
    DELTAX           = abs(dens_dat.x[2]-dens_dat.x[1])
    FACT             = sum(dens_dat.density)*DELTAX
    dens_dat.density = dens_dat.density./FACT  # Produce area=1, normalizing the histogram
    #dens_dat.density = normalize(dens_dat.density.*deltax,2)#./sig
    return dens_dat,sig,mu
end





"""
    xhi_fct_p(pvec,struct_lag,lagvec,y0,confidinterv)   

Fit the model y=A*(x)^B for multiple order and lag of structure functions.
"""
function xhi_fct_p(pvec,struct_lag,lagvec,confidinterv)
    zeta = Array{Float64}(undef,size(pvec)[1],2)
    for ix=1:size(pvec)[1]
        zeta[ix,:] .= fit_fctstruct(struct_lag[ix,:],struct_lag[3,:],confidinterv,confinterv=false)
    end
    return(zeta)
end

function xhi_fct_p!(pvec,struct_lag,lagvec,y0,confidinterv,zeta)
    for ix in eachindex(pvec)
        zeta[ix,:] .= fit_fctstruct(struct_lag[ix,:],struct_lag[3,:],confidinterv,confinterv=false)
    end
    return(zeta)
end



end #module
